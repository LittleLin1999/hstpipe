"""TweakReg candidate fits are accepted only after measured quality checks."""
from pathlib import Path
import json,shutil,tempfile,os
import numpy as np
from astropy.io import fits

def align(inputs,options=None):
    """Low-level DrizzlePac adapter; callers own acceptance of candidate files."""
    from drizzlepac import tweakreg
    options=dict(options or {})
    if 'updatehdr' in options and not options['updatehdr']:raise ValueError('Alignment must update working-copy WCS')
    options.update(updatehdr=True,interactive=False)
    tweakreg.TweakReg([str(p) for p in inputs],**options)

def acceptable_fit(nmatch,xrms,yrms,dx,dy,options):
    return (nmatch>=options['minobj'] and np.isfinite([xrms,yrms,dx,dy]).all()
            and max(abs(xrms),abs(yrms))<=options['max_rms']
            and np.hypot(dx,dy)<=options['max_shift'])

def register(inputs,directory,diagnostics,options,*,paths=None):
    from photutils.detection import DAOStarFinder
    from astropy.stats import sigma_clipped_stats
    from .qa import image_plot
    report=dict(status='original_wcs',reference=str(inputs[0]),fits=[],warnings=[])
    mode=options['mode']
    if mode=='original' or len(inputs)==1:
        report['warnings'].append('Relative alignment not measured; inspect calibrated WCS')
        return report
    counts=[];catalogs={}
    from .astrometry.matching import detect_sources
    for p in inputs:
        cat=detect_sources(p,dict(options,snr=5,fwhm_pixels=2.,fwhm_arcsec=None),paths=paths,label=p.stem+'_legacy_points')
        counts.append(len(cat));catalogs[p.name]=[]
        with fits.open(p) as hd:
            for h in hd:
                if h.name!='SCI':continue
                v=h.header.get('EXTVER',1);c=cat[cat['extver']==v]
                catalog=paths.exposure_file('intermediate','alignment',hd[0].header.get('ROOTNAME',p.stem.removesuffix('_flc')),f'sci{v}_sources','cat') if paths else directory/f'{p.stem}_sci{v}_sources.cat'
                np.savetxt(catalog,np.column_stack([c['x']+1,c['y']+1,c['flux']]),fmt='%.8f')
                catalogs[p.name].append(catalog)
    report['detected_sources']=counts
    if min(counts)<options['minobj'] and mode!='user':
        report.update(status='sparse_original_wcs');report['warnings'].append('Too few PSF-like sources; no forced TweakReg solution')
        if mode in ('relative','external'):raise RuntimeError('Requested alignment cannot be constrained: insufficient PSF-like sources')
        return report
    previous=Path.cwd()
    with paths.scratch('candidate-') as tmp:
        work=Path(tmp);candidates=[]
        for p in inputs:
            dest=work/p.name;shutil.copy2(p,dest);candidates.append(dest)
        try:
            os.chdir(work)
            if mode=='user':
                from drizzlepac import updatehdr
                from stwcs.wcsutil import HSTWCS
                from .inventory import validate_flc
                first_version=validate_flc(candidates[0])['imsets'][0]
                reference=HSTWCS(str(candidates[0]),ext=('SCI',first_version))
                transforms=options.get('transformations',{})
                for p in candidates:
                    name=validate_flc(p)['ROOTNAME'];t=transforms.get(name)
                    if t is None:raise ValueError(f'Missing user transformation for {name}')
                    dx,dy=float(t['dx']),float(t['dy'])
                    if not np.isfinite([dx,dy]).all() or np.hypot(dx,dy)>options['max_shift']:raise ValueError('User shift exceeds acceptance limit')
                    updatehdr.updatewcs_with_shift(str(p),reference,wcsname='HSTPIPE',xsh=dx,ysh=dy,
                        rot=float(t.get('rotation',0)),scale=float(t.get('scale',1)),force=True)
                    report['fits'].append(dict(exposure=name,dx=dx,dy=dy,status='user_supplied_unverified'))
                for p,c in zip(inputs,candidates):
                    paths.copy(c,p) if paths else shutil.copy2(c,p)
                report['status']='user_supplied_unverified'
                report['warnings'].append('User transformations applied by STScI updatehdr; relative residuals unmeasured')
                return report
            catfile=work/'catalogs.txt'
            catfile.write_text('\n'.join(str(p)+' '+ ' '.join(str(c) for c in catalogs[p.name]) for p in candidates))
            kwargs=dict(catfile=str(catfile),refimage=str(candidates[0]),updatewcs=False,updatehdr=True,interactive=False,clean=False,writecat=True,
                shiftfile=True,outshifts='shifts.txt',outwcs='reference_wcs.fits',
                minobj=options['minobj'],fitgeometry=options['fitgeometry'],
                residplot='No plot',see2dplot=False,runfile='tweakreg.log')
            if mode=='external':
                if not options.get('refcat'):raise ValueError('external alignment requires refcat (RA/Dec degrees)')
                kwargs['refcat']=str(Path(options['refcat']).resolve());report['reference']=kwargs['refcat']
            align(candidates,kwargs)
            rows=[]
            for line in (work/'shifts.txt').read_text().splitlines():
                if line.strip() and not line.startswith('#'):
                    fields=line.split();rows.append((Path(fields[0]).name,[float(x) for x in fields[1:7]]))
            accepted=len(rows)==len(inputs)
            for name,values in rows:
                dx,dy,rot,scale,xrms,yrms=values
                catalogs=sorted(set(work.glob(Path(name).stem+'*fitmatch*'))|set(work.glob(Path(name).stem+'*_catalog_fit.match')))
                nmatch=0
                for cat in catalogs:
                    nmatch+=sum(bool(s.strip()) and not s.startswith('#') for s in cat.read_text().splitlines())
                reference=name==inputs[0].name and mode!='external'
                ok=reference or acceptable_fit(nmatch,xrms,yrms,dx,dy,options)
                accepted&=ok
                report['fits'].append(dict(exposure=name,nmatch=nmatch,dx=dx,dy=dy,rotation=rot,scale=scale,
                                           xrms=xrms,yrms=yrms,accepted=bool(ok)))
            if accepted:
                for p,c in zip(inputs,candidates):
                    paths.copy(c,p) if paths else shutil.copy2(c,p)
                report['status']='accepted'
            else:
                report['status']='rejected_original_wcs';report['warnings'].append('Unreliable or unmeasured fit; all original WCS retained')
                if mode in ('relative','external'):raise RuntimeError('Requested alignment fit failed acceptance limits')
        except Exception as exc:
            report.update(status='failed_original_wcs');report['warnings'].append(str(exc))
            if mode in ('relative','external','user'):raise RuntimeError('Requested alignment failed; no drizzle attempted') from exc
        finally:
            os.chdir(previous)
            for f in work.iterdir():
                if f.is_file() and f.suffix!='.fits':paths.copy(f,paths.file('intermediate','alignment',f.stem,f.suffix.lstrip('.'))) if paths else shutil.copy2(f,directory/f.name)
                elif f.name=='reference_wcs.fits':paths.copy(f,paths.file('intermediate','alignment',f.stem,f.suffix.lstrip('.'))) if paths else shutil.copy2(f,directory/f.name)
    if report['fits']:
        image_plot([np.array([[r['dx'],r['dy'],r['xrms'],r['yrms']] for r in report['fits']])],
                   ['dx, dy, x RMS, y RMS (reference pixels)'],paths.file('diagnostics','alignment','alignment','png') if paths else diagnostics/'fits.png')
    return report
