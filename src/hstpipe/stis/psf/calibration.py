"""Recalibration of working copies; archive products are never modified."""
from pathlib import Path
import os,shutil,logging
from astropy.io import fits
from .common import *


def recalibrate(c,limit=None):
    from stistools import calstis,x1d,x2d
    base=Path('data/crds').resolve();os.environ['CRDS_PATH']=str(base);os.environ['CRDS_SERVER_URL']='https://hst-crds.stsci.edu'
    import crds
    try:context=c.get('crds_context') or crds.get_default_context('hst')
    except Exception as exc:
        context='ARCHIVE_HEADER_REFERENCES';logging.warning('Current CRDS context unavailable: %s',exc)
    rows=load_json('data/catalogs/observations.json');records=[]
    # Raw-only CALSPEC datasets first; then one monitor provides a baseline check.
    work=[r for r in rows if r.get('raw_header') and r.get('paths',{}).get('raw')]
    work += [r for r in rows if r['program']==18169 and r.get('paths',{}).get('raw')][:1]
    if limit:work=work[:limit]
    for r in work:
        raw=Path(r['paths']['raw']);dest=Path('data/recalibrated')/r['root'];dest.mkdir(parents=True,exist_ok=True);copy=dest/raw.name
        rec=dict(root=r['root'],raw=str(raw),context=context,raw_sha256=sha(raw))
        try:
            with fits.open(raw) as h:
                header=dict(h[0].header);header.update({k:h[1].header[k] for k in ['DATE-OBS','TIME-OBS','EXPSTART','EXPTIME'] if k in h[1].header});times=[hd.header.get('EXPTIME',0) for hd in h if hd.name=='SCI']
            if not any(t>0 for t in times):
                r['status']='REJECT_ZERO_EXPTIME';r['reason']='Actual raw SCI EXPTIME is zero';rec['status']='ZERO_EXPTIME';records.append(rec);continue
            shutil.copyfile(raw,copy);copy.chmod(0o644)
            try:
                if context=='ARCHIVE_HEADER_REFERENCES':raise RuntimeError('No current CRDS context')
                refs=crds.getreferences(header,observatory='hst',context=context)
            except Exception as exc:
                import urllib.request,time
                rec['crds_limitation']=str(exc);rec['reference_selection']='ARCHIVED_HEADER_FALLBACK'
                refs={}
                refdir=Path('data/reference_files').resolve();refdir.mkdir(exist_ok=True)
                for key,value in header.items():
                    if not isinstance(value,str) or not value.startswith('oref$'):continue
                    name=value.split('$')[-1];existing=Path(c['reference_dir'])/name;target=refdir/name
                    if existing.is_file():refs[key.lower()]=str(existing);continue
                    if not target.is_file():
                        for attempt in range(3):
                            try:
                                urllib.request.urlretrieve('https://hst-crds.stsci.edu/unchecked_get/references/hst/'+name,str(target)+'.part')
                                with fits.open(str(target)+'.part') as hh:hh.verify('exception')
                                Path(str(target)+'.part').replace(target);target.chmod(0o444);break
                            except Exception as error:logging.warning('Reference %s: %s',name,error);time.sleep(1+attempt)
                    if target.is_file():refs[key.lower()]=str(target)
            rec['reference_selection']=rec.get('reference_selection','CURRENT_CRDS')
            rec['references']={k:dict(path=v,sha256=sha(v)) for k,v in refs.items() if Path(v).is_file()}
            # Resolve reference names to absolute paths and explicitly request first-order products.
            with fits.open(copy,mode='update') as h:
                for key,value in refs.items():
                    if Path(value).is_file():h[0].header[key.upper()]='oref$'+Path(value).name
                for key in ['X1DCORR','X2DCORR','DISPCORR','FLUXCORR','HELCORR','BACKCORR']:h[0].header[key]='PERFORM'
            # Other archived tables not selected by CRDS can be found in the user's existing reference cache.
            refview=Path('data/reference_view').resolve();refview.mkdir(exist_ok=True)
            for value in refs.values():
                if Path(value).is_file():
                    link=refview/Path(value).name
                    if not link.exists():link.symlink_to(Path(value).resolve())
            for value in header.values():
                if isinstance(value,str) and value.startswith('oref$'):
                    name=value.split('$')[-1];src=Path(c['reference_dir'])/name;link=refview/name
                    if src.is_file() and not link.exists():link.symlink_to(src.resolve())
            os.environ['oref']=str(refview)+'/'
            trailer=str(dest/'calstis.log')
            prefix=str(dest/r['root'])
            if not Path(prefix+'_flt.fits').exists():
                code=calstis.calstis(str(copy),outroot=str(dest)+'/',wavecal=str(Path(r['paths']['wav']).resolve()) if r['paths'].get('wav') else '',trailer=trailer)
                rec['calstis_returncode']=code
                if code:raise RuntimeError(f'calstis returned {code}; see {trailer}')
            outputs={k:prefix+'_'+k+'.fits' for k in ['flt','x1d','x2d'] if Path(prefix+'_'+k+'.fits').exists()}
            if 'flt' in outputs and 'x1d' not in outputs:
                x1d.x1d(outputs['flt'],output=prefix+'_x1d.fits',trailer=str(dest/'x1d.log'))
                if Path(prefix+'_x1d.fits').exists():outputs['x1d']=prefix+'_x1d.fits'
            if 'flt' in outputs and 'x2d' not in outputs:
                x2d.x2d(outputs['flt'],output=prefix+'_x2d.fits',trailer=str(dest/'x2d.log'))
                if Path(prefix+'_x2d.fits').exists():outputs['x2d']=prefix+'_x2d.fits'
            rec['outputs']=outputs;rec['status']='COMPLETE' if all(k in outputs for k in ['flt','x1d','x2d']) else 'PARTIAL'
            r['recalibrated_paths']=outputs
            if 'flt' not in r['paths'] and 'x1d' in outputs:r['paths'].update(outputs);r['status']='RECALIBRATED_FOR_ANALYSIS';r['aperture']=header['APERTURE']
        except Exception as e:rec['status']='FAILED';rec['error']=str(e);logging.exception('Recalibration %s',r['root'])
        records.append(rec)
    save_json('data/catalogs/recalibration.json',records);save_json('data/catalogs/observations.json',rows)
    return records
