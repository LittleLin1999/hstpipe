"""Thin orchestration of the existing, pinned hstpipe STIS algorithms.

The existing hstpipe installation is used; no source is patched. Archive FLTs/X2Ds are symlinked into compatible
stellar input groups. Native geometry/background/flux/drizzle/extraction come
from hstpipe; this module only configures, measures registration, and catalogs.
"""
from pathlib import Path
from collections import defaultdict
import contextlib,logging,sys,os,time,urllib.request,json
import numpy as np
from astropy.io import fits
from astropy.table import Table
from .common import *


def use_hstpipe():
    from hstpipe.stis import reduce_stis
    return reduce_stis


def reference_view(r,c):
    target=Path('data/hstpipe_references').resolve();target.mkdir(parents=True,exist_ok=True)
    with fits.open(r['paths']['flt']) as h:header=h[0].header
    keys=['PHOTTAB','TDSTAB','APDESTAB','SDCTAB','SPTRCTAB','DISPTAB','INANGTAB','PCTAB','XTRACTAB','APERTAB']
    for key in keys:
        name=str(header.get(key,'' )).split('$')[-1]
        if not name or name in ['N/A','NONE']:continue
        name=Path(name).name;link=target/name
        if link.exists():continue
        sources=[Path(c['reference_dir'])/name,Path('data/reference_files')/name,*Path('data/crds/references').rglob(name)]
        src=next((p for p in sources if p.is_file()),None)
        if src:link.symlink_to(src.resolve());continue
        for attempt in range(3):
            try:
                tmp=target/(name+'.part');url='https://hst-crds.stsci.edu/unchecked_get/references/hst/'+name
                urllib.request.urlretrieve(url,tmp)
                with fits.open(tmp) as hh:hh.verify('exception')
                tmp.replace(link);link.chmod(0o444);break
            except Exception as e:logging.warning('Reference %s attempt %s: %s',name,attempt+1,e);time.sleep(attempt+1)
        if not link.exists():raise FileNotFoundError(key+': '+name)
    return target


def reduce_stars(c,only_roots=None):
    from .extraction import refresh_core_quality
    refresh_core_quality(c)
    reduce=use_hstpipe()
    from hstpipe.common.config import ReductionConfig
    from hstpipe.stis.config import Config
    from hstpipe.stis.io import read_inputs,validate
    from hstpipe.stis.reference import prepare,operator
    rows=catalog_apertures(c);groups=defaultdict(list);out=load_json('data/catalogs/hstpipe_reductions.json',[]);previous={r['group']:r for r in out}
    for r in rows:
        if only_roots and r['root'] not in only_roots:continue
        if r.get('psf_status')!='CORE_USABLE' or not r.get('point_source_prior'):continue
        if not all(k in r.get('paths',{}) for k in ['flt','x2d']):continue
        with fits.open(r['paths']['flt']) as h:
            hdr=h[0].header;pa=float(h[1].header.get('PA_APER',0))
            if hdr.get('WAVECORR')!='COMPLETE':
                previous[r['root']]=dict(group=r['root'],roots=[r['root']],target=r['target'],status='UNSUPPORTED_WAVECORR',reason='Existing hstpipe requires WAVECORR=COMPLETE; nominal-dispersion raw-only calibrations are native-only comparisons')
                continue
        # Same star, exact configuration, visit, PA, and coarse detector position.
        key=(r['target'],r.get('configuration_aperture',r['aperture']),r['program'],r['root'][:6],round(pa,2),int(np.rint(r['trace_y']/c.get('detector_position_group_pixels',100))))
        groups[key].append(r)
    for key,rr in sorted(groups.items(),key=lambda kv:([18169,17644,17649,18229,18113,16966,17804,17223].index(kv[0][2]) if kv[0][2] in [18169,17644,17649,18229,18113,16966,17804,17223] else 99,kv[0])):
        setup=safe(f'{key[1]}_p{key[2]}_{key[3]}_y{key[5]}');group=safe(key[0]+'_'+setup);rec=dict(group=group,roots=[r['root'] for r in rr],target=key[0],aperture=key[1],program=key[2],trace_y_min=min(r['trace_y'] for r in rr),trace_y_max=max(r['trace_y'] for r in rr))
        log=Path('logs')/(group+'_hstpipe.log');view=Path('data/reduction_inputs')/group;view.mkdir(parents=True,exist_ok=True)
        try:
            refs=None
            for r in rr:
                refs=reference_view(r,c)
                for suffix in ['flt','x2d']:
                    src=Path(r['paths'][suffix]).resolve();link=view/src.name
                    if not link.exists():link.symlink_to(src)
            controls=Path('data/controls').resolve();controls.mkdir(exist_ok=True)
            directory=Path('redux').resolve()/key[0]/'stis'/'G140L'/setup
            options=dict(c.get('reduction_options',{}));options.update(reference_dir=str(refs),prepare_controls=True,control_dir=str(controls))
            # Use unchanged hstpipe controls and operator to measure source offsets
            # in its own target-relative geometry. No FLT/WCS header is edited.
            regfile=Path('data/catalogs')/(group+'_registration.json')
            registration=load_json(regfile,{})
            with log.open('a') as stream,contextlib.redirect_stdout(stream),contextlib.redirect_stderr(stream):
                native=Config(project_dir=str(Path('.').resolve()),input_glob=str(view.resolve()/'*_flt.fits'),mast_dir=str(view.resolve()),output_dir=str(directory),**options)
                exp=read_inputs(sorted(view.glob('*_flt.fits')));validate(exp,native);prepare(exp,native)
                offsets={}
                for e in exp:
                    if e.path.name not in registration:
                        op=operator(e,native);root=e.path.name[:9];trace=np.load(f'products/individual_profiles/{root}_trace.npz')['measured_trace'];cols=np.arange(e.counts.shape[1]);yy=np.clip(np.rint(trace).astype(int),0,e.counts.shape[0]-1)
                        wave=op['observed_wavelength'][yy,cols];mask=(wave>1330)&(wave<1680)&op['valid'][yy,cols]
                        # Subpixel spatial correction from local mapping gradient.
                        sy=op['spatial'][yy,cols]+(trace-yy)*np.gradient(op['spatial'],axis=0)[yy,cols]
                        value=float(np.median(sy[mask]));scatter=float(1.4826*np.median(abs(sy[mask]-value)))
                        registration[e.path.name]=dict(source_offset_arcsec=value,trace_scatter_arcsec=scatter,n_columns=int(mask.sum()),postarg1=e.p.get('POSTARG1'),postarg2=e.p.get('POSTARG2'),shifta1=e.h.get('SHIFTA1'),shifta2=e.h.get('SHIFTA2'),pa_aper=e.h.get('PA_APER'),spectral_method='Archive wavecal plus HSTCAL wavelength mapping; no additional spectral shift applied')
                    offsets[e.path.name]=registration[e.path.name]['source_offset_arcsec']
                save_json(regfile,registration);options['source_offsets']=offsets
                cfg=ReductionConfig(str(view.resolve()),str(Path('redux').resolve()),'stis',key[0],'G140L',setup,options=options,resume=True)
                configfile=Path('config/reductions')/(group+'.yaml');configfile.parent.mkdir(exist_ok=True);configfile.write_text(yaml.safe_dump(cfg.effective(),sort_keys=False))
                products=reduce(cfg)
            rec.update(status='COMPLETE',directory=str(directory),products=[str(p) for p in (directory/'products').glob('*.fits')],configuration=str(configfile),log=str(log),registration=str(regfile))
        except Exception as e:rec.update(status='FAILED',error=str(e),log=str(log));logging.exception('hstpipe %s',group)
        previous[group]=rec;save_json('data/catalogs/hstpipe_reductions.json',list(previous.values()));logging.info('hstpipe %s: %s',group,rec['status'])
    return list(previous.values())
