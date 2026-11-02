"""Reproducible UVIS stages, with immutable upstream stage checkpoints."""
import json,shutil,os,tempfile
from pathlib import Path
from astropy.io import fits
import numpy as np
from ..common.runner import run
from ..common.cache import StageCache
from ..common.provenance import identity,versions
from .config import settings
from .qa import inspect_flc,validate_drizzle,cr_report,final_report
from .background import prepare_background,restore_background
from .storage import stage_keys,registration_table
from .alignment import register
from .drizzle import drizzle,annotate_product,apply_cr_masks

def write_json(path,value):
    # STScI reference-image fit residuals are undefined (NaN), not zero.
    def finite(obj):
        if isinstance(obj,(float,np.floating)) and not np.isfinite(obj):return None
        if isinstance(obj,dict):return {k:finite(v) for k,v in obj.items()}
        if isinstance(obj,(list,tuple)):return [finite(v) for v in obj]
        return obj
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(finite(value),indent=2,default=str,allow_nan=False))

def _engine(config,paths,exposures,*,prepare_only=False):
    opts=settings(config.options);two_pass=opts['astrometry']['method']=='two_pass';report=dict(inventory={e.path.name:e.metadata for e in exposures},grouping=[],excluded=[],exposures={},background={},alignment={},cr={})
    from .dq import policy,mdriztab_recommendation
    report['dq_policy']=policy(opts)
    report['mdriztab']=mdriztab_recommendation(exposures[0].path,len(exposures)) if exposures else {}
    rows=[]
    for e in exposures:
        report['grouping'].append(dict(path=str(e.path),reason=e.metadata['grouping_reason'],geometry=e.metadata['geometry']))
        if e.path.name in opts['exclude'] or e.metadata['ROOTNAME'] in opts['exclude']:
            report['excluded'].append(dict(path=str(e.path),reason='explicit user exclusion'));continue
        if any(e.metadata[k]!='COMPLETE' for k in ('BIASCORR','BLEVCORR','DARKCORR','FLATCORR','PCTECORR')) and not opts['allow_incomplete_calibration']:
            raise ValueError(f'Incomplete calibration: {e.path.name}; inspect calibration switches')
        if any(not g['has_distortion'] for g in e.metadata['geometry']):
            raise ValueError(f'Missing embedded distortion WCS: {e.path.name}; update working-copy WCS with STScI references before reduction')
        rows.append(e)
    if not rows:raise ValueError('All group exposures explicitly excluded')
    report['dq_policy']=policy(opts,len(rows))
    report['sampling']=dict(final_scale=opts['drizzle']['final_scale'],final_pixfrac=opts['drizzle']['final_pixfrac'],driz_sep_pixfrac=opts['drizzle']['driz_sep_pixfrac'],nimages=len(rows),commanded_offsets={e.metadata['ROOTNAME']:[e.metadata.get('POSTARG1'),e.metadata.get('POSTARG2')] for e in rows},rationale='Conservative default 0.04 arcsec/pixel and pixfrac 1; 0.0333 is configurable after evaluating dither sampling/coverage, not universally superior')
    background=paths.area('intermediate','background',create=True)
    alignment=paths.area('intermediate','alignment')
    work=paths.area('intermediate','drizzle')
    diagnostics=paths.area('diagnostics','drizzle')
    cache=paths.cache()
    copies=[paths.exposure_file('intermediate','alignment',e.metadata['ROOTNAME'],'flc') for e in rows]
    skyfiles=[]
    for e in rows:
        for kind in ('sky','mask'):
            skyfiles.append(paths.exposure_file('intermediate','background',e.metadata['ROOTNAME'],kind))
        amp=paths.exposure_file('intermediate','background',e.metadata['ROOTNAME'],'amp')
        if amp.exists():skyfiles.append(amp)
    background_key,alignment_key,_=stage_keys(cache,rows,opts,copies)
    bg_summary=paths.file('intermediate','background','background')
    alignment_summary=paths.file('intermediate','alignment','alignment')
    registration=paths.file('intermediate','alignment','registration','ecsv')
    bg_valid=config.resume and cache.valid('uvis_background',background_key)
    if bg_valid:
        report.update(json.loads(bg_summary.read_text()))
    else:
        for e,p in zip(rows,copies):
            report['exposures'][e.path.name]=inspect_flc(e.path,paths.area('diagnostics','alignment'),paths=paths)
            paths.detach(p)
            paths.detach(paths.exposure_file('intermediate','background',e.metadata['ROOTNAME'],'amp'))
            report['background'][e.path.name]=prepare_background(e.path,p,background,opts['background'],paths=paths)
        write_json(bg_summary,{k:report[k] for k in ('exposures','background')})
        skyfiles=[p for p in paths._files if p.parent==background and p.suffix=='.fits' and p.exists()]
        cache.record('uvis_background',background_key,skyfiles+[bg_summary])
    if bg_valid and not all(p.exists() for p in copies):
        for e,p in zip(rows,copies):restore_background(e.path,p,paths)
    sky_rates=[(r['measured']-r['applied'])/e.metadata['EXPTIME'] for e in rows for r in report['background'][e.path.name] if r.get('reliable')]
    report['background_discontinuity_range']=float(np.ptp(sky_rates)) if sky_rates else None
    report['background_discontinuity_caveat']='Reliable source-masked chip estimates only; may retain astrophysical extended emission'
    from .astrometry.identification import identify
    identification=identify([e.path for e in rows],paths,opts['astrometry'])
    report['source_identification']={k:v for k,v in identification.items() if k!='catalogs'}
    _,alignment_key,_=stage_keys(cache,rows,opts,copies,identification)
    if config.resume and cache.valid('uvis_alignment',alignment_key):
        report['alignment']=json.loads(alignment_summary.read_text())
    else:
        if bg_valid:
            for e,p in zip(rows,copies):
                paths.detach(p)
                restore_background(e.path,p,paths)
        if two_pass:
            from .astrometry.relative import inspect_wcs
            report['input_wcs_solutions']=inspect_wcs(copies)
        relative_options=dict(opts['alignment'])
        if two_pass:relative_options['mode']=opts['astrometry']['relative']['mode'] if opts['astrometry']['relative']['refine_if_needed'] else 'original'
        if opts['alignment']['mode']=='external':raise ValueError('External catalogs cannot initialize FLC-relative registration; use the absolute stage after provisional drizzle')
        if opts['alignment']['mode']=='user':
            report['requested_alignment']=register(copies,alignment,paths.area('diagnostics','alignment'),opts['alignment'],paths=paths)
        from .astrometry.relative_engine import refine_relative
        report['alignment']=refine_relative(copies,paths,opts['astrometry']['relative'],identification=identification)
        if 'requested_alignment' in report:report['alignment']['requested_initialization']=report['requested_alignment']
        write_json(alignment_summary,report['alignment'])
        registration_table(registration,report['alignment'])
        catalog_files=[Path(c['path']) for c in report['alignment'].get('detector_catalogs',[])]
        quality_files=[p for p in paths._files if p.parent==alignment and (p.name.endswith('_quality.ecsv') or p.name.endswith('_rejected_cutouts.fits')) and p.exists()]
        cache.record('uvis_alignment',alignment_key,list(dict.fromkeys(copies+[alignment_summary,registration]+catalog_files+quality_files)))
    from .chips import report_summary
    report['multi_chip']=report_summary(copies,rows,report['alignment'])
    report['stages']=dict(source_detection=identification.get('full_field',dict(status='fallback')),source_identification=dict(status=identification['status'],external_positions_used_for_relative_fit=False),relative_alignment=dict(status=report['alignment']['status'],cache_key=alignment_key,inputs=[str(e.path) for e in rows],outputs=[str(p) for p in copies]+[str(alignment_summary)]))
    if not report['alignment'].get('may_combine',False):
        report.update(status='relative_astrometry_unvalidated',relative_astrometry=report['alignment'],astrometry=dict(status='unvalidated',reason='Relative registration failed; no absolute calibration or TweakBack'))
        for stage in ('provisional_drizzle','cross_filter_alignment','absolute_alignment','tweakback','final_drizzle','final_validation'):report['stages'][stage]=dict(status='blocked',reason='Relative registration unvalidated')
        write_json(paths.file('diagnostics','drizzle','qa_summary'),report)
        return []
    masks=[paths.exposure_file('intermediate','drizzle',e.metadata['ROOTNAME'],'crmask') for e in rows]
    cr_options={k:v for k,v in opts['drizzle'].items() if not k.startswith('final_')}
    _,_,cr_key=stage_keys(cache,rows,opts,copies)
    cr_summary=paths.file('intermediate','drizzle','cr')
    provisional=paths.file('intermediate','drizzle','provisional','fits')
    provisional_inputs=paths.file('intermediate','drizzle','provisional_inputs')
    if two_pass:
        from .astrometry.relative import validate_relative,inspect_wcs
        from .astrometry.workflow import finalize
        report['input_wcs_solutions']=inspect_wcs(copies)
        relative=report['alignment']
        if not relative['may_combine']:
            report.update(status='astrometry_unvalidated',relative_astrometry=relative,
                          astrometry=dict(status='unvalidated',reason='Mandatory relative registration could not be established'))
            write_json(paths.file('diagnostics','drizzle','qa_summary'),report)
            return []
        provisional_key=cache.key(copies+[Path(__file__).parent/'drizzle.py'],dict(cr=cr_options,provisional=opts['astrometry']['provisional'],versions=versions()))
        def finish_two_pass():
            expected=json.loads(provisional_inputs.read_text())
            if [identity(p) for p in copies]!=expected:raise RuntimeError('Relative FLC WCS changed since provisional drizzle')
            report['relative_astrometry']=relative
            report['provisional_image']=dict(path=str(provisional),relative_status=relative['status'])
            from .chips import coverage_report
            report['provisional_chips']=coverage_report(provisional,copies)
            report['stages']['provisional_drizzle']=dict(status='complete',cache_key=provisional_key,inputs=expected,outputs=[str(provisional)])
            if prepare_only:
                report.update(status='provisional_ready',astrometry=dict(status='pending',reason='Internal provisional preparation; identity priors are separate from absolute fitting'))
                write_json(paths.file('diagnostics','drizzle','qa_summary'),report)
                cache.record('uvis_provisional',provisional_key,[provisional,provisional_inputs])
                return []
            cross=paths.file('intermediate','alignment','cross_filter')
            if cross.exists():report['cross_filter_alignment']=json.loads(cross.read_text())
            products,assessment=finalize(provisional,copies,rows,paths,opts['astrometry'],opts['drizzle'],relative,opts['target_cr_fraction'])
            report.update(assessment);report['effective_parameters']=opts
            report['stages']['cross_filter_alignment']=report.get('cross_filter_alignment',dict(status='not_applicable'))
            report['stages']['absolute_alignment']=dict(status=assessment.get('astrometry',{}).get('validation',{}).get('status','unvalidated'),outputs=[str(paths.file('intermediate','alignment','absolute'))])
            report['stages']['tweakback']=dict(status=assessment.get('tweakback',{}).get('status','blocked'),outputs=[h['path'] for h in assessment.get('tweakback',{}).get('headerlets',[])])
            report['stages']['final_drizzle']=dict(status='complete' if products else 'blocked',outputs=[str(p) for p in products])
            report['stages']['final_validation']=assessment.get('astrometry',{}).get('final_validation',dict(status='blocked'))
            write_json(paths.file('diagnostics','drizzle','qa_summary'),report)
            cache.record('uvis_provisional',provisional_key,[provisional,provisional_inputs])
            return products
        if config.resume and cache.valid('uvis_provisional',provisional_key):return finish_two_pass()
    previous=Path.cwd()
    with paths.scratch('run-') as tmp:
        scratch=Path(tmp);inputs=[]
        for p in copies:
            dest=scratch/p.name;shutil.copy2(p,dest);inputs.append(dest)
        try:
            os.chdir(scratch)
            if config.resume and cache.valid('uvis_cr',cr_key):
                report['cr']=json.loads(cr_summary.read_text())
                from drizzlepac.resetbits import reset_dq_bits
                for p,mask in zip(inputs,masks):
                    reset_dq_bits(str(p),4096,extname='DQ')
                    with fits.open(p,mode='update') as hd,fits.open(mask) as mh:
                        for h in mh:
                            if h.name=='SCI':hd['DQ',h.header.get('EXTVER',1)].data[h.data!=0]|=4096
                        hd[0].header['HPCRRUN']=True
            else:
                if len(inputs)>1:
                    drizzle(inputs,scratch/'cr',dict(cr_options,driz_combine=False,runfile='cr_astrodrizzle.log'))
                    apply_cr_masks(inputs)
                for e,p,mask in zip(rows,inputs,masks):
                    report['cr'][e.path.name]=cr_report(e.path,p,mask,diagnostics,opts['target_cr_fraction'],paths=paths)
                write_json(cr_summary,report['cr']);cache.record('uvis_cr',cr_key,masks+[cr_summary])
            prefix=scratch/paths.group_id
            final_options=dict(opts['drizzle'],static=True,driz_separate=False,median=False,blot=False,driz_cr=False,driz_combine=True,runfile='final_astrodrizzle.log')
            if two_pass:final_options.update(opts['astrometry']['provisional'])
            drizzle(inputs,prefix,final_options)
            product=prefix.with_name(prefix.name+'_drc.fits')
            if not product.exists():raise FileNotFoundError(f'AstroDrizzle did not create {product}')
            annotate_product(product,opts['drizzle']['final_wht_type'],len(inputs),report['alignment']['status'],opts['drizzle']['final_bits'])
            validate_drizzle(product)
            final=provisional if two_pass else paths.product('drc');final.parent.mkdir(parents=True,exist_ok=True);paths.copy(product,final)
            if two_pass:
                with fits.open(provisional,mode='update') as hd:hd[0].header['HPPASS']='PROVISIONAL'
                write_json(provisional_inputs,[identity(p) for p in copies])
                for log in scratch.glob('*.log'):paths.copy(log,paths.file('logs',None,'provisional_'+log.stem,'log'))
                return finish_two_pass()
            report['final']=final_report(final,copies,diagnostics,paths=paths)
            if len(rows)==1:report['final']['limitations'].append('Single exposure: no reliable multi-exposure CR rejection; NOT fully CR cleaned')
            residuals=[]
            for e in rows:
                for sky in report['background'][e.path.name]:residuals.append((sky['measured']-sky['applied'])/e.metadata['EXPTIME'])
            report['background_residual_rates']=residuals
            report['background_discontinuity_range']=float(np.ptp(residuals)) if residuals else None
            report['status']='review_required'
            report['effective_parameters']=opts
            write_json(paths.file('diagnostics','drizzle','qa_summary'),report)
            if opts['keep_drizzle_intermediates']:
                for p in scratch.iterdir():
                    if p.is_file() and p not in inputs and p!=product:paths.copy(p,paths.file('intermediate','drizzle',p.stem,p.suffix.lstrip('.')))
            # Always retain numerical logs, not duplicate calibrated exposures.
            for p in scratch.glob('*.log'):paths.copy(p,paths.file('logs',None,p.stem,'log'))
            return [final]
        finally:os.chdir(previous)

def reduce_uvis(config,*,dry_run=False):
    from ..common.config import ReductionConfig
    if isinstance(config,(str,Path)):config=ReductionConfig.read(config)
    elif isinstance(config,dict):config=ReductionConfig(**config)
    if config.instrument!='uvis':raise ValueError('reduce_uvis requires instrument=uvis')
    config.options=settings(config.options)
    from .astrometry.catalogs import _FAILED_REQUESTS
    _FAILED_REQUESTS.clear()
    if dry_run or config.options['astrometry']['method']=='relative_only':return run(config,_engine,dry_run)
    from .astrometry.orchestration import relative_then_absolute
    return relative_then_absolute(config,_engine)
