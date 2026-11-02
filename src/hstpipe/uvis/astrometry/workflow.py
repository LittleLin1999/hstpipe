"""Validated provisional WCS -> official TweakBack -> new detector-frame drizzle."""
import json,tempfile
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from .matching import science_header,science_shape
from ...common.cache import StageCache
from ...common.provenance import identity,versions
from ...common.target_runner import save
from ..drizzle import drizzle,annotate_product,apply_cr_masks
from ..qa import cr_report,final_report,validate_drizzle,image_plot
from .alignment import absolute_alignment,observation_epoch
from .matching import detect_sources
from .absolute_matching import validate_absolute as validate
from .catalogs import gaia_at_epoch,legacy_at_epoch
from .wcs import detector_wcs


def _science_versions():
    from importlib.metadata import version
    return dict(versions(),stwcs=version('stwcs'))


def tweakback_inputs(provisional,relative_inputs,working_inputs,solution,paths,options):
    from drizzlepac import tweakback
    from stwcs.wcsutil import headerlet
    cache=paths.cache()
    provenance=paths.file('intermediate','alignment','tweakback')
    headerlets=[paths.exposure_file('intermediate','alignment',fits.getheader(p)['ROOTNAME'],'absolute_hdrlet') for p in relative_inputs]
    import ast,hashlib,inspect
    code=hashlib.sha256(ast.dump(ast.parse(inspect.getsource(tweakback_inputs))).encode()).hexdigest()
    key=cache.key([provisional]+relative_inputs,dict(solution=solution,options=options,versions=_science_versions(),code=code))
    if cache.valid('uvis_tweakback',key):
        for p,h in zip(working_inputs,headerlets):headerlet.apply_headerlet_as_primary(str(p),str(h),attach=True,archive=True,force=False)
        return dict(json.loads(provenance.read_text()),restored_from_cache=True)
    before=[];chips_after=[]
    for source,p in zip(relative_inputs,working_inputs):
        if identity(source)['sha256']!=identity(p)['sha256']:raise ValueError('TweakBack must begin with unchanged relative-alignment working copies')
        with fits.open(p) as hd:
            if hd[0].header.get('HPTBID'):raise ValueError('Astrometric correction already applied to working copy')
            before.append(dict(exposure=p.name,extensions=[h.header.get('EXTVER',1) for h in hd if h.name=='SCI']))
    if solution['correction_applied']:
        if not options['enabled']:raise ValueError('An absolute correction requires TweakBack before final drizzle')
        if not hasattr(tweakback,'apply_tweak'):raise RuntimeError('Two-pass astrometry requires DrizzlePac >=3.11 with apply_tweak')
        version=science_header(provisional).get('EXTVER',1)
        tweakback.apply_tweak(f'{provisional}[SCI,{version}]',orig_wcs_name=solution['original_wcs_name'],orig_wcs_key=solution['original_wcs_key'],
            tweaked_wcs_name=solution['corrected_wcs_name'],output_wcs_name=solution['corrected_wcs_name'],
            input_files=','.join(str(p) for p in working_inputs),default_extname='SCI')
    for p,hdr in zip(working_inputs,headerlets):
        with fits.open(p) as hd:
            for h in hd:
                if h.name!='SCI':continue
                world=detector_wcs(hd,h.header.get('EXTVER',1)).all_pix2world([[0,0],[h.shape[1]/2,h.shape[0]/2]],0)
                if not np.isfinite(world).all():raise ValueError('Invalid FLC WCS after TweakBack')
                if solution['correction_applied'] and h.header.get('WCSNAME')!=solution['corrected_wcs_name']:raise ValueError('TweakBack did not update every SCI WCS')
        from ..chips import geometry
        chips_after.append(dict(exposure=p.name,chips=geometry(p)))
        hdr.parent.mkdir(parents=True,exist_ok=True)
        hlet=headerlet.create_headerlet(str(p),hdrname=p.stem+'-HPABS',sciext='SCI')
        if hlet is None:raise ValueError('STWCS could not preserve corrected exposure headerlet')
        hlet.tofile(str(hdr),clobber=True)
    record=dict(status='applied' if solution['correction_applied'] else 'not_needed_validated_existing',
        api='drizzlepac.tweakback.apply_tweak',original_wcs=solution['original_wcs_name'],corrected_wcs=solution['corrected_wcs_name'],
        inputs=[identity(p) for p in relative_inputs],chip_wcs_after=chips_after,updated_extensions=before,headerlets=[identity(p) for p in headerlets],applied_once=True)
    save(provenance,record);cache.record('uvis_tweakback',key,headerlets+[provenance])
    return record


def validate_relative_transfer(relative_inputs,working_inputs,relative,options,paths):
    """Recheck fixed internal source pairs after TweakBack, without external fits."""
    from .sources import sources,project
    from .relative_fit import match,catalog_match,quality
    if len(working_inputs)==1:
        return dict(status='single_exposure',pairs=[],external_catalog_used=False)
    entries={r['exposure']:r['path'] for r in relative.get('detector_catalogs',[])}
    catalogs=[Table.read(entries[p.name],format='ascii.ecsv') if p.name in entries else sources(p,options['detection'],True,6) for p in relative_inputs]
    reference=relative.get('reference') or relative.get('reference_selection',{}).get('selected')
    order=sorted(range(len(relative_inputs)),key=lambda i:relative_inputs[i].name!=reference)
    before=[relative_inputs[i] for i in order];after=[working_inputs[i] for i in order];cats=[catalogs[i] for i in order]
    with fits.open(before[0]) as hd:old=WCS(hd['SCI'].header,hd)
    with fits.open(after[0]) as hd:new=WCS(hd['SCI'].header,hd)
    a=[project(c,p,old) for c,p in zip(cats,before)];b=[project(c,p,new) for c,p in zip(cats,after)];pairs=[]
    for i in range(1,len(before)):
        try:
            ri,ii=catalog_match(a[0],a[i],cats[0],cats[i],options,True)
            errors=np.hypot(cats[0]['centroid_error_pix'][ri],cats[i]['centroid_error_pix'][ii])
            initial=quality(a[0][ri],a[i][ii],errors,options,identities_verified='common_source_id' in cats[0].colnames and 'common_source_id' in cats[i].colnames)
            current=quality(b[0][ri],b[i][ii],errors,options,identities_verified='common_source_id' in cats[0].colnames and 'common_source_id' in cats[i].colnames)
            pairs.append(dict(reference=before[0].name,exposure=before[i].name,before=initial,after=current))
        except Exception as exc:pairs.append(dict(exposure=before[i].name,after=dict(status='unvalidated'),error=str(exc)))
    result=dict(status='validated' if all(p['after']['status']=='validated' for p in pairs) else 'unvalidated',pairs=pairs,
        external_catalog_used=False,method='Same detector centroids and fixed internal pairs, projected through updated distortion-aware WCS')
    save(paths.file('diagnostics','alignment','relative_after_tweakback'),result)
    return result


def annotate_astrometry(path,solution,validation,relative):
    with fits.open(path,mode='update') as hd:
        ok=validation['status']=='validated'
        for h in hd:
            if h.name not in ('PRIMARY','SCI','WHT','CTX'):continue
            h.header['ASTSTAT']='VALIDATED' if ok else 'REVIEW'
            h.header['ASTCORR']=(solution.get('correction_applied',False),'Absolute WCS correction was propagated')
            h.header['ASTREF']=solution['reference'];h.header['ASTRELS']=','.join(solution.get('reference_releases',[]))
            h.header['ASTMODE']=solution['reference_mode']
            h.header['ASTNREF']=validation['matched'];h.header['ASTEPCH']=(solution['observation_epoch'],'Julian year of effective HST epoch')
            h.header['ASTRMS']=(validation['rms_arcsec'],'RMS residual in arcsec; see ASTNREF')
            h.header['ASTWCS']=solution['corrected_wcs_name'];h.header['HPPASS']='FINAL'
            h.header['ASTREL']=relative['status']
            if solution.get('correction_applied',False):h.header.add_history('Absolute reference correction propagated to FLC copies with STScI TweakBack before final AstroDrizzle')
            else:h.header.add_history('Existing exposure WCS validated against external reference; no TweakBack correction was needed')
            if solution['reference']=='legacy_dr10':h.header.add_history('Legacy reference may include DR9-North; actual catalog releases are preserved in alignment catalog')
        for h in hd:h.add_checksum()


def finalize(provisional,relative_inputs,rows,paths,options,drizzle_options,relative,cr_threshold=.02):
    if relative.get('status') not in ('validated','single_exposure'):
        return [],dict(status='relative_astrometry_unvalidated',relative_astrometry=relative,astrometry=dict(status='unvalidated',reason='Relative WCS gate failed; no absolute alignment or TweakBack'))
    solution,reference=absolute_alignment(provisional,rows,paths,options)
    result=dict(astrometry=solution,relative_astrometry=relative)
    if solution['status']!='validated':
        result.update(status='astrometry_unvalidated',final=dict(validation_status='review_required',limitations=['No accepted absolute reference solution; provisional image retained; no new science-ready final DRC']))
        return [],result
    cache=paths.cache();final_destination=paths.product('drc')
    with paths.scratch('final-pass-') as tmp:
        work=Path(tmp);inputs=[];final=work/(paths.group_id+'_candidate_drc.fits')
        for source in relative_inputs:
            p=work/source.name;paths.copy(source,p);inputs.append(p)
        result['tweakback']=tweakback_inputs(provisional,relative_inputs,inputs,solution,paths,options['tweakback'])
        transfer=validate_relative_transfer(relative_inputs,inputs,relative,options['relative'],paths)
        result['relative_after_tweakback']=transfer
        if transfer['status'] not in ('validated','single_exposure'):
            result.update(status='relative_astrometry_unvalidated',final=dict(validation_status='review_required',limitations=['Internal registration failed after TweakBack; no final drizzle']))
            return [],result
        headerlets=[Path(h['path']) for h in result['tweakback']['headerlets']]
        final_options=dict(drizzle_options)
        if options['cross_filter']['common_grid'] and solution.get('cross_filter_anchor'):
            for name in ('final_scale','final_rot','final_ra','final_dec','final_outnx','final_outny'):final_options.pop(name,None)
            anchor=solution['cross_filter_anchor'];version=science_header(anchor).get('EXTVER',1)
            final_options.update(final_refimage=f'{anchor}[SCI,{version}]',final_wcs=True)
        grid_inputs=[Path(solution['cross_filter_anchor'])] if options['cross_filter']['common_grid'] and solution.get('cross_filter_anchor') else []
        key=cache.key(relative_inputs+headerlets+grid_inputs+[Path(__file__).parents[1]/'drizzle.py'],dict(parameters=final_options,versions=versions()))
        crfile=paths.file('intermediate','drizzle','final_cr')
        final_masks=[paths.exposure_file('intermediate','drizzle',e.metadata['ROOTNAME'],'final_crmask') for e in rows]
        final_cached=cache.valid('uvis_final',key)
        if final_cached:paths.copy(final_destination,final)
        if not final_cached:
            previous=Path.cwd()
            try:
                import os
                os.chdir(work)
                prefix=work/paths.group_id
                drizzle(inputs,prefix,dict(final_options,driz_separate=len(inputs)>1,median=len(inputs)>1,blot=len(inputs)>1,
                    driz_cr=len(inputs)>1,driz_combine=True,runfile='final_astrometric_astrodrizzle.log'))
                if len(inputs)>1:apply_cr_masks(inputs)
                rejected={}
                for e,p,mask in zip(rows,inputs,final_masks):
                    rejected[e.path.name]=cr_report(e.path,p,mask,paths.area('diagnostics','drizzle'),cr_threshold,paths=paths)
                save(crfile,rejected)
                product=prefix.with_name(prefix.name+'_drc.fits')
                annotate_product(product,final_options['final_wht_type'],len(inputs),relative['status'],final_options['final_bits']);validate_drizzle(product)
                paths.copy(product,final)
                for f in work.glob('*.log'):paths.copy(f,paths.file('logs',None,f.stem,'log'))
            finally:os.chdir(previous)
        if solution['reference_mode']=='cross_filter':
            from .sources import sources
            detected=sources(final,options['relative']['detection'],True,6,paths=paths,label='final_sources')
        else:detected=detect_sources(final,options['detection'],paths=paths,label='final_sources')
        validation,matched=validate(detected,reference,solution['thresholds'],science_shape(final))
        validation['independence']='Final image residuals can reuse source-identification/fit objects; no held-out absolute accuracy is claimed'
        validation['wcs_identity']=science_header(final).get('WCSNAME')
        if solution['reference_mode']=='cross_filter':
            validation['basis']='relative to validated optical anchor; direct absolute accuracy inherited'
            validation['anchor_absolute_validation']=solution['inherited_absolute_validation']
            validation['relative_to_anchor_rms_arcsec']=validation['rms_arcsec'];validation['absolute_rms_arcsec']=None
            anchor_error=(solution['inherited_absolute_validation'] or {}).get('uncertainty_arcsec')
            if anchor_error is not None and validation['uncertainty_arcsec'] is not None:validation['uncertainty_arcsec']=float(np.hypot(anchor_error,validation['uncertainty_arcsec']))
        else:
            validation['absolute_rms_arcsec']=validation['rms_arcsec']
            validation['basis']='final DRC against independently retrieved external reference; fit-star scatter is not independent holdout accuracy'
        matched.write(paths.file('diagnostics','alignment','final_matches','ecsv'),format='ascii.ecsv',overwrite=True)
        if len(matched):image_plot([np.column_stack([matched['dx_arcsec'],matched['dy_arcsec']])],['Final residuals: east/north arcsec'],paths.file('diagnostics','alignment','final_residuals','png'))
        epoch_checks=[]
        for e,p in zip(rows,inputs):
            epoch=observation_epoch([e])[0]
            refs_at_epoch=gaia_at_epoch(reference,epoch) if solution['reference_mode']=='gaia_dr3' else (legacy_at_epoch(reference,epoch) if solution['reference_mode']=='legacy_dr10' else reference)
            measured,_=validate(detect_sources(p,options['detection'],paths=paths,label=p.stem+'_final_sources'),refs_at_epoch,solution['thresholds'])
            measured.update(exposure=e.path.name,observation_epoch=epoch,reference_epoch_handling=('Gaia propagated to this exposure epoch' if solution['reference_mode']=='gaia_dr3' else 'Legacy/anchor epoch limitations retained'))
            epoch_checks.append(measured)
        validation['per_exposure_epoch_checks']=epoch_checks
        validation['relative_registration']=transfer
        validation['exposure_registration_status']=transfer['status']
        validation['individual_absolute_status']='validated' if all(r['status']=='validated' for r in epoch_checks) else 'review_required'
        validation['warnings']=[]
        if validation['individual_absolute_status']!='validated':
            if all(r['status'] in ('validated','diagnostic_single_reference','insufficient_verified_references') for r in epoch_checks):
                validation['warnings'].append('Some FLCs have insufficient external objects for an independent absolute check; their WCS is supported by validated internal registration and the calibrated mosaic')
            else:validation['warnings'].append('Some individual exposure/reference residuals failed validation; inspect relative registration and source centroids')
        if relative['status'] not in ('validated','single_exposure'):validation['warnings'].append('Exposure-to-exposure relative registration remains unresolved')
        if solution['exposure_epoch_span_years']>options['cross_filter']['max_epoch_difference_years'] and any(r['status']!='validated' for r in epoch_checks):
            validation['status']='unvalidated_cross_epoch_registration'
        solution['final_validation']=validation
        if validation['status']!='validated':solution['status']='final_validation_failed'
        annotate_astrometry(final,solution,validation,relative)
        # Annotation changes the checksum; keep the final-stage receipt current.
        result['astrometry']=solution
        result['cr']=json.loads(crfile.read_text())
        result['final']=final_report(final,relative_inputs,paths.area('diagnostics','drizzle'),paths=paths)
        if len(rows)==1:result['final']['limitations'].append('Single exposure: no reliable multi-exposure CR rejection; NOT fully CR cleaned')
        result['status']='review_required' if validation['status']!='validated' or relative['status'] not in ('validated','single_exposure') else 'astrometry_validated_image_review_required'
        save(paths.file('diagnostics','alignment','final_astrometry'),solution)
        from ..chips import coverage_report
        result['final_chips']=coverage_report(final,inputs)
        from ..source_qa import validate_source_rejection,publish_validated_image
        source_validation=validate_source_rejection(final,inputs,paths,relative,final_options)
        result['source_validation']=source_validation
        for chips in result['cr'].values():
            for chip in chips:chip['persistent_source_validation']=source_validation['status']
        excessive_rejection=any((c.get('current_cr_fraction') or 0)>.1 for exposure in result['cr'].values() for c in exposure)
        image_warning=(source_validation['status'] not in ('validated','not_applicable_single_exposure') or
            excessive_rejection or bool(result['final'].get('warnings')))
        if image_warning:
            result['status']='image_quality_review_required'
            result['final']['limitations'].append('Persistent-source flux/shape/coverage or excessive-CR check requires review; previous production product preserved')
        else:
            result['final']['validation_status']='validated'
            result['status']=('science_ready' if source_validation['status']=='validated' else 'complete_with_limitations') if validation['status']=='validated' else 'review_required'
        if validation['status']!='validated' or image_warning:
            review=paths.file('intermediate','drizzle','unvalidated_final','fits');paths.copy(final,review)
            result['unvalidated_final']=str(review)
            return [],result
        if source_validation['status']=='validated':publish_validated_image(final,paths,source_validation)
        else:paths.copy(final,final_destination)
        cache.record('uvis_final',key,[final_destination,crfile]+final_masks)
        return [final_destination],result
