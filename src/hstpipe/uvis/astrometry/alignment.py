"""Absolute correction through TweakReg, validated before committing WCS updates."""
import json,os
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from .matching import science_header,science_shape
from ...common.target_runner import save,merge
from .catalogs import query_gaia_dr3,query_legacy_dr10,query_footprint,CatalogQueryError,CatalogSchemaError,legacy_at_epoch,cached_footprint
from .matching import detect_sources
from .absolute_matching import validate_absolute as validate,reference_audit,plot_reference_audit


def observation_epoch(rows):
    from astropy.time import Time
    times=[]
    for e in rows:
        value=e.metadata.get('EXPSTART')
        if value is not None:times.append(Time(float(value),format='mjd').tcb.jyear)
        elif e.metadata.get('DATE-OBS'):times.append(Time(e.metadata['DATE-OBS']+'T'+str(e.metadata.get('TIME-OBS') or '00:00:00')).tcb.jyear)
        else:raise ValueError('Astrometric propagation requires an actual HST observation date')
    return float(np.mean(times)),float(np.ptp(times))


from .wcs import signature


def criteria(options,reference):
    base={k:options['absolute'][k] for k in ('min_matches','min_span_fraction','search_radius_arcsec','clip_sigma','fitgeometry')}
    base['match_tolerance_arcsec']=options['absolute'].get('match_tolerance_arcsec',.2)
    return merge(base,options['catalogs'][reference]['validation'])


def _anchor(paths,epoch,options):
    state=paths.directory/'manifest.json'
    if not state.exists() or not options['cross_filter']['enabled']:return None
    records=json.loads(state.read_text()).get('filters',{})
    requested=options['cross_filter']['reference_filter']
    own=records.get(paths.group_key,{}).get('scientific_assessment',{}).get('cross_filter_alignment',{})
    if own.get('status') not in ('validated','anchor'):return None
    rank=lambda item:(item[1].get('scientific_assessment',{}).get('astrometry',{}).get('final_validation',{}).get('uncertainty_arcsec') or float('inf'),-item[1].get('scientific_assessment',{}).get('astrometry',{}).get('final_validation',{}).get('matched',0))
    for key,r in sorted(records.items(),key=rank):
        if key==paths.group_key or (requested!='auto' and r.get('filter')!=requested):continue
        assessment=r.get('scientific_assessment',{})
        other=assessment.get('cross_filter_alignment',{})
        if not (own.get('anchor_group')==key or (other.get('status') in ('validated','anchor') and own.get('anchor_group')==other.get('anchor_group'))):continue
        relative=assessment.get('relative_astrometry',assessment.get('alignment',{}))
        if relative.get('status') not in ('validated','single_exposure'):continue
        ast=assessment.get('astrometry',{})
        if ast.get('status')!='validated' or not r.get('products'):continue
        if abs(epoch-ast.get('observation_epoch',epoch))>options['cross_filter']['max_epoch_difference_years']:continue
        image=Path(r['products'][0]['path'])
        if image.exists():
            provisional=r.get('scientific_assessment',{}).get('provisional_image',{}).get('path')
            return image,dict(ast,provisional_reference=provisional) if provisional and Path(provisional).exists() else ast
    return None


def _sparse_absolute_shift(candidate,sources,refs,matched,corrected_name):
    """Native STScI shift fitter without TweakReg's degenerate catalog polygon.

    The installed TweakReg wrapper cannot construct a two-point bounding
    polygon. Use its maintained fitting/WCS-update stack, never fake references
    to make a polygon. Source positions remain measured, zero-origin pixels.
    """
    from tweakwcs.linearfit import iter_linear_fit
    from drizzlepac.updatehdr import updatewcs_with_shift
    from .wcs import detector_wcs
    with fits.open(candidate) as hd:
        w=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
        xy=np.column_stack([sources['x'],sources['y']])
        target=w.all_world2pix(np.column_stack([refs['ra'],refs['dec']]),0)
    errors=np.asarray(matched['positional_error_arcsec'],float)
    fit=iter_linear_fit(target,xy,wxy=1/np.maximum(errors,1e-8)**2,fitgeom='shift',nclip=0)
    shift=np.asarray(fit['shift'],float)
    if not np.isfinite(shift).all():raise ValueError('Nonfinite sparse absolute transformation')
    # updatehdr maps reference-frame offsets in the inverse direction to the
    # fitted image->reference shift; identity matrix fixes scale and rotation.
    updatewcs_with_shift(str(candidate),w,wcsname=corrected_name,fitgeom='shift',fit=np.eye(2),sciext='SCI',
        xsh=float(-shift[0]),ysh=float(-shift[1]),force=False)
    residual=target-xy-shift
    return dict(dx=float(-shift[0]),dy=float(-shift[1]),rotation=0.,scale=1.,
        xrms=float(np.sqrt(np.mean(residual[:,0]**2))),yrms=float(np.sqrt(np.mean(residual[:,1]**2))))


def absolute_alignment(provisional,rows,paths,options,*,reference_order=None):
    from stwcs.wcsutil import altwcs
    from drizzlepac import tweakreg
    epoch,epoch_span=observation_epoch(rows);footprint=query_footprint(provisional,options['catalogs']['margin_arcsec'])
    if footprint['radius']>options['catalogs']['max_radius_deg']:raise ValueError('Provisional WCS requires an unexpectedly large catalog cone')
    base=dict(status='unvalidated',observation_epoch=epoch,exposure_epoch_span_years=epoch_span,footprint=footprint,attempts=[],warnings=[],correction_applied=False)
    catalog_options=options['catalogs'];cache_dir=paths.catalog_directory()
    paths.area('intermediate','alignment',create=True)
    if epoch_span>options['cross_filter']['max_epoch_difference_years']:
        base['warnings'].append('Exposure epochs differ substantially: catalog is propagated to the effective mean epoch; inspect motion across epochs')
    attempts=[];anchor=_anchor(paths,epoch,options)
    base['cross_filter_status']='anchor_available' if anchor else 'no_validated_epoch_compatible_anchor'
    if anchor:attempts.append(('cross_filter',anchor))
    attempts.extend((r,None) for r in (reference_order or ('gaia_dr3','legacy_dr10')))
    for reference,anchor_data in attempts:
        attempt=dict(reference=reference,status='query_pending');base['attempts'].append(attempt)
        try:
            if reference=='cross_filter':
                anchor_image,anchor_report=anchor_data
                reference_image=Path(anchor_report.get('provisional_reference',anchor_image))
                from .sources import sources as cross_sources
                refs=cross_sources(reference_image,options['relative']['detection'],True,6)
                refs['source_id']=[f'anchor-{i}' for i in range(len(refs))]
                refs.meta.update(anchor=str(anchor_image),epoch=epoch)
                catalog_file=paths.file('intermediate','alignment','anchor_sources','ecsv');refs.write(catalog_file,format='ascii.ecsv',overwrite=True)
                reference_kind=anchor_report['reference'];thresholds=criteria(options,reference_kind)
                pixel_scale=float(np.sqrt(abs(np.linalg.det(WCS(science_header(provisional)).pixel_scale_matrix)))*3600)
                thresholds.update(min_matches=1,max_rms_arcsec=min(thresholds['max_rms_arcsec'],options['relative']['max_scatter_pix']*pixel_scale),max_offset_arcsec=options['relative']['max_offset_pix']*pixel_scale)
                dependencies=list(dict.fromkeys([catalog_file,anchor_image,reference_image,paths.file('intermediate','alignment','cross_filter')]))
                anchor_abs=anchor_report.get('final_validation',{}).get('rms_arcsec')
                attempt['anchor_absolute_rms_arcsec']=anchor_abs
            else:
                provider=options['catalogs'][reference];kwargs=dict(quality=provider['quality'],cache_dir=cache_dir,refresh=(catalog_options['refresh'] or not catalog_options['cache']),offline=catalog_options['offline'],timeout=catalog_options['timeout'],retries=catalog_options['retries'])
                if reference=='gaia_dr3':kwargs['quality']=dict(kwargs['quality'],apply_proper_motion=provider['apply_proper_motion'])
                cached=None
                if catalog_options['cache'] and not catalog_options['refresh']:
                    try:
                        candidate_ref=cached_footprint(reference,footprint['ra'],footprint['dec'],footprint['radius'],cache_dir=cache_dir,quality=kwargs['quality'],observation_epoch=epoch)
                        cached=candidate_ref
                        if not candidate_ref.meta['complete_query_coverage']:
                            attempt['coverage_warning']='Partial verified cache: fitting may use covered sources; missing sources outside cache are not nondetections'
                    except Exception:pass
                if cached is not None:refs=cached if reference=='gaia_dr3' else legacy_at_epoch(cached,epoch)
                elif reference=='gaia_dr3':refs=query_gaia_dr3(footprint['ra'],footprint['dec'],footprint['radius'],epoch,**kwargs)
                else:refs=legacy_at_epoch(query_legacy_dr10(footprint['ra'],footprint['dec'],footprint['radius'],tap_url=provider['tap_url'],**kwargs),epoch)
                catalog_file=Path(refs.meta['cache_file']) if refs.meta.get('cache_file') else paths.file('intermediate','alignment',reference+'_reference','ecsv')
                if not refs.meta.get('cache_file'):refs.write(catalog_file,format='ascii.ecsv',overwrite=True)
                dependencies=[catalog_file];thresholds=criteria(options,reference);reference_kind=reference
            attempt.update(reference_sources=len(refs),catalog=str(catalog_file),query_metadata=dict(refs.meta))
            if not len(refs):
                attempt['status']='no_catalog_coverage' if reference=='legacy_dr10' and refs.meta.get('coverage_source_count')==0 else 'no_catalog_sources'
                continue
            cache=paths.cache()
            source_code=[Path(__file__).parent/name for name in ('alignment.py','matching.py','absolute_matching.py','catalogs.py','sources.py','wcs.py','centroids.py')]
            key=cache.key(dependencies+source_code,dict(image=signature(provisional),options=options['absolute'],thresholds=thresholds,reference=reference,epoch=epoch,detection=options['detection'],cross_filter=options['cross_filter'],internal_detection=options['relative']['detection']))
            solution_file=paths.file('intermediate','alignment','absolute')
            if not catalog_options['refresh'] and cache.valid('uvis_absolute',key):return json.loads(solution_file.read_text()),refs
            with paths.scratch('absolute-') as tmp:
                temp=Path(tmp);candidate=temp/provisional.name;paths.copy(provisional,candidate)
                with fits.open(candidate) as hd:
                    v=hd['SCI'].header.get('EXTVER',1);original=hd['SCI'].header.get('HPORIGK')
                if original:altwcs.restoreWCS(str(candidate),ext=[('SCI',v)],wcskey=original)
                if reference=='cross_filter':
                    from .sources import sources as internal_sources
                    sources=internal_sources(candidate,options['relative']['detection'],True,6,paths=paths,label='absolute_cross_filter_sources')
                else:sources=detect_sources(candidate,options['detection'],paths=paths,label='absolute_'+reference+'_sources')
                shape=science_shape(candidate)
                thresholds['pixel_scale_arcsec']=float(np.sqrt(abs(np.linalg.det(WCS(science_header(candidate)).pixel_scale_matrix)))*3600)
                audit=reference_audit(candidate,sources,refs,thresholds)
                audit_file=paths.file('intermediate','alignment','absolute_'+reference+'_references','ecsv');audit.write(audit_file,format='ascii.ecsv',overwrite=True)
                attempt.update(eligible_catalog_objects=len(refs),covered_references=int(np.count_nonzero(audit['covered'])),reference_audit=str(audit_file),existing_wcs_solutions=altwcs.wcsnames(str(candidate),ext=('SCI',v)))
                initial,matched=validate(sources,refs,thresholds,shape)
                attempt['initial_validation']=initial
                matched_file=paths.file('intermediate','alignment','absolute_'+reference+'_matches','ecsv')
                matched.write(matched_file,format='ascii.ecsv',overwrite=True)
                attempt.update(verified_associations=initial['verified_associations'],fit_inliers=initial['fit_inliers'],rejected_associations=initial['rejected'],independent_validation_sources=0)
                for row in matched:
                    audit['association_status'][int(row['reference_index'])]='verified_candidate' if row['accepted'] else 'rejected_incoherent_association'
                audit.write(audit_file,format='ascii.ecsv',overwrite=True)
                plot_reference_audit(candidate,sources,refs,audit,matched,paths,reference)
                if initial['matched']<thresholds['min_matches'] or initial['matched']<2:
                    attempt['status']=initial['status'] if initial['matched'] else 'insufficient_verified_references'
                    attempt['diagnostic_shift']=initial.get('fitted_shift_arcsec');continue
                if initial['postfit_rms_arcsec']>thresholds['max_rms_arcsec'] or (initial.get('leave_one_out_rms_arcsec') or 0)>thresholds['max_rms_arcsec']:
                    attempt['status']='excessive_postfit_residuals';continue
                current=science_header(candidate).get('WCSNAME','')
                existing=reference=='gaia_dr3' and 'GAIA' in current.upper() and options['absolute']['allow_existing_gaia_wcs'] and initial['status']=='validated'
                original_name='HP-PROVISIONAL';corrected_name='HP-GAIA-DR3' if reference_kind=='gaia_dr3' else 'HP-LEGACY'
                orig_key,orig_name=altwcs.archive_wcs(str(candidate),ext=[('SCI',v)],wcsname=original_name,mode=altwcs.ArchiveMode.AUTO_RENAME)
                fit_file=None;transformations=[]
                if not existing:
                    fit_matches=matched[matched['accepted']];fit_sources=sources[np.asarray(fit_matches['source_index'],int)];fit_refs=refs[np.asarray(fit_matches['reference_index'],int)]
                    np.savetxt(temp/'sources.cat',np.column_stack([fit_sources['x']+1,fit_sources['y']+1,fit_sources['flux']]))
                    (temp/'catalogs.txt').write_text(str(candidate)+' '+str(temp/'sources.cat')+'\n')
                    np.savetxt(temp/'reference.cat',np.column_stack([fit_refs['ra'],fit_refs['dec']]))
                    if len(fit_matches)==2 and options['absolute']['fitgeometry']=='shift':
                        transformations=[_sparse_absolute_shift(candidate,fit_sources,fit_refs,fit_matches,corrected_name)]
                        attempt.update(fit_source_count=2,native_fit_catalog_inputs=2,
                            fit_backend='tweakwcs.linearfit.iter_linear_fit + drizzlepac.updatehdr.updatewcs_with_shift',
                            wrapper_limitation='Installed TweakReg cannot construct a two-point source polygon')
                        fit_file=paths.file('intermediate','alignment','absolute_'+reference+'_shifts','txt')
                        fit_file.write_text('# Native STScI sparse shift; dx dy rotation scale xrms yrms (pixels)\n'+str(candidate)+' '+' '.join(str(transformations[0][k]) for k in ('dx','dy','rotation','scale','xrms','yrms'))+'\n')
                    else:
                        previous=Path.cwd()
                        try:
                            os.chdir(temp)
                            tweakreg.TweakReg([str(candidate)],updatewcs=False,updatehdr=True,wcsname=corrected_name,
                                interactive=False,see2dplot=False,residplot='No plot',clean=False,writecat=True,
                                refimage=str(candidate),enforce_user_order=True,expand_refcat=False,
                                catfile=str(temp/'catalogs.txt'),refcat=str(temp/'reference.cat'),refxyunits='degrees',
                                fitgeometry=options['absolute']['fitgeometry'],minobj=thresholds['min_matches'],searchrad=thresholds['search_radius_arcsec'],searchunits='arcseconds',
                                tolerance=thresholds.get('match_tolerance_pixels',2 if reference_kind=='gaia_dr3' else 5),
                                nclip=0,shiftfile=True,outshifts='absolute_shifts.txt',outwcs='absolute_reference_wcs.fits',runfile='absolute_tweakreg.log')
                        finally:os.chdir(previous)
                        fit_file=paths.file('intermediate','alignment','absolute_'+reference+'_shifts','txt')
                        if (temp/'absolute_shifts.txt').exists():paths.copy(temp/'absolute_shifts.txt',fit_file)
                        nfit=0
                        for f in sorted(set(temp.glob('*fitmatch*'))|set(temp.glob('*_catalog_fit.match'))):
                            nfit+=sum(bool(line.strip()) and not line.lstrip().startswith('#') for line in f.read_text().splitlines())
                            paths.copy(f,paths.file('intermediate','alignment','absolute_'+reference+'_'+f.stem,f.suffix.lstrip('.')))
                        attempt['fit_source_count']=nfit
                        attempt['native_fit_catalog_inputs']=len(fit_matches)
                        if nfit<thresholds['min_matches']:raise RuntimeError('TweakReg fit retained too few sources')
                        if (temp/'absolute_tweakreg.log').exists():paths.copy(temp/'absolute_tweakreg.log',paths.file('logs',None,'absolute_'+reference+'_tweakreg','log'))
                        if science_header(candidate).get('WCSNAME')!=corrected_name:raise RuntimeError('TweakReg did not produce the requested corrected WCS')
                        if fit_file is None or not fit_file.exists():raise RuntimeError('Missing measured transformation')
                        for line in fit_file.read_text().splitlines():
                            if line.strip() and not line.startswith('#'):
                                fields=line.split();values=[float(x) for x in fields[1:7]]
                                if len(values)!=6 or not np.isfinite(values).all():raise RuntimeError('Unconstrained TweakReg transformation')
                                transformations.append(dict(zip(('dx','dy','rotation','scale','xrms','yrms'),values)))
                        if not transformations:raise RuntimeError('No fitted astrometric transformation')
                    corrected=internal_sources(candidate,options['relative']['detection'],True,6) if reference=='cross_filter' else detect_sources(candidate,options['detection'])
                    assessment,matched=validate(corrected,refs,thresholds,shape)
                else:assessment=initial
                # A finite correction below the ceiling is not an improvement.
                # Retain an independently acceptable original WCS when the
                # native proposal worsens external residuals on this same field.
                if not existing and initial['status']=='validated' and assessment.get('rms_arcsec',float('inf'))>initial['rms_arcsec']+1e-6:
                    attempt['rejected_transformation']=dict(reason='Corrected external RMS exceeds the already acceptable original WCS',
                        transformations=transformations,validation=assessment)
                    altwcs.restoreWCS(str(candidate),ext=[('SCI',v)],wcskey=orig_key)
                    existing=True;transformations=[];assessment=initial
                    matched=validate(sources,refs,thresholds,shape)[1]
                    attempt['accepted_wcs']='validated_existing';attempt['applied_fit_inliers']=0
                else:attempt['accepted_wcs']='validated_existing' if existing else 'validated_refined'
                matched.write(matched_file,format='ascii.ecsv',overwrite=True);attempt.update(status=assessment['status'],validation=assessment)
                if assessment['status']!='validated':continue
                if (assessment.get('leave_one_out_rms_arcsec') or 0)>thresholds['max_rms_arcsec']:
                    attempt['status']='unstable_leave_one_out_solution';continue
                with fits.open(candidate,mode='update') as hd:
                    hd['SCI'].header['HPORIGK']=orig_key
                    hd[0].header['ASTREF']=reference_kind;hd[0].header['ASTSTAT']='PROVISIONAL'
                paths.copy(candidate,provisional)
                releases=sorted(set(str(x) for x in refs['release'])) if 'release' in refs.colnames else (anchor_data[1].get('reference_releases',[]) if anchor_data else [])
                base.update(reference_releases=releases,status='validated',reference=reference_kind,reference_mode=reference,reference_catalog=str(catalog_file),
                    validation=assessment,original_wcs_key=orig_key,original_wcs_name=orig_name,
                    corrected_wcs_name=science_header(provisional).get('WCSNAME',current),correction_applied=not existing,
                    reference_dependencies=[str(p) for p in dependencies],thresholds=thresholds,fit_geometry=options['absolute']['fitgeometry'],
                    transformation_units='reference pixels, rotation degrees, dimensionless scale',
                    cross_filter_anchor=str(anchor_data[0]) if anchor_data else None,
                    inherited_absolute_validation=anchor_data[1].get('final_validation') if anchor_data else None,
                    transformation_file=str(fit_file) if fit_file else None,transformations=transformations)
                save(solution_file,base)
                matched.write(paths.file('intermediate','alignment','absolute_matches','ecsv'),format='ascii.ecsv',overwrite=True)
                cache.record('uvis_absolute',key,[provisional,solution_file,matched_file]+dependencies)
                receipt=cache.receipt('uvis_absolute');data=json.loads(receipt.read_text());data['dependencies']=[str(p) for p in dependencies];save(receipt,data)
                return base,refs
        except CatalogQueryError as exc:attempt.update(status='external_catalog_unavailable',error=str(exc))
        except CatalogSchemaError as exc:attempt.update(status='malformed_catalog',error=str(exc))
        except Exception as exc:
            import traceback
            attempt.update(status='fit_failed',error=f'{type(exc).__name__}: {exc}',failure_location=traceback.format_exc(limit=6))
    save(paths.file('intermediate','alignment','absolute'),base)
    return base,None
