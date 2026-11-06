"""Measured, internal registration of provisional images before absolute fitting.

Catalogs contain native HST centroids. Pairwise connections, including rejected
ones, are retained; a disconnected filter never inherits an unmeasured frame.
"""
import json
import tempfile
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from stwcs.wcsutil import altwcs
from .sources import project
from .matching import detect_sources
from .relative_fit import match, quality, select_model
from .relative_engine import apply_native_fit
from .wcs import signature, detector_wcs
from ...common.paths import ReductionPaths
from ...common.target_runner import save, target_lock


def cross_sources(image, options, *, paths=None):
    """Reuse stable compact centroids, with all existing masked-core checks."""
    detection = dict(options['detection'], **options['cross_filter'].get('detection', {}))
    detection.setdefault('compact_kernel_widths', [3.5, 6.])
    detection['method'] = options['cross_filter'].get('source_method', 'compact_stable')
    def measure():
        if detection['method'] != 'dao_compact':return detect_sources(image,detection,paths=paths,label='cross_filter_sources')
        from .full_field import merge_detectors
        compact=detect_sources(image,dict(detection,method='compact_stable'),paths=paths,label='cross_filter_compact')
        stars=detect_sources(image,dict(detection,method='dao'),paths=paths,label='cross_filter_DAO')
        return merge_detectors([stars,compact],options['relative'].get('dedup_radius_pix',2.))
    if paths is None:
        return measure()
    output = paths.file('intermediate', 'alignment', 'cross_filter_sources', 'ecsv')
    cache = paths.cache()
    deps = [Path(__file__).with_name(n) for n in ('cross_filter.py','sources.py', 'matching.py', 'centroids.py', 'wcs.py','full_field.py')]
    key = cache.key(deps, dict(image=signature(image), detection=detection))
    if cache.valid('uvis_cross_sources', key):
        return Table.read(output, format='ascii.ecsv')
    cat = measure()
    cat.write(output, format='ascii.ecsv', overwrite=True)
    cache.record('uvis_cross_sources', key, [output])
    return cat


def assess_pair(reference, image, a, b, options):
    """STScI histogram associations and existing constrained model selection."""
    with fits.open(reference) as hd:
        tangent = detector_wcs(hd, hd['SCI'].header.get('EXTVER', 1))
    with fits.open(image) as hd:
        target = detector_wcs(hd, hd['SCI'].header.get('EXTVER', 1))
    x, y = project(a, reference, tangent), project(b, image, tangent)
    ri, ii = match(x, y, options, True)
    ratio = np.sqrt(abs(np.linalg.det(target.pixel_scale_matrix) / np.linalg.det(tangent.pixel_scale_matrix)))
    error = np.hypot(a['centroid_error_pix'][ri], b['centroid_error_pix'][ii] * ratio)
    comparison = select_model(x[ri], y[ii], error, options)
    fit = comparison['selected']
    # Blind histogram candidates are not verified identities. Exact/minimal
    # fits need predictive redundancy, not a successful numerical solution.
    supported = fit['status'] in ('validated', 'refinement_supported') and fit.get('loo_rms_pix') is not None
    pairs = Table(dict(reference_index=ri, source_index=ii,
        reference_x=np.asarray(a['x'][ri]), reference_y=np.asarray(a['y'][ri]),
        source_x=np.asarray(b['x'][ii]), source_y=np.asarray(b['y'][ii]),
        reference_error_pix=np.asarray(a['centroid_error_pix'][ri]), source_error_pix=np.asarray(b['centroid_error_pix'][ii]),
        dx_pix=(x[ri]-y[ii])[:, 0], dy_pix=(x[ri]-y[ii])[:, 1], centroid_error_pix=error,
        accepted=np.asarray(fit.get('fit_inlier_mask', np.zeros(len(ri), bool)))))
    return dict(status='supported' if supported else 'unvalidated', before=fit,
        model_comparison=comparison, reason=None if supported else 'No stable, redundant common-source transformation'), pairs, tangent


def connection_plan(keys, edges, counts, requested='auto'):
    """Choose anchors by reliable overlap; retain separate disconnected frames."""
    pending=set(keys); records={}
    while pending:
        score=lambda k:(sum(edges.get((k,j),{}).get('status')=='supported' for j in pending-{k}),
                        sum(edges.get((k,j),{}).get('before',{}).get('robust_inliers',0) for j in pending-{k} if edges.get((k,j),{}).get('status')=='supported'), counts.get(k,0))
        anchor=requested if requested in pending else max(sorted(pending),key=score)
        records[anchor]=dict(status='anchor',anchor_group=anchor,network_root=anchor,depth=0)
        pending.remove(anchor); reached={anchor}
        while pending:
            candidates=[(r,i) for r in reached for i in pending if edges.get((r,i),{}).get('status')=='supported']
            if not candidates:break
            def rank(pair):
                r,i=pair; f=edges[pair]['before']
                return (records[r]['depth'], f.get('transform_uncertainty_pix',np.inf), -f.get('robust_inliers',0), r,i)
            parent,child=min(candidates,key=rank)
            records[child]=dict(status='validated',anchor_group=parent,network_root=anchor,depth=records[parent]['depth']+1)
            pending.remove(child);reached.add(child)
        for k in reached:records[k]['connected_groups']=sorted(reached)
    return records


def align_cross_filters(directory, selected=None, options=None):
    """Measure all usable filter pairs on copies, preserving original SCI/WCS."""
    from ..config import settings
    options=settings({'astrometry':options or {}})['astrometry']
    directory=Path(directory)
    with target_lock(directory):
        state=json.loads((directory/'manifest.json').read_text())
        selected=list(state['filters']) if selected is None else list(selected)
        if set(selected)-set(state['filters']):raise ValueError('Unknown filter/setup group')
        paths={k:ReductionPaths.from_directory(directory,r['filter'],r.get('setup'),target=state.get('target')) for k,r in state['filters'].items()}
        available={}; catalogs={}; stellar={}; errors={}; normalized={}; effective={}; edges={}; pairs={}; outputs={}
        with tempfile.TemporaryDirectory(prefix='hstpipe-cross-filter-') as tmp:
            for key,r in state['filters'].items():
                effective[key]=settings(r.get('configuration',{}).get('options',{}))['astrometry'] if r.get('configuration') else options
                a=r.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}))
                p=paths[key].file('intermediate','drizzle','provisional','fits')
                if rel.get('status') not in ('validated','single_exposure') or not p.exists():continue
                if not effective[key]['cross_filter']['enabled']:continue
                available[key]=p
                try:
                    copy=Path(tmp)/p.name;paths[key].copy(p,copy)
                    with fits.open(copy) as hd:v=hd['SCI'].header.get('EXTVER',1);original=hd['SCI'].header.get('HPORIGK')
                    if original:altwcs.restoreWCS(str(copy),ext=[('SCI',v)],wcskey=original)
                    normalized[key]=copy;catalogs[key]=cross_sources(copy,effective[key],paths=paths[key])
                except Exception as exc:errors[key]=f'{type(exc).__name__}: {exc}'
            keys=sorted(k for k in catalogs if len(catalogs[k]))
            signatures={k:signature(p) for k,p in available.items()}
            deps=[Path(__file__).with_name(n) for n in ('cross_filter.py','relative_fit.py','relative_engine.py','spatial_residuals.py','matching.py','centroids.py')]
            for ref in keys:
                for key in keys:
                    if key==ref:continue
                    output=paths[key].file('intermediate','alignment','cross_to_'+ref)
                    tablefile=paths[key].file('intermediate','alignment','cross_to_'+ref+'_matches','ecsv')
                    cache=paths[key].cache();stage='uvis_cross_to_'+ref
                    cachekey=cache.key(deps+[paths[k].file('intermediate','alignment','cross_filter_sources','ecsv') for k in (ref,key)],dict(reference=signatures[ref],target=signatures[key],relative=effective[key]['relative']))
                    if cache.valid(stage,cachekey):
                        edge=json.loads(output.read_text());table=Table.read(tablefile,format='ascii.ecsv')
                    else:
                        try:
                            edge,table,_=assess_pair(normalized[ref],normalized[key],catalogs[ref],catalogs[key],effective[key]['relative'])
                            edge['source_selection']='stable_compact'
                            if edge['status']!='supported':
                                original=edge.copy()
                                for k in (ref,key):
                                    if k not in stellar:
                                        stellar[k]=detect_sources(normalized[k],dict(effective[k]['relative']['detection'],method='dao'),paths=paths[k],label='cross_filter_DAO')
                                        stellar[k].write(paths[k].file('intermediate','alignment','cross_filter_DAO_sources','ecsv'),format='ascii.ecsv',overwrite=True)
                                alternative,alternate_table,_=assess_pair(normalized[ref],normalized[key],stellar[ref],stellar[key],effective[key]['relative'])
                                if alternative['status']=='supported':
                                    edge,table=alternative,alternate_table
                                    edge.update(source_selection='DAO_stellar_fallback',compact_attempt=original)
                                else:edge['stellar_attempt']=alternative
                            edge['source_catalog_counts']=[len(stellar[k]) if edge['source_selection']=='DAO_stellar_fallback' else len(catalogs[k]) for k in (ref,key)]
                        except Exception as exc:edge=dict(status='unvalidated',reason=f'{type(exc).__name__}: {exc}');table=Table(names=['reference_index','source_index'],dtype=['i8','i8'])
                        table.write(tablefile,format='ascii.ecsv',overwrite=True);save(output,edge);cache.record(stage,cachekey,[output,tablefile])
                    edges[ref,key]=edge;pairs[ref,key]=table
            requested=options['cross_filter']['reference_filter']
            plan=connection_plan(keys,edges,{k:len(catalogs[k]) for k in keys},requested)
            for key in selected:
                record=dict(stage='cross_filter_alignment',precedes='absolute_alignment',external_catalog_used=False,status='skipped')
                if key not in available:record['reason']='Within-filter registration unvalidated, provisional unavailable, or cross-filter disabled'
                elif key in errors:record.update(status='unvalidated',reason='Source catalog unavailable',error=errors[key])
                elif key not in plan:record.update(status='unvalidated',reason='No reliable HST source detections')
                else:
                    record.update(plan[key],source_count=len(catalogs[key]),target_signature=signatures[key])
                    record['pairwise_attempts']={ref:edges[ref,key] for ref in keys if ref!=key}
                    if len(available)<2:record.update(status='not_applicable',reason='Only one validated provisional image')
                    elif record['status']=='anchor' and len(record['connected_groups'])==1:
                        record.update(status='unvalidated',reason='No reliable cross-filter connection; separate absolute solution required')
                        attempts=[(ref,e) for (ref,dest),e in edges.items() if dest==key and e.get('before')]
                        if attempts:
                            ref,best=max(attempts,key=lambda item:item[1]['before'].get('robust_inliers',0))
                            record.update(before=best['before'],best_attempt_reference=ref,model_comparison=best.get('model_comparison'))
                    elif record['status']=='validated':
                        ref=record['anchor_group'];edge=edges[ref,key];table=pairs[ref,key]
                        record.update(edge);record.update(status='validated',anchor_group=ref,anchor_filter=state['filters'][ref]['filter'],anchor_signature=signatures[ref],reference_provisional=str(available[ref]))
                        output=paths[key].file('intermediate','alignment','cross_filter_matches','ecsv');table.write(output,format='ascii.ecsv',overwrite=True)
                        record.update(matches=str(output),geometry=edge['before']['geometry'],method='tweakwcs_internal_catalogs',correction_required=edge['before']['status']!='validated',
                            validation=edge['before'],validation_scope='Matched compact sources with leave-one-out prediction; wavelength-dependent morphology remains a limitation')
                        # Verify the WCS application itself on a disposable copy.
                        fit=edge['before'];copy=Path(tmp)/('fit_'+normalized[key].name);paths[key].copy(normalized[key],copy)
                        with fits.open(normalized[ref]) as hd:tangent=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
                        try:
                            if record['correction_required']:apply_native_fit(copy,tangent,fit)
                            x=np.column_stack([table['reference_x'],table['reference_y']])
                            with fits.open(copy) as hd:w=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
                            y=tangent.wcs_world2pix(w.all_pix2world(np.column_stack([table['source_x'],table['source_y']]),0),0)
                            validation=quality(x,y,table['centroid_error_pix'],effective[key]['relative'])
                            record.update(validation=validation,status='validated' if validation['status']=='validated' else 'unvalidated',transformations=[fit] if record['correction_required'] else [])
                            diagnostic=paths[key].file('diagnostics','alignment','cross_filter_sources','png')
                            plot_pairs(normalized[ref],normalized[key],copy,table,diagnostic,ref,key)
                            record['diagnostic']=str(diagnostic)
                        except Exception as exc:record.update(status='unvalidated',reason=f'WCS application check failed: {exc}')
                outputs[key]=record
            # A child cannot depend on an unsuccessful WCS application.
            for key in sorted(outputs,key=lambda k:outputs[k].get('depth',0)):
                record=outputs[key];parent=record.get('anchor_group')
                if record['status']=='validated' and parent in outputs and outputs[parent]['status'] not in ('validated','anchor'):
                    record.update(status='unvalidated',reason='Parent cross-filter frame failed validation')
                p=paths[key];r=state['filters'][key];output=p.file('intermediate','alignment','cross_filter');save(output,record)
                a=r.setdefault('scientific_assessment',{});a['cross_filter_alignment']=record
                a.setdefault('stages',{})['cross_filter_alignment']=dict(status=record['status'],outputs=[str(output)],external_catalog_used=False)
                qa=p.file('diagnostics','drizzle','qa_summary')
                if qa.exists():save(qa,a)
                from ...common.provenance import identity
                artifacts={item['path']:item for item in r.get('artifacts',[]) if Path(item['path']).exists()}
                for artifact in p._files:
                    if artifact.is_file() and artifact.is_relative_to(directory):artifacts[str(artifact)]=identity(artifact)
                r['artifacts']=list(artifacts.values())
                r['intermediate_files']=[item for item in artifacts.values() if '/intermediate/' in item['path']]
            save(directory/'manifest.json',state)
    return outputs


def plot_pairs(reference,before,after,pairs,output,reference_name,filter_name):
    """Native-pixel cutouts; crosses are one fixed reference sky position."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm
    selected=pairs[np.asarray(pairs['accepted'],bool)][:3]
    if not len(selected):return
    output.parent.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(len(selected),3,figsize=(8,2.5*len(selected)),squeeze=False)
    with fits.open(reference) as rh, fits.open(before) as bh, fits.open(after) as ah:
        rw=detector_wcs(rh,rh['SCI'].header.get('EXTVER',1))
        for i,row in enumerate(selected):
            sky=rw.all_pix2world([[row['reference_x'],row['reference_y']]],0)
            for j,(hd,title) in enumerate([(rh,reference_name+' reference'),(bh,filter_name+' before'),(ah,filter_name+' candidate')]):
                w=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1));x,y=w.all_world2pix(sky,0)[0];ix,iy=int(round(x)),int(round(y));r=12
                data=hd['SCI'].data;sy=slice(max(0,iy-r),min(data.shape[0],iy+r+1));sx=slice(max(0,ix-r),min(data.shape[1],ix+r+1));cut=data[sy,sx]
                ax=axes[i,j]
                if cut.size and np.isfinite(cut).any():ax.imshow(cut,origin='lower',extent=(sx.start-.5,sx.stop-.5,sy.start-.5,sy.stop-.5),cmap='gray',norm=simple_norm(cut[np.isfinite(cut)],'asinh',percent=99))
                ax.plot(x,y,'+',color='cyan');ax.set_title(title,fontsize=9)
    fig.suptitle('Same reference sky position; native pixels, no coadd or image recentering',fontsize=10);fig.tight_layout();fig.savefig(output,dpi=100);plt.close(fig)


internal_cross_filter=align_cross_filters


def transfer_absolute_frame(provisional, rows, paths, options, anchor_data):
    """Refit the saved HST pairs to the now-absolute parent WCS, then TweakBack.

    Only provisional WCS headers change here. The original drizzle WCS stays
    archived so the downstream official TweakBack applies the complete combined
    correction once to every contributing FLC SCI extension.
    """
    from .alignment import observation_epoch, criteria
    from .absolute_matching import validate_absolute
    from .matching import science_shape, science_header
    from ...common.provenance import identity
    checkpoint=paths.file('intermediate','alignment','cross_filter')
    record=json.loads(checkpoint.read_text())
    if record.get('status')!='validated' or signature(provisional)!=record['target_signature']:
        raise ValueError('Missing or stale validated cross-filter relationship')
    anchor_image,anchor_report=anchor_data
    reference=Path(anchor_report.get('provisional_reference') or anchor_image)
    if signature(reference)!=record['anchor_signature']:
        raise ValueError('Parent provisional pixels or original WCS changed')
    pairfile=Path(record['matches']);pairs=Table.read(pairfile,format='ascii.ecsv')
    pairs=pairs[np.asarray(pairs['accepted'],bool)]
    if len(pairs)<2:raise ValueError('No redundant cross-filter constraints')
    cache=paths.cache();output=paths.file('intermediate','alignment','absolute')
    refout=paths.file('intermediate','alignment','anchor_sources','ecsv')
    dependencies=[checkpoint,pairfile,reference,anchor_image]
    deps=dependencies+[Path(__file__),Path(__file__).with_name('relative_fit.py'),Path(__file__).with_name('relative_engine.py')]
    key=cache.key(deps,dict(image=signature(provisional),options=options['cross_filter'],relative=options['relative'],anchor_validation=anchor_report.get('final_validation')))
    if cache.valid('uvis_absolute',key):return json.loads(output.read_text()),Table.read(refout,format='ascii.ecsv')
    epoch,span=observation_epoch(rows)
    with paths.scratch('cross-transfer-') as tmp:
        candidate=Path(tmp)/provisional.name;paths.copy(provisional,candidate)
        with fits.open(candidate) as hd:v=hd['SCI'].header.get('EXTVER',1);original=hd['SCI'].header.get('HPORIGK')
        if original:altwcs.restoreWCS(str(candidate),ext=[('SCI',v)],wcskey=original)
        with fits.open(reference) as hd:tangent=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
        with fits.open(candidate) as hd:w=detector_wcs(hd,v)
        refxy=np.column_stack([pairs['reference_x'],pairs['reference_y']])
        xy=np.column_stack([pairs['source_x'],pairs['source_y']])
        world=tangent.all_pix2world(refxy,0);target=tangent.wcs_world2pix(w.all_pix2world(xy,0),0)
        errors=np.asarray(pairs['centroid_error_pix'],float)
        fit=quality(refxy,target,errors,options['relative'],record['geometry'])
        if fit['status'] not in ('validated','refinement_supported') or fit.get('loo_rms_pix') is None:
            raise ValueError('Cross-filter transfer fit failed: '+str(fit))
        orig_key,orig_name=altwcs.archive_wcs(str(candidate),ext=[('SCI',v)],wcsname='HP-PROVISIONAL',mode=altwcs.ArchiveMode.AUTO_RENAME)
        corrected_name='HP-CROSS-ABS'
        names=set(altwcs.wcsnames(str(candidate),ext=('SCI',v)).values())
        while corrected_name in names:corrected_name+='X'
        # Even a below-tolerance shift is applied to establish the same measured
        # reference frame; its significance/fit validity was checked above.
        apply_native_fit(candidate,tangent,fit,wcsname=corrected_name)
        with fits.open(candidate) as hd:updated=detector_wcs(hd,v);mapped=tangent.wcs_world2pix(updated.all_pix2world(xy,0),0)
        after=quality(refxy,mapped,errors,options['relative'])
        if after['status']!='validated':raise ValueError('WCS propagation failed fixed-pair validation: '+str(after))
        refs=Table(dict(source_id=[f'{record["anchor_group"]}-{i}' for i in range(len(pairs))],ra=world[:,0],dec=world[:,1]))
        refs.meta.update(anchor=str(anchor_image),reference_provisional=str(reference),epoch=epoch,reference_frame='Validated HST cross-filter centroids; absolute frame inherited')
        measured=Table(dict(x=xy[:,0],y=xy[:,1],ra=updated.all_pix2world(xy,0)[:,0],dec=updated.all_pix2world(xy,0)[:,1],centroid_error_pix=np.asarray(pairs['source_error_pix'])))
        thresholds=criteria(options,anchor_report['reference']);scale=np.sqrt(abs(np.linalg.det(tangent.pixel_scale_matrix)))*3600
        thresholds.update(min_matches=1,min_span_fraction=0.,max_rms_arcsec=min(thresholds['max_rms_arcsec'],options['relative']['max_scatter_pix']*scale),max_offset_arcsec=options['relative']['max_offset_pix']*scale,pixel_scale_arcsec=float(scale))
        validation,matched=validate_absolute(measured,refs,thresholds,science_shape(candidate))
        if validation['status']!='validated':raise ValueError('Absolute-frame transfer residual check failed: '+str(validation))
        validation['fixed_pair_validation']=after
        refs.write(refout,format='ascii.ecsv',overwrite=True)
        with fits.open(candidate,mode='update') as hd:
            hd['SCI'].header['HPORIGK']=orig_key;hd[0].header['ASTREF']=anchor_report['reference'];hd[0].header['ASTSTAT']='PROVISIONAL'
        solution=dict(status='validated',provisional_image_path=str(provisional),reference=anchor_report['reference'],reference_mode='cross_filter',reference_catalog=str(refout),reference_releases=anchor_report.get('reference_releases',[]),
            observation_epoch=epoch,exposure_epoch_span_years=span,validation=validation,thresholds=thresholds,original_wcs_key=orig_key,original_wcs_name=orig_name,
            corrected_wcs_name=corrected_name,correction_applied=True,fit_geometry=record['geometry'],transformations=[fit],cross_filter_anchor=str(anchor_image),
            inherited_absolute_validation=anchor_report.get('final_validation',anchor_report.get('validation')),reference_dependencies=[str(p) for p in dependencies],
            cross_filter_status='validated_saved_pairs',cross_filter_alignment=record,attempts=[dict(reference='cross_filter',status='validated',initial_validation=fit,validation=validation)],warnings=[],
            transformation_units='Parent tangent-plane pixels; counterclockwise rotation degrees; dimensionless scale',
            independence='Saved HST associations; leave-one-out stability is not an independent external accuracy measurement')
        paths.copy(candidate,provisional);save(output,solution)
        matchout=paths.file('intermediate','alignment','absolute_matches','ecsv');matched.write(matchout,format='ascii.ecsv',overwrite=True)
        cache.record('uvis_absolute',key,[provisional,output,refout,matchout]+dependencies)
        return solution,refs


def validate_existing_frame(product, paths, options, anchor_data):
    """Retain an already-reviewed final image if it meets the measured frame.

    This is a read-only final-image check, not acceptance of a new candidate's
    image-QA failures. Published SCI pixels and their review remain unchanged.
    """
    from .alignment import criteria
    from .absolute_matching import validate_absolute
    from .matching import science_shape
    from ...common.provenance import identity
    product=Path(product)
    with fits.open(product) as hd:
        if hd[0].header.get('HPQASTAT') not in ('VALIDATED','REVIEWED') or hd['SCI'].header.get('ASTSTAT')!='VALIDATED':
            return dict(status='unvalidated',reason='Published product lacks validated astrometry and image QA')
    checkpoint=paths.file('intermediate','alignment','cross_filter');record=json.loads(checkpoint.read_text())
    if record.get('status')!='validated' or not record.get('matches'):return dict(status='unvalidated',reason='No validated internal filter connection')
    reference=Path(anchor_data[1].get('provisional_reference') or anchor_data[0])
    if signature(reference)!=record['anchor_signature']:return dict(status='unvalidated',reason='Stale parent provisional frame')
    pairs=Table.read(record['matches'],format='ascii.ecsv');pairs=pairs[np.asarray(pairs['accepted'],bool)]
    with fits.open(reference) as hd:w=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
    sky=w.all_pix2world(np.column_stack([pairs['reference_x'],pairs['reference_y']]),0)
    refs=Table(dict(source_id=[f'parent-{i}' for i in range(len(pairs))],ra=sky[:,0],dec=sky[:,1]))
    scale=float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600)
    thresholds=criteria(options,anchor_data[1]['reference'])
    thresholds.update(min_matches=1,min_span_fraction=0.,max_rms_arcsec=min(thresholds['max_rms_arcsec'],options['relative']['max_scatter_pix']*scale),max_offset_arcsec=options['relative']['max_offset_pix']*scale,pixel_scale_arcsec=scale)
    method=record.get('source_selection');detection_options=options
    if method=='DAO_stellar_fallback':detection_options=dict(options,cross_filter=dict(options['cross_filter'],source_method='dao'))
    cache=paths.cache();output=paths.file('intermediate','alignment','cross_filter_existing_final')
    key=cache.key([product,reference,checkpoint,Path(record['matches']),Path(__file__),*[Path(__file__).with_name(n) for n in ('absolute_matching.py','matching.py','sources.py','centroids.py','wcs.py','full_field.py')]],dict(thresholds=thresholds,detection=options['detection'],cross_filter=detection_options['cross_filter']))
    if cache.valid('uvis_cross_existing_final',key):return json.loads(output.read_text())
    measured=cross_sources(product,detection_options)
    validation,matches=validate_absolute(measured,refs,thresholds,science_shape(product))
    validation.update(method='validated_existing_product' if validation['status']=='validated' else 'existing_product_unvalidated',product=identity(product),parent=identity(reference),
        relative_to_anchor_rms_arcsec=validation.get('rms_arcsec'),preserved_without_resampling=True,
        independence='Direct final-image check against measured HST parent sources; not independent absolute accuracy')
    matchfile=paths.file('intermediate','alignment','cross_filter_existing_matches','ecsv');matches.write(matchfile,format='ascii.ecsv',overwrite=True)
    save(output,validation);cache.record('uvis_cross_existing_final',key,[output,matchfile])
    return validation
