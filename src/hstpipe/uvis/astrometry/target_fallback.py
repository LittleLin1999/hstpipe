"""Bounded local recovery using the existing detector, QC and native fitter."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from .centroids import QUALITY_VERSION
from .identification import associate
from ...common.target_runner import save

def target_alignment(inputs,paths,options,identification=None,catalogs=None,previous_joint=None):
    from .relative_engine import single_source_fallback,native_graph_fit,evaluate,internal_mosaic,image_quality,coadd_reference_catalog
    settings=options.get('local_fallback',{})
    if not settings.get('enabled',True):return dict(status='not_run',reason='Target-specific fallback disabled'),[]
    positions=settings.get('reference_positions',[])
    if not positions and catalogs and len(catalogs)==len(inputs) and all('common_source_id' in c.colnames for c in catalogs):
        counts={}
        for c in catalogs:
            for sid in set(map(str,c['common_source_id'])):counts[sid]=counts.get(sid,0)+1
        common={sid for sid,count in counts.items() if count>=2}
        if common:
            def score(sid):
                rows=[c[np.asarray(c['common_source_id']).astype(str)==sid][0] for c in catalogs if sid in set(map(str,c['common_source_id']))]
                stellar=all(str(r['detection_method']) in ('PSFPhotometry','aperture_profile') for r in rows) if 'detection_method' in catalogs[0].colnames else False
                return (not stellar,max(float(r['centroid_error_pix']) for r in rows),sid)
            selected=sorted(common,key=score)
            if len(selected)>=2:
                cats=[c[np.isin(np.asarray(c['common_source_id']).astype(str),selected)] for c in catalogs]
                return _refine_catalogs(inputs,cats,paths,options,dict(status='unvalidated',method='automatic_field_sources',source_counts=[len(c) for c in cats],available_common_sources=selected,external_coordinates_used_for_fit=False),previous_joint=previous_joint)
            # Sparse references constrain local translation; do not manufacture
            # a global fit with an underconstrained source distribution.
            if len(selected)==1 and counts[selected[0]]==len(inputs):
                sid=selected[0];cats=[c[np.asarray(c['common_source_id']).astype(str)==sid] for c in catalogs]
                summary=(identification or {}).get('stellar_resolution',[])
                unique={r['exposure'] for r in summary if r.get('source_id')==sid and r.get('selected') and r.get('raw_candidates')==1}
                anchor=min(range(len(inputs)),key=lambda i:(inputs[i].name not in unique,float(cats[i]['centroid_error_pix'][0])))
                order=[anchor]+[i for i in range(len(inputs)) if i!=anchor]
                record=single_source_fallback([inputs[i] for i in order],[cats[i] for i in order],options,paths,
                    reference_source=dict(source_id=sid,selection='Measured field-source quality; unique HST profile preferred as anchor',stellar=not score(sid)[0]))
                record.update(method='automatic_field_source',source_counts=[len(c) for c in cats],available_common_sources=selected,
                    result_state='local_shift_only' if record.get('status')=='local_only' else 'unvalidated')
                record['additional_source_checks']=_additional_checks([inputs[i] for i in order],[catalogs[i] for i in order],sid,record,options)
                checks=record['additional_source_checks']
                if checks:
                    record['held_out_sources']=len({c['source_id'] for c in checks})
                    record['held_out_offset_rms_pix']=float(np.sqrt(np.mean([c['radial_residual_pix']**2 for c in checks])))
                    record['all_additional_checks_consistent']=all(c['consistent'] for c in checks)
                    record['caveat']='Primary-source translation only; additional objects are held-out checks with partial exposure coverage. Rotation, scale and full-field accuracy are not independently fitted.'
                    if not record['all_additional_checks_consistent']:
                        record['existing_wcs_locally_consistent']=False;record['local_state']='additional_sources_disagree'
                save(paths.file('intermediate','alignment','target_alignment'),record)
                return record,cats
    if not positions:
        header=fits.getheader(inputs[0]);positions=[dict(ra=header['RA_TARG'],dec=header['DEC_TARG'],label='science_target')]
    if len(positions)>12:raise ValueError('Target-specific source list must remain bounded (maximum 12)')
    from .config import DEFAULTS
    local=dict(DEFAULTS['source_identification'],search_radius_arcsec=settings.get('search_radius_arcsec',2.))
    catalogs=[];measurements=[]
    for image in inputs:
        rows=[];contexts={}
        for index,position in enumerate(positions):
            sid=f'LOCAL{index+1:03d}';candidate=dict(source_id=sid,classification='point_candidate',gaia_id='',legacy_id='',
                ra=float(position['ra']),dec=float(position['dec']),release='HST local identification only')
            if not np.isfinite([candidate['ra'],candidate['dec']]).all() or not -90<=candidate['dec']<=90:raise ValueError('Invalid manual reference coordinates')
            result=associate(image,candidate,local,options['detection'],paths,label=image.stem+'_'+sid,image_contexts=contexts)
            measurements.append(result)
            if result['status']!='valid_detection':continue
            rows.append(dict(x=result['x'],y=result['y'],extver=result['extver'],source_id=sid,common_source_id=sid,
                centroid_error_pix=result['centroid_error_pix'],width_pix=result['width_pix'],snr=result['snr'],ellipticity=result['ellipticity'],
                flux=result['flux_rate']*fits.getheader(image)['EXPTIME'],centroid_valid=True))
        cat=Table(rows=rows) if rows else Table(names=['x','y','extver','source_id','common_source_id','centroid_error_pix','width_pix','snr','ellipticity','flux','centroid_valid'],dtype=['f8','f8','i8','U32','U32','f8','f8','f8','f8','f8','bool'])
        cat.meta['centroid_quality']=dict(version=QUALITY_VERSION,purpose='Independent local HST centroids passed shared pixel-level QC')
        cat.write(paths.file('intermediate','alignment',image.stem+'_target_sources','ecsv'),format='ascii.ecsv',overwrite=True)
        catalogs.append(cat)
    output=paths.file('intermediate','alignment','target_identifications','ecsv')
    Table(rows=measurements).write(output,format='ascii.ecsv',overwrite=True)
    record=dict(status='unvalidated',method='local_centroid_recovery',positions=positions,source_counts=[len(c) for c in catalogs],external_coordinates_used_for_fit=False)
    common={sid for c in catalogs for sid in map(str,c['common_source_id'])
            if sum(sid in set(map(str,other['common_source_id'])) for other in catalogs)>=2}
    if len(common)<2:
        local_result=single_source_fallback(inputs,catalogs,options,paths)
        record.update(local_result);record['result_state']='local_shift_only' if record.get('status')=='local_only' else 'unvalidated'
        save(paths.file('intermediate','alignment','target_alignment'),record);return record,catalogs
    return _refine_catalogs(inputs,catalogs,paths,options,record)


def _refine_catalogs(inputs,catalogs,paths,options,record,previous_joint=None):
    """Use the standard joint fitter, including partially observed sources.

    Reuse an already evaluated joint proposal for identical input catalogs.
    A failed image-quality assessment is never bypassed by the fallback.
    """
    from .relative_engine import joint_registration
    try:
        assessed=previous_joint if previous_joint is not None else joint_registration(inputs,catalogs,options,paths)
        record.update(validation=assessed,transformations=assessed.get('fits',[]),
                      image_quality=assessed.get('image_quality'),reused_joint_assessment=previous_joint is not None)
        if assessed.get('may_combine'):
            record.update(status='validated',result_state='validated_existing' if not assessed.get('corrections_applied') else 'validated_refined',
                          may_combine=True,scope=assessed.get('scope','common_source_constrained'))
        else:record['reason']=assessed.get('reason',assessed.get('status','Joint registration unvalidated'))
    except Exception as exc:record['error']=f'{type(exc).__name__}: {exc}'
    save(paths.file('intermediate','alignment','target_alignment'),record)
    return record,catalogs


def _additional_checks(inputs,catalogs,primary,record,options):
    """Use partial-coverage sources as independent checks of a local guide.

    No object is required in every exposure, and no secondary pixel position
    participates in the primary-star translation fit.
    """
    from .sources import project
    from .wcs import detector_wcs
    if not record.get('transformations'):return []
    with fits.open(inputs[0]) as hd:reference=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
    xy=[project(c,p,reference) for c,p in zip(catalogs,inputs)];checks=[]
    lookups=[{str(s):k for k,s in enumerate(c['common_source_id'])} for c in catalogs]
    corrections=[np.zeros(2)]+[np.asarray(next(f for f in record['transformations'] if f['exposure']==p.name)['shift']) for p in inputs[1:]]
    for i in range(len(inputs)):
        for j in range(i+1,len(inputs)):
            for sid in set(lookups[i])&set(lookups[j])-{primary}:
                a=lookups[i][sid];b=lookups[j][sid]
                residual=xy[i][a]+corrections[i]-xy[j][b]-corrections[j]
                guide_variance=sum(float(catalogs[k]['centroid_error_pix'][lookups[k][primary]])**2 for k in (i,j))
                uncertainty=float(np.sqrt(catalogs[i]['centroid_error_pix'][a]**2+catalogs[j]['centroid_error_pix'][b]**2+guide_variance))
                limit=min(options['max_scatter_pix'],max(options['scatter_floor_pix'],options['centroid_sigma']*uncertainty))
                checks.append(dict(source_id=sid,reference=inputs[i].name,exposure=inputs[j].name,residual_pix=residual.tolist(),radial_residual_pix=float(np.linalg.norm(residual)),uncertainty_pix=uncertainty,limit_pix=limit,consistent=bool(np.linalg.norm(residual)<=limit),held_out_from_primary_fit=True))
    return checks
