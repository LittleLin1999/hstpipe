"""Common HST identities from external labels and native HAP pair matching."""
import hashlib
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack
from .wcs import detector_wcs
from .sources import project
from .relative_fit import match,quality
from .centroids import QUALITY_VERSION

def build_common(inputs,catalogs,external,options,association_options=None):
    nodes=[(i,j) for i,c in enumerate(catalogs) for j in range(len(c))];parent={n:n for n in nodes}
    def root(n):
        while parent[n]!=n:parent[n]=parent[parent[n]];n=parent[n]
        return n
    def join(a,b):parent[root(b)]=root(a)
    labels={}
    for r in external:
        if r['status'] not in ('valid_detection','confirmed_common') or 'hst_source_id' not in r.colnames or np.ma.is_masked(r['hst_source_id']):continue
        i=next((i for i,p in enumerate(inputs) if p.name==r['exposure']),None)
        if i is None:continue
        idx=np.flatnonzero(np.asarray(catalogs[i]['hst_source_id']).astype(str)==str(r['hst_source_id']))
        if len(idx)==1:labels[(i,int(idx[0]))]=str(r['source_id'])
    def stellar(node):
        i,k=node
        return 'detection_method' in catalogs[i].colnames and str(catalogs[i]['detection_method'][k])=='PSFPhotometry'
    labelgroups={label:[n for n,k in labels.items() if k==label] for label in set(labels.values())}
    guided=any(len(group)>=2 and all(stellar(n) for n in group) for group in labelgroups.values())
    for label in set(labels.values()):
        group=[n for n,k in labels.items() if k==label]
        if not guided or all(stellar(n) for n in group):
            for n in group[1:]:join(group[0],n)
    with fits.open(inputs[0]) as hd:tangent=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
    xy=[project(c,p,tangent) for p,c in zip(inputs,catalogs)];pairs=[];tentative=[]
    for i in range(len(inputs)):
        for j in range(i+1,len(inputs)):
            # An initial offset from shared identities uses measured HST positions only.
            li={labels[(i,a)]:a for a in range(len(catalogs[i])) if (i,a) in labels}
            lj={labels[(j,b)]:b for b in range(len(catalogs[j])) if (j,b) in labels}
            known=[(a,lj[sid]) for sid,a in li.items() if sid in lj]
            trusted=[(a,b) for a,b in known if stellar((i,a)) and stellar((j,b))]
            if trusted:known=trusted
            guess=np.median([xy[i][a]-xy[j][b] for a,b in known],axis=0) if known else None
            try:
                ri,ii=match(xy[i],xy[j],options,histogram=guess is None,guess=guess)
                err=np.hypot(catalogs[i]['centroid_error_pix'][ri],catalogs[j]['centroid_error_pix'][ii])
                q=quality(xy[i][ri],xy[j][ii],err,options)
                q.update(reference=inputs[i].name,exposure=inputs[j].name,external_identity_seeds=len(known))
                if trusted:
                    guide_error=float(np.median([np.hypot(catalogs[i]['centroid_error_pix'][a],catalogs[j]['centroid_error_pix'][b]) for a,b in trusted]))
                    residual=xy[i][ri]-xy[j][ii]-guess
                    limits=np.minimum(options['match_tolerance_pix'],np.maximum(2*options['scatter_floor_pix'],options['centroid_sigma']*np.hypot(err,guide_error)))
                    allowed=np.linalg.norm(residual,axis=1)<=limits
                    for a,b in zip(ri[allowed],ii[allowed]):tentative.append(((i,int(a)),(j,int(b))))
                    q['HST_guide']=dict(shift_pix=guess.tolist(),uncertainty_pix=guide_error,consistent_candidates=int(allowed.sum()),total_candidates=len(ri),external_coordinates_used=False)
                # No blind peak may contradict measured common external identities.
                if known and 'shift' in q:
                    delta=np.asarray([xy[i][a]-xy[j][b] for a,b in known])-q['shift']
                    if np.max(np.linalg.norm(delta,axis=1))>options['max_current_rms_pix']:
                        q['status']='contradicts_common_HST_identifications'
                if q['status'] in ('validated','refinement_supported'):
                    for a,b,keep in zip(ri,ii,q['fit_inlier_mask']):
                        if keep:join((i,int(a)),(j,int(b)))
                pairs.append(q)
            except Exception as exc:pairs.append(dict(reference=inputs[i].name,exposure=inputs[j].name,status='matching_failed',error=str(exc)))
    # Preserve repeated HST detections even when artifacts spoil the *whole*
    # blind fit. These are identity candidates, not accepted WCS solutions.
    candidate_parent={n:n for n in nodes}
    def candidate_root(n):
        while candidate_parent[n]!=n:n=candidate_parent[n]
        return n
    for u,v in tentative:candidate_parent[candidate_root(v)]=candidate_root(u)
    candidate_groups={}
    for n in nodes:candidate_groups.setdefault(candidate_root(n),[]).append(n)
    repeat_count=0
    selection=association_options or {};max_flux=selection.get('max_flux_ratio',5.);max_width=selection.get('max_width_ratio',2.);max_ellip=selection.get('max_ellipticity_difference',.3)
    for group in candidate_groups.values():
        if len(group)<2 or len({i for i,k in group})!=len(group):continue
        group_set=set(group);edges=[(u,v) for u,v in tentative if u in group_set and v in group_set]
        degrees={n:sum(n in edge for edge in edges) for n in group}
        ids={labels[n] for n in group if n in labels}
        if len(ids)>1:continue
        repeated=len(group)>=3 and min(degrees.values())>=2
        identified_pair=len(group)==2 and all(n in labels for n in group) and len(ids)==1
        if not (repeated or identified_pair):continue
        rows=[catalogs[i][k] for i,k in group]
        widths=np.asarray([r['width_pix'] for r in rows]);flux=np.asarray([r['flux']/fits.getheader(inputs[i])['EXPTIME'] for (i,k),r in zip(group,rows)])
        ellipticity=np.asarray([r['ellipticity'] for r in rows])
        if not np.isfinite(np.r_[widths,flux,ellipticity]).all() or flux.min()<=0 or widths.min()<=0:continue
        if widths.max()/widths.min()>max_width or flux.max()/flux.min()>max_flux or np.ptp(ellipticity)>max_ellip:continue
        for n in group[1:]:join(group[0],n)
        repeat_count+=1
    groups={}
    for n in nodes:groups.setdefault(root(n),[]).append(n)
    selected=[{} for p in inputs];long=[]
    for group in groups.values():
        exposures={i for i,j in group};external_ids={labels[n] for n in group if n in labels}
        ambiguous=len(exposures)!=len(group) or len(external_ids)>1
        common=len(exposures)>=2 and not ambiguous
        token=hashlib.sha256('|'.join(sorted(str(catalogs[i]['hst_source_id'][j]) for i,j in group)).encode()).hexdigest()[:12]
        sid=next(iter(external_ids)) if len(external_ids)==1 else 'HST-'+token
        confidence='catalog-associated HST measurements' if external_ids else 'blind HST pair support; pending image QA'
        for i,j in group:
            c=catalogs[i][j]
            long.append(dict(source_id=sid,hst_source_id=str(c['hst_source_id']),exposure=inputs[i].name,extver=int(c['extver']),x=float(c['x']),y=float(c['y']),
                centroid_error_pix=float(c['centroid_error_pix']),snr=float(c['snr']),width_pix=float(c['width_pix']),external_id=next(iter(external_ids)) if len(external_ids)==1 else '',
                classification='point-like' if c['width_pix']<=2 else 'compact',confidence=confidence,status='common_candidate' if common else ('ambiguous_component' if ambiguous else 'single_exposure')))
            if common:selected[i][j]=sid
    result=[]
    for c,chosen in zip(catalogs,selected):
        idx=sorted(chosen);out=c[idx].copy();out['common_source_id']=[chosen[j] for j in idx]
        out.meta['centroid_quality']=dict(version=QUALITY_VERSION,accepted=len(out),purpose='Native HST centroids; common IDs select pairs, never external-coordinate constraints')
        result.append(out)
    report=dict(status='candidates_established' if any(selected) else 'unresolved',pairs=pairs,matched_nodes=[len(c) for c in result],HST_guide_used=guided,HST_repeat_components=repeat_count,external_positions_used_for_fit=False,
        validation='Common candidates require relative fit, closure and independent no-CR image QA before acceptance')
    return result,Table(rows=long) if long else Table(names=['source_id','exposure','status'],dtype=['U64']*3),report


def resolve_repeated_counterparts(inputs,catalogs,common,associations,options):
    """Identify an unambiguous *repeated HST object* inside a catalog window.

    Missing exposures are not filled with projected positions. Catalog windows
    identify objects only; all returned positions remain measured HST pixels.
    """
    if not len(common):return []
    associations['reason']=np.asarray(associations['reason']).astype('U512')
    associations['centroid_method']=np.asarray(associations['centroid_method']).astype('U64')
    associations['hst_source_id']=np.asarray(associations['hst_source_id']).astype('U128')
    supported={(str(r['exposure']),str(r['hst_source_id'])) for r in common if r['status']=='common_candidate'}
    for row in associations:
        if str(row['status'])=='confirmed_common' and (str(row['exposure']),str(row['hst_source_id'])) not in supported:
            associations['status'][row.index]='valid_not_common'
            associations['reason'][row.index]='Pixel detection retained; HST reference-assisted matching did not corroborate the cross-exposure identity'
    decisions=[];byname={p.name:i for i,p in enumerate(inputs)}
    ids=np.asarray(common['source_id']).astype(str)
    for sid in set(str(r['source_id']) for r in associations):
        group=[r for r in associations if str(r['source_id'])==sid];possible={}
        for row in group:
            i=byname[str(row['exposure'])];v=int(row['extver'])
            if v<0:continue
            with fits.open(inputs[i]) as hd:scale=np.sqrt(abs(np.linalg.det(detector_wcs(hd,v).pixel_scale_matrix)))*3600
            radius=min(options['max_search_radius_pix'],options['search_radius_arcsec']/scale)
            for r in common:
                if r['status']!='common_candidate' or str(r['exposure'])!=inputs[i].name or int(r['extver'])!=v:continue
                if str(r['external_id']) not in ('',sid):continue
                if np.hypot(float(r['x'])-float(row['predicted_x']),float(r['y'])-float(row['predicted_y']))<=radius:possible.setdefault(str(r['source_id']),[]).append((row,r))
        viable={key:rows for key,rows in possible.items() if len(rows)>=3}
        if len(viable)!=1:continue
        key,rows=next(iter(viable.items()));changed=[]
        for row,r in rows:
            if str(row['status'])=='confirmed_common':continue
            i=byname[str(row['exposure'])];c=catalogs[i];idx=np.flatnonzero(np.asarray(c['hst_source_id']).astype(str)==str(r['hst_source_id']))
            if len(idx)!=1:continue
            measured=c[int(idx[0])]
            for name in ('x','y','extver','snr','width_pix','ellipticity','centroid_error_pix'):associations[name][row.index]=measured[name]
            associations['hst_source_id'][row.index]=str(r['hst_source_id'])
            associations['flux_rate'][row.index]=float(measured['flux'])/fits.getheader(inputs[i])['EXPTIME']
            associations['status'][row.index]='confirmed_common'
            associations['centroid_method'][row.index]='HST_repeat_association'
            associations['reason'][row.index]='Repeated native HST detections consistent with measured HST reference shifts; catalog window identifies object only'
            changed.append(inputs[i].name)
        if changed:decisions.append(dict(source_id=sid,HST_component=key,measured_exposures=changed,external_coordinates_used_for_fit=False,missing_exposures_not_fabricated=True))
    return decisions
