"""External-reference association and validation, independent of FLC relatives."""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.stats import mad_std

def reference_audit(image,sources,reference,options):
    """Audit every eligible row, including cone objects outside actual coverage."""
    with fits.open(image) as hd:
        w=WCS(hd['SCI'].header,hd);xy=w.all_world2pix(np.column_stack([reference['ra'],reference['dec']]),0)
        shape=hd['SCI'].shape;weight=hd['WHT'].data if 'WHT' in hd else np.ones(shape)
        covered=[]
        for x,y in xy:
            ix=int(round(x)) if np.isfinite(x) else -1;iy=int(round(y)) if np.isfinite(y) else -1
            covered.append(0<=ix<shape[1] and 0<=iy<shape[0] and weight[iy,ix]>0)
    table=reference.copy();table['predicted_x']=xy[:,0];table['predicted_y']=xy[:,1];table['covered']=covered
    table['association_status']=np.where(covered,'covered_undetected','outside_coverage').astype('U64')
    table['candidate_count']=np.zeros(len(table),int);table['hst_source_index']=np.full(len(table),-1)
    if len(sources) and len(reference):
        src=SkyCoord(np.asarray(sources['ra'])*u.deg,np.asarray(sources['dec'])*u.deg)
        for i,row in enumerate(reference):
            distance=SkyCoord(float(row['ra'])*u.deg,float(row['dec'])*u.deg).separation(src).arcsec
            near=np.flatnonzero(distance<=options['search_radius_arcsec']);table['candidate_count'][i]=len(near)
            if len(near)==1 and covered[i]:table['association_status'][i]='candidate';table['hst_source_index'][i]=near[0]
            elif len(near)>1 and covered[i]:table['association_status'][i]='ambiguous_candidates'
    return table


def plot_reference_audit(image,sources,reference,audit,matched,paths,label):
    """Small native-pixel counterpart/coverage diagnostic for the offline report."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm
    indices=np.unique(np.asarray(matched['reference_index'],int))[:12] if len(matched) else np.flatnonzero(audit['covered'])[:12]
    if not len(indices):return
    paths.area('diagnostics',create=True)
    fig,axes=plt.subplots(len(indices),2,figsize=(8,max(2.4,2.4*len(indices))),squeeze=False)
    with fits.open(image) as hd:
        data=hd['SCI'].data;weight=hd['WHT'].data if 'WHT' in hd else np.ones(data.shape)
        for panels,i in zip(axes,indices):
            x,y=audit['predicted_x'][i],audit['predicted_y'][i];radius=20;ix,iy=int(round(x)),int(round(y))
            sx=slice(max(0,ix-radius),min(data.shape[1],ix+radius+1));sy=slice(max(0,iy-radius),min(data.shape[0],iy+radius+1))
            rows=matched[np.asarray(matched['reference_index'])==i]
            for ax,array,title in zip(panels,(data,weight),('SCI','WHT')):
                stamp=array[sy,sx];finite=stamp[np.isfinite(stamp)]
                if not len(finite):continue
                ax.imshow(stamp,origin='lower',extent=(sx.start,sx.stop,sy.start,sy.stop),cmap='gray',norm=simple_norm(finite,'asinh',percent=99))
                ax.scatter([x],[y],marker='+',color='orange',label='catalog prediction')
                for row in rows:ax.scatter(row['x'],row['y'],facecolors='none',edgecolors='lime' if row['accepted'] else 'red',s=80)
                ax.set_xlim(sx.start,sx.stop);ax.set_ylim(sy.start,sy.stop)
                ax.set_title(str(reference['source_id'][i])+' '+title+' '+str(audit['association_status'][i]),fontsize=8)
    fig.suptitle(label+': orange prediction; green coherent HST association; red rejected',fontsize=9)
    fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment','absolute_references_'+label,'png'),dpi=95);plt.close(fig)


def absolute_matches(sources,reference,options):
    """STScI histogram candidates, uniqueness checks, robust offset rejection.

    All coordinates here are external-stage measurements only. This function is
    never used to introduce catalog constraints in the FLC relative solution.
    """
    names=['source_index','reference_index','x','y','ra','dec','ref_ra','ref_dec','dx_arcsec','dy_arcsec','accepted']
    if not len(sources) or not len(reference):return Table(names=names,dtype=['i8','i8']+['f8']*8+['bool'])
    sky=SkyCoord(np.asarray(sources['ra'])*u.deg,np.asarray(sources['dec'])*u.deg)
    ref=SkyCoord(np.asarray(reference['ra'])*u.deg,np.asarray(reference['dec'])*u.deg)
    index,distance,_=sky.match_to_catalog_sky(ref);reverse,_,_=ref.match_to_catalog_sky(sky)
    good=(distance.arcsec<=options['search_radius_arcsec'])&(reverse[index]==np.arange(len(sky)))
    if len(ref)>1:
        _,second,_=sky.match_to_catalog_sky(ref,nthneighbor=2)
        good&=second.arcsec>np.maximum(options.get('ambiguity_arcsec',.2),distance.arcsec*options.get('ambiguity_ratio',2))
    selected=np.flatnonzero(good);ri=index[good];method='mutual_unique_WCS'
    # Use the maintained TweakReg matcher to recover coherent nonzero offsets.
    # Keep the initial WCS associations for the audit, including rejected rows.
    histogram_pairs=set()
    if min(len(sky),len(ref))>=3:
        from tweakwcs.matchutils import XYXYMatch
        origin=ref[0];sx,sy=origin.spherical_offsets_to(sky);rx,ry=origin.spherical_offsets_to(ref)
        sc=Table(dict(TPx=sx.arcsec,TPy=sy.arcsec));rc=Table(dict(TPx=rx.arcsec,TPy=ry.arcsec))
        try:
            ir,ii=XYXYMatch(searchrad=options['search_radius_arcsec'],separation=options.get('ambiguity_arcsec',.2),
                tolerance=options.get('match_tolerance_arcsec',.2),use2dhist=True)(rc,sc,tp_pscale=.04,tp_units='arcsec')
            histogram_pairs=set(zip(map(int,ii),map(int,ir)));method='STScI_XYXYMatch_2D_histogram'
            existing=set(zip(map(int,selected),map(int,ri)));extra=sorted(histogram_pairs-existing)
            if extra:selected=np.r_[selected,[p[0] for p in extra]];ri=np.r_[ri,[p[1] for p in extra]]
        except (ValueError,RuntimeError):pass
    dx,dy=sky[selected].spherical_offsets_to(ref[ri]);dx=dx.arcsec;dy=dy.arcsec
    # Median/MAD, rather than ordinary standard deviation (which large false
    # associations inflate). Never clip below catalog+centroid measurement noise.
    errors=np.full(len(selected),options.get('systematic_floor_arcsec',.01))
    for axis in ('ra','dec'):
        if axis+'_error_mas' in reference.colnames:errors=np.hypot(errors,np.asarray(reference[axis+'_error_mas'],float)[ri]/1000)
    if 'centroid_error_pix' in sources.colnames:
        errors=np.hypot(errors,np.asarray(sources['centroid_error_pix'],float)[selected]*options.get('pixel_scale_arcsec',.04))
    keep=np.isfinite(dx+dy+errors)
    for _ in range(5):
        if keep.sum()<2:break
        center=np.array([np.median(dx[keep]),np.median(dy[keep])]);scatter=np.hypot(mad_std(dx[keep]),mad_std(dy[keep]))
        limit=options.get('clip_sigma',3)*np.maximum(scatter,errors)
        updated=keep&(np.hypot(dx-center[0],dy-center[1])<=np.maximum(limit,options.get('numeric_floor_arcsec',1e-8)))
        if np.array_equal(updated,keep):break
        keep=updated
    # Reject duplicate source/reference identities; a histogram alternative must
    # not silently create two measurements of the same astronomical object.
    for values in (selected,ri):
        for value in np.unique(values[keep]):
            group=np.flatnonzero(keep&(values==value))
            if len(group)>1:keep[group]=False
    table=Table([selected,ri,np.asarray(sources['x'])[selected],np.asarray(sources['y'])[selected],
        np.asarray(sources['ra'])[selected],np.asarray(sources['dec'])[selected],np.asarray(reference['ra'])[ri],np.asarray(reference['dec'])[ri],dx,dy,keep],names=names)
    table['positional_error_arcsec']=errors
    table['rejection_reason']=np.where(keep,'','incoherent_offset_or_ambiguous_identity')
    table['reference_id']=np.asarray(reference['source_id']).astype(str)[ri] if 'source_id' in reference.colnames else ri.astype(str)
    table.meta['matching_method']=method
    return table


def validate_absolute(sources,reference,options,shape=None):
    table=absolute_matches(sources,reference,options);kept=table[table['accepted']]
    result=dict(status='insufficient_verified_references',matched=len(kept),fit_inliers=len(kept),initial_associations=len(table),
                verified_associations=len(kept),independent_validation_sources=0,rejected=len(table)-len(kept),detected=len(sources),reference_count=len(reference),
                all_associations_initial_rms_arcsec=float(np.sqrt(np.mean(np.asarray(table['dx_arcsec'])**2+np.asarray(table['dy_arcsec'])**2))) if len(table) else None,
                thresholds=dict(options),source_quality=sources.meta.get('centroid_quality'),source_selection=sources.meta.get('compact_selection'),rms_arcsec=None,robust_scatter_arcsec=None,offset_arcsec=None,uncertainty_arcsec=None)
    if not len(kept):return result,table
    from tweakwcs.linearfit import iter_linear_fit
    offsets=np.column_stack([kept['dx_arcsec'],kept['dy_arcsec']]);error=np.asarray(kept['positional_error_arcsec'])
    fit=iter_linear_fit(offsets,np.zeros_like(offsets),fitgeom='shift',wxy=1/np.maximum(error,1e-8)**2,nclip=0)
    offset=np.asarray(fit['shift'],float);res=offsets-offset
    rms=float(np.sqrt(np.mean(np.sum(offsets**2,axis=1))));postfit=float(np.sqrt(np.mean(np.sum(res**2,axis=1))))
    scatter=float(np.hypot(mad_std(res[:,0]),mad_std(res[:,1])))
    span=float(np.hypot(np.ptp(kept['x']),np.ptp(kept['y']))/np.hypot(*shape)) if shape is not None else None
    loo=[]
    if len(kept)>1:
        for i in range(len(kept)):
            use=np.arange(len(kept))!=i;weights=1/np.maximum(error[use],1e-8)**2
            center=np.average(offsets[use],axis=0,weights=weights);loo.append(float(np.linalg.norm(offsets[i]-center)))
    uncertainty=float(np.hypot(max(scatter,postfit)/np.sqrt(len(kept)),options['systematic_floor_arcsec']))
    result.update(rms_arcsec=rms,initial_offset_rms_arcsec=rms,postfit_rms_arcsec=postfit,robust_scatter_arcsec=scatter,
        fitted_shift_arcsec=offset.tolist(),offset_arcsec=offset.tolist(),span_fraction=span,uncertainty_arcsec=uncertainty,
        shift_uncertainty_arcsec=np.sqrt(1/np.sum(1/np.maximum(error,1e-8)**2)).tolist(),
        leave_one_out_rms_arcsec=float(np.sqrt(np.mean(np.square(loo)))) if loo else None,
        leave_one_out_max_arcsec=max(loo) if loo else None,
        matching_method=table.meta['matching_method'],
        uncertainty_caveat='Leave-one-out prediction uses alternate fitted subsets, not independent withheld objects; catalog systematics and morphology remain')
    minimum=max(options.get('min_matches',1),1 if options.get('fitgeometry','shift')=='shift' else 2)
    if len(kept)<minimum:result['status']='software_match_count_restriction';return result,table
    if len(kept)==1:result['status']='diagnostic_single_reference';return result,table
    passed=(rms<=options['max_rms_arcsec'] and np.linalg.norm(offset)<=options['max_offset_arcsec'] and
            (not options.get('require_leave_one_out',True) or (loo and np.sqrt(np.mean(np.square(loo)))<=options['max_rms_arcsec'])) and
            (span is None or span>=options.get('min_span_fraction',.1)))
    result['status']='validated' if passed else 'excessive_postfit_residuals_or_distribution'
    return result,table

