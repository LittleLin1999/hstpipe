"""Detector-space evidence for relative matches; no WCS changes or acceptance."""
from pathlib import Path
import json
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from tweakwcs.linearfit import iter_linear_fit
from .sources import project
from .relative_fit import match


def diagnose_pair(reference,image,refcat,cat,options,paths,label):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval
    with fits.open(reference) as hd:tangent=WCS(hd['SCI'].header,hd)
    a=project(refcat,reference,tangent);b=project(cat,image,tangent)
    ri,ii=match(a,b,options,True);delta=a[ri]-b[ii]
    table=Table(dict(candidate=np.arange(len(ri))+1,reference_index=ri,source_index=ii,
        ref_x=refcat['x'][ri],ref_y=refcat['y'][ri],ref_extver=refcat['extver'][ri],
        src_x=cat['x'][ii],src_y=cat['y'][ii],src_extver=cat['extver'][ii],
        tan_x=b[ii,0],tan_y=b[ii,1],dx_pix=delta[:,0],dy_pix=delta[:,1],
        ref_width_pix=refcat['width_pix'][ri],src_width_pix=cat['width_pix'][ii],
        centroid_error_pix=np.hypot(refcat['centroid_error_pix'][ri],cat['centroid_error_pix'][ii]),
        visual_class=np.full(len(ri),'unreviewed',dtype='U32')))
    if len(ri)>=4:
        fit=iter_linear_fit(a[ri],b[ii],wxy=1/np.maximum(table['centroid_error_pix'],1e-6)**2,fitgeom='shift',nclip=3,sigma=(3.,'rmse'))
        center=np.asarray(fit['center'],float);matrix=np.asarray(fit['matrix'],float);shift=np.asarray(fit['shift'],float)
        residual=a[ri]-((b[ii]-center)@matrix.T+center+shift)
        table['residual_x_pix']=residual[:,0];table['residual_y_pix']=residual[:,1];table['geometric_inlier']=fit['fitmask']
    else:
        table['residual_x_pix']=np.full(len(ri),np.nan);table['residual_y_pix']=np.full(len(ri),np.nan);table['geometric_inlier']=np.zeros(len(ri),bool)
    # A histogram of all pair offsets is evidence, not an alternate fitting algorithm.
    all_delta=(a[:,None,:]-b[None,:,:]).reshape(-1,2)
    radius=options['search_radius_pix'];all_delta=all_delta[(abs(all_delta)<radius).all(axis=1)]
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    axes[0,0].hist2d(all_delta[:,0],all_delta[:,1],bins=int(radius*2),range=[[-radius,radius]]*2)
    axes[0,0].set_title('All pair offsets / one-pixel bins')
    axes[0,1].scatter(table['dx_pix'],table['dy_pix']);axes[0,1].set_title('Matched offsets before fit')
    axes[0,2].quiver(table['tan_x'],table['tan_y'],table['residual_x_pix'],table['residual_y_pix'],angles='xy',scale_units='xy',scale=.02)
    axes[0,2].set_title('Residual vectors after shift fit (50x)')
    for ax,coordinate in zip(axes[1,:2],('tan_x','tan_y')):
        ax.scatter(table[coordinate],table['residual_x_pix'],label='dx');ax.scatter(table[coordinate],table['residual_y_pix'],label='dy');ax.legend();ax.set_xlabel(coordinate)
    axes[1,2].scatter(table['centroid_error_pix'],np.hypot(table['residual_x_pix'],table['residual_y_pix']))
    axes[1,2].set_xlabel('Estimated centroid error');axes[1,2].set_ylabel('Radial residual')
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment',label+'_residual_fields','png'),dpi=130);plt.close(fig)
    with fits.open(reference) as ref,fits.open(image) as src:
        for start in range(0,len(table),5):
            fig,axes=plt.subplots(min(5,len(table)-start),4,figsize=(12,3*min(5,len(table)-start)),squeeze=False)
            for row,axesrow in zip(table[start:start+5],axes):
                for frame,prefix,ax1,ax2 in ((ref,'ref',axesrow[0],axesrow[1]),(src,'src',axesrow[2],axesrow[3])):
                    v=int(row[prefix+'_extver']);x=int(round(row[prefix+'_x']));y=int(round(row[prefix+'_y']));r=16
                    h=frame['SCI',v];bounds=(slice(max(0,y-r),min(h.shape[0],y+r+1)),slice(max(0,x-r),min(h.shape[1],x+r+1)))
                    cut=h.data[bounds];dq=frame['DQ',v].data[bounds]
                    lo,hi=ZScaleInterval().get_limits(cut)
                    ax1.imshow(cut,origin='lower',vmin=lo,vmax=hi,cmap='gray');ax1.axvline(x-max(0,x-r),color='red',alpha=.5);ax1.axhline(y-max(0,y-r),color='red',alpha=.5)
                    ax1.set_title(f"#{row['candidate']} {prefix} SCI{v} ({row[prefix+'_x']:.1f},{row[prefix+'_y']:.1f})")
                    ax2.imshow(dq!=0,origin='lower',cmap='gray');ax2.set_title('DQ != 0')
                    table[prefix+'_dq_fraction' if prefix+'_dq_fraction' in table.colnames else prefix+'_dq_fraction']=table[prefix+'_dq_fraction'] if prefix+'_dq_fraction' in table.colnames else np.zeros(len(table))
                    table[prefix+'_dq_fraction'][int(row['candidate'])-1]=np.mean(dq!=0)
            fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment',label+f'_pairs_{start+1:02d}','png'),dpi=120);plt.close(fig)
    table.meta.update(reference=str(reference),image=str(image),label=label,
        caveat='Geometric inliers and visually genuine objects are different classifications')
    output=paths.file('intermediate','alignment',label+'_candidates','ecsv');table.write(output,format='ascii.ecsv',overwrite=True)
    return table


def image_offsets(inputs,paths,options):
    """Exploratory masked correlation on native distortion-corrected grids only."""
    import tempfile
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter,shift as array_shift
    from skimage.registration import phase_cross_correlation
    from .relative_engine import internal_mosaic
    from ...common.target_runner import save
    from ...common.provenance import identity
    results=[];before=[identity(p) for p in inputs]
    with paths.scratch('image-offsets-') as tmp:
        work=Path(tmp);copies=[];products=[]
        for p in inputs:
            dest=work/p.name;paths.copy(p,dest);copies.append(dest)
        reference=internal_mosaic([copies[0]],work/'reference')
        products=[reference]+[internal_mosaic([p],work/p.stem,reference) for p in copies[1:]]
        with fits.open(reference) as ref:
            a=ref['SCI'].data.astype(float);mask_a=(ref['WHT'].data>0)&np.isfinite(a);w=WCS(ref['SCI'].header,ref)
        a=np.where(mask_a,a,0);a=gaussian_filter(a,1.)
        with fits.open(inputs[0]) as hd:
            base_header=hd[0].header.copy();h=hd['SCI'];center=[h.shape[1]/2,h.shape[0]/2]
            fixed_world=WCS(h.header,hd).all_pix2world([center],0);fixed_reference=w.all_world2pix(fixed_world,0)[0]
        for p,product in zip(inputs[1:],products[1:]):
            with fits.open(product) as hd:
                b=hd['SCI'].data.astype(float);mask_b=(hd['WHT'].data>0)&np.isfinite(b)
            b=gaussian_filter(np.where(mask_b,b,0),1.)
            # Clip unmasked extreme events for correlation only; no FITS SCI is edited.
            clip_a=np.percentile(a[mask_a],99.5);clip_b=np.percentile(b[mask_b],99.5)
            aa=np.minimum(a,clip_a);bb=np.minimum(b,clip_b)
            shifts=[]
            for region in (slice(None),slice(0,a.shape[0]//2),slice(a.shape[0]//2,None)):
                try:
                    value=phase_cross_correlation(aa[region],bb[region],reference_mask=mask_a[region],moving_mask=mask_b[region])
                    if isinstance(value,tuple):value=value[0]
                    shifts.append(np.asarray(value,float).tolist())
                except Exception as exc:shifts.append(dict(error=str(exc)))
            with fits.open(p) as hd:
                header=hd[0].header;world=WCS(hd['SCI'].header,hd).all_pix2world([center],0)
                static_shift=fixed_reference-w.all_world2pix(world,0)[0]
                dither=[float(header.get(k,0))-float(base_header.get(k,0)) for k in ('POSTARG1','POSTARG2')]
            result=dict(exposure=p.name,reference=inputs[0].name,correlation_shift_yx_pix=shifts,
                detector_fixed_pattern_prediction_xy_pix=static_shift.tolist(),commanded_delta_postarg_arcsec=dither,
                status='candidate_only',validation='Requires independently verified astronomical matches')
            results.append(result)
            if isinstance(shifts[0],list):
                moved=array_shift(bb,shifts[0],order=1,cval=np.nan)
                fig,axes=plt.subplots(1,2,figsize=(12,5))
                valid=mask_a&mask_b;scale=np.percentile(abs((aa-bb)[valid]),99)
                for ax,array,title in zip(axes,(aa-bb,aa-moved),('Original WCS residual','Exploratory correlation residual')):
                    ax.imshow(array[::4,::4],origin='lower',cmap='RdBu_r',vmin=-scale,vmax=scale);ax.set_title(title)
                fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment','image_offset_'+p.stem,'png'),dpi=110);plt.close(fig)
    if [identity(p) for p in inputs]!=before:raise RuntimeError('Diagnostic correlation modified an input')
    save(paths.file('intermediate','alignment','image_offsets'),results)
    return results


def compare_pair_candidate(reference,image,refcat,catalog,options,paths,label):
    """Exploratory no-CR image QA for a source-supported pair; never commits WCS.

    This is diagnostic evidence, not authorization to combine an unvalidated
    exposure group. All native WCS/drizzle operations use temporary copies.
    """
    import tempfile
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval
    from .relative_engine import internal_mosaic,native_catalog_fit,image_quality
    from ...common.provenance import identity
    from ...common.target_runner import save
    before_inputs=[identity(p) for p in (reference,image)]
    with paths.scratch('pair-quality-') as tmp:
        work=Path(tmp);copies=[]
        for p in (reference,image):
            copy=work/p.name;paths.copy(p,copy);copies.append(copy)
        before=internal_mosaic(copies,work/'before')
        transforms=native_catalog_fit([copies[1]],[catalog],options,reference_image=copies[0],reference_catalog=refcat)
        after=internal_mosaic(copies,work/'after',before)
        qa=image_quality(before,after,refcat,options,paths,label)
        arrays=[];weights=[]
        with fits.open(reference) as hd:ra=hd[0].header['RA_TARG'];dec=hd[0].header['DEC_TARG']
        stamps=[fits.PrimaryHDU()]
        for p,tag in ((before,'BEFORE'),(after,'AFTER')):
            with fits.open(p) as hd:
                w=WCS(hd['SCI'].header,hd);x,y=w.all_world2pix([[ra,dec]],0)[0];x,y=int(round(x)),int(round(y));radius=60
                bounds=(slice(max(0,y-radius),min(hd['SCI'].shape[0],y+radius+1)),slice(max(0,x-radius),min(hd['SCI'].shape[1],x+radius+1)))
                arrays.append(hd['SCI'].data[bounds].copy());weights.append(hd['WHT'].data[bounds].copy())
                for name,array in [('SCI',arrays[-1]),('WHT',weights[-1])]:
                    h=fits.ImageHDU(array,name=tag+'_'+name);h.header['HPPASS']='EXPLORATORY';stamps.append(h)
        fits.HDUList(stamps).writeto(paths.file('intermediate','alignment',label+'_target_cutouts','fits'),overwrite=True,checksum=True)
        fig,axes=plt.subplots(1,4,figsize=(16,4))
        lo,hi=ZScaleInterval().get_limits(arrays[0])
        for ax,array,title in zip(axes,arrays+weights,['Original WCS SCI','Candidate WCS SCI','Original WHT','Candidate WHT']):
            ax.imshow(array,origin='lower',cmap='gray',vmin=lo if 'SCI' in title else None,vmax=hi if 'SCI' in title else None);ax.set_title(title)
        fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment',label+'_target','png'),dpi=130);plt.close(fig)
        report=dict(status='diagnostic_only',image_quality=qa,transformations=transforms,cr_rejection=False,group_validation='not_authorized_by_pair_test',
            caveat='Sources are QC candidates, not confirmed astronomical objects; pair QA cannot establish whole-group registration')
        save(paths.file('intermediate','alignment',label+'_image_quality'),report)
    if [identity(p) for p in (reference,image)]!=before_inputs:raise RuntimeError('Pair diagnostic modified an input')
    return report
