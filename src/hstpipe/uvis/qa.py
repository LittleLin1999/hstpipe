"""Scientific review diagnostics, separate from structural FITS validation."""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from ..common.fits import validate_product
from .fits import image_data

def image_plot(arrays,titles,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import ImageNormalize,AsinhStretch
    path.parent.mkdir(parents=True,exist_ok=True)
    fig,axs=plt.subplots(1,len(arrays),figsize=(4*len(arrays),4),squeeze=False)
    for ax,a,title in zip(axs[0],arrays,titles):
        a=np.asarray(a,dtype=float);a=a[::max(1,a.shape[0]//700),::max(1,a.shape[1]//700)]
        finite=a[np.isfinite(a)]
        lo,hi=np.percentile(finite,[0.5,99.5]) if finite.size else (0,1)
        norm=ImageNormalize(vmin=lo,vmax=hi,stretch=AsinhStretch()) if hi>lo else None
        ax.imshow(a,origin='lower',cmap='gray',norm=norm);ax.set_title(title)
    fig.tight_layout();fig.savefig(path,dpi=120);plt.close(fig)

def target_slice(h,hd,radius=1.0):
    w=WCS(h.header,hd);p=hd[0].header
    try:x,y=w.all_world2pix([[p['RA_TARG'],p['DEC_TARG']]],0)[0]
    except (KeyError,ValueError):return None
    if not np.isfinite([x,y]).all() or not (0<=x<h.data.shape[1] and 0<=y<h.data.shape[0]):return None
    radius=max(1,int(np.ceil(radius/(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600))))
    x,y=int(round(x)),int(round(y))
    return np.s_[max(0,y-radius):min(h.data.shape[0],y+radius+1),max(0,x-radius):min(h.data.shape[1],x+radius+1)]

def morphology(data):
    a=np.array(data,dtype=float);a[~np.isfinite(a)]=0
    sky=float(np.median(a));a=np.maximum(a-sky,0);total=a.sum()
    if total<=0:return dict(flux=0,elongation=None)
    yy,xx=np.indices(a.shape);x=(xx*a).sum()/total;y=(yy*a).sum()/total
    moments=np.array([[((xx-x)**2*a).sum(),((xx-x)*(yy-y)*a).sum()],
                      [((xx-x)*(yy-y)*a).sum(),((yy-y)**2*a).sum()]])/total
    ev=np.linalg.eigvalsh(moments)
    return dict(flux=float(total),centroid=[float(x),float(y)],elongation=float(np.sqrt(ev[1]/max(ev[0],1e-10))),
        caveat='Positive-weight target-cutout moments; extended morphology and CRs can mimic elongation')

def inspect_flc(path,diagnostics=None,*,paths=None):
    from .inventory import validate_flc
    metadata=validate_flc(path);result=[]
    with fits.open(path) as hd:
        for h in hd:
            if h.name!='SCI':continue
            v=h.header.get('EXTVER',1);dq=image_data(hd['DQ',v]);err=image_data(hd['ERR',v])
            good=(dq==0)&np.isfinite(h.data)&np.isfinite(err)&(err>0)
            _,med,std=sigma_clipped_stats(h.data,mask=~good)
            cut=target_slice(h,hd);warnings=list(metadata['warnings'])
            if not good.any():raise ValueError('No valid SCI/ERR/DQ pixels')
            if not np.isfinite(err).all() or np.any(err<0):warnings.append('Invalid ERR pixels')
            from .dq import bit_statistics,masks
            _,_,astrometric_mask=masks(hd,v,{})
            row=dict(dq_bits=bit_statistics(dq),effective_astrometric_mask_fraction=float(astrometric_mask.mean()),extver=v,shape=h.shape,flagged_fraction=float(np.mean(dq!=0)),finite_fraction=float(np.isfinite(h.data).mean()),
                saturation_fraction=float(((dq&256)!=0).mean()),existing_cr_fraction=float(((dq&4096)!=0).mean()),
                hot_pixel_fraction=float(((dq&16)!=0).mean()),sky=float(med),sky_rms=float(std),
                high_outlier_fraction=float((h.data[good]>med+10*max(std,1e-10)).mean()),warnings=warnings)
            if cut is not None:
                row['target']=morphology(h.data[cut]);row['target']['detected']=bool(np.max(h.data[cut])>med+5*std)
                if row['target']['elongation'] and row['target']['elongation']>1.5:warnings.append('Target elongation: inspect morphology and guiding; WCS shifts cannot undo smearing')
            else:warnings.append('Target outside this chip')
            row['artifact_checks']='Conditional visual QA: CTE trails, chip edges, hot pixels, ghosts/scattered light/crosstalk and guiding; no mandatory artifact correction'
            row['cte_background_context']=dict(calibrated_sky_electrons=float(med),postflash_removed=hd[0].header.get('FLSHCORR')=='COMPLETE',flash_duration=hd[0].header.get('FLASHDUR'),flash_level=hd[0].header.get('FLASHLVL'),total_detector_background='Cannot infer pre-calibration sky+dark+postflash solely from calibrated FLC sky')
            if diagnostics:
                arrays=[h.data[::max(1,h.shape[0]//512),::max(1,h.shape[1]//512)],dq[::max(1,h.shape[0]//512),::max(1,h.shape[1]//512)]]
                titles=['SCI including chip boundaries','Original detector DQ']
                if cut is not None:arrays.append(h.data[cut]);titles.append('Target')
                image_plot(arrays,titles,paths.exposure_file('diagnostics','alignment',path.stem.replace('_flc',''),f'sci{v}','png',group=True) if paths else diagnostics/f'{path.stem}_sci{v}.png')
            result.append(row)
    return result

def validate_drizzle(path):
    validate_product(path,('SCI','WHT','CTX'))
    with fits.open(path) as hd:
        sci=hd['SCI'];wht=hd['WHT'].data;ctx=hd['CTX'].data
        if sci.shape!=wht.shape or ctx.shape[-2:]!=sci.shape:raise ValueError('SCI/WHT/CTX geometry mismatch')
        if not np.isfinite(wht).all() or np.any(wht<0) or not np.any(wht>0):raise ValueError('Invalid drizzle coverage')
        if not np.isfinite(sci.data[wht>0]).all():raise ValueError('Nonfinite covered SCI')
        if not WCS(sci.header).has_celestial:raise ValueError('Missing output WCS')
        if sci.header.get('BUNIT','').upper() not in ('ELECTRONS/S','ELECTRONS / S','ELECTRON/S'):raise ValueError('Unexpected output units')

def cr_report(original,working,maskfile,diagnostics,threshold,*,paths=None):
    masks=[fits.PrimaryHDU()];report=[]
    with fits.open(original) as raw,fits.open(working) as work:
        for h in raw:
            if h.name!='SCI':continue
            v=h.header.get('EXTVER',1);before=image_data(raw['DQ',v]);after=image_data(work['DQ',v])
            required=before&np.uint16(65535^4096)
            if np.any((after&required)!=required):raise ValueError('Detector/calwf3 DQ bits were reset beyond the explicit AstroDrizzle 4096 reassessment')
            inherited=(before&4096)!=0;current=(after&4096)!=0
            new=current&~inherited
            # Retain the newly validated current mask; raw inherited flags remain untouched.
            mh=fits.ImageHDU(current.astype('u1'),name='SCI');mh.header['EXTVER']=v;mh.header['BUNIT']='UNITLESS';mh.header['MASKONE']='CURRENT_DRIZ_CR';masks.append(mh)
            cut=target_slice(h,raw);row=dict(extver=v,inherited_4096_fraction=float(inherited.mean()),calwf3_8192_fraction=float(((before&8192)!=0).mean()),reconfirmed_4096_fraction=float((current&inherited).mean()),cleared_inherited_4096_fraction=float((inherited&~current).mean()),current_cr_fraction=float(current.mean()),new_cr_fraction=float(new.mean()),target_cr_fraction=None,warnings=[])
            if cut is not None:
                row['target_cr_fraction']=float(current[cut].mean())
                row['new_target_cr_fraction']=float(new[cut].mean())
                _,med,std=sigma_clipped_stats(h.data[cut]);bright=h.data[cut]>med+5*std
                core=float(current[cut][bright].mean()) if bright.any() else 0.0
                row['bright_source_cr_fraction']=core
                row['raw_bright_pixel_cr_fraction']=core
                row['raw_bright_pixel_metric_caveat']='Raw high pixels and whole-aperture area include cosmic rays and sky; persistent-source QA determines real-emission loss'
                image_plot([h.data[cut],current[cut]],['Original target SCI','Current AstroDrizzle CR flags (1=rejected)'],paths.exposure_file('diagnostics','drizzle',original.stem.replace('_flc',''),f'cr{v}','png',group=True) if paths else diagnostics/f'{original.stem}_cr{v}.png')
            from scipy.ndimage import maximum_filter,binary_dilation
            _,level,noise=sigma_clipped_stats(h.data)
            peaks=(h.data==maximum_filter(h.data,size=7))&(h.data>level+20*max(noise,1e-10))
            bright_regions=binary_dilation(peaks,iterations=3)
            row['bright_source_region_cr_fraction']=float(new[bright_regions].mean()) if bright_regions.any() else None
            row['source_validation_required']=True
            if row['new_cr_fraction']>0.1:row['warnings'].append('Excessive rejection; inspect alignment')
            report.append(row)
    # A fresh reduction has no drizzle directory until its first retained output.
    from pathlib import Path
    Path(maskfile).parent.mkdir(parents=True,exist_ok=True)
    fits.HDUList(masks).writeto(maskfile,overwrite=True,checksum=True)
    return report

def final_report(product,inputs,diagnostics,*,paths=None):
    report={};cutouts=[];titles=[];fluxes=[]
    for p in inputs:
        with fits.open(p) as hd:
            for h in hd:
                if h.name!='SCI':continue
                cut=target_slice(h,hd)
                if cut is not None:
                    a=h.data[cut]/hd[0].header['EXPTIME'];cutouts.append(a.copy());titles.append(p.stem)
                    fluxes.append(morphology(a)['flux'])
    with fits.open(product) as hd:
        sci=hd['SCI'].data;wht=hd['WHT'].data;good=wht>0
        _,med,std=sigma_clipped_stats(sci,mask=~good)
        report.update(coverage_fraction=float(good.mean()),sky=float(med),sky_rms=float(std),weight_min=float(wht[good].min()),weight_max=float(wht.max()))
        ctx=hd['CTX'].data
        planes=ctx[None,...] if ctx.ndim==2 else ctx
        contributions=np.zeros(sci.shape,dtype='u2')
        for plane in planes:
            for bit in range(32):contributions+=((plane.astype('uint32')>>bit)&1).astype('u2')
        report['contribution_histogram']={str(k):int((contributions==k).sum()) for k in np.unique(contributions)}
        tiles=[]
        for iy in range(0,sci.shape[0],128):
            for ix in range(0,sci.shape[1],128):
                a=sci[iy:iy+128,ix:ix+128];g=good[iy:iy+128,ix:ix+128]
                if g.sum()>100:tiles.append(float(sigma_clipped_stats(a,mask=~g)[1]))
        report['sky_tile_medians']=tiles
        report['sky_structure_peak_to_peak']=float(np.ptp(tiles)) if tiles else None
        # Source-masked adjacent residual correlation is an indicator, not a covariance correction.
        blank=good&(abs(sci-med)<3*max(std,1e-15));valid=blank[:,1:]&blank[:,:-1]
        report['adjacent_pixel_correlation']=float(np.corrcoef(sci[:,1:][valid],sci[:,:-1][valid])[0,1]) if valid.sum()>30 else None
        corr=report['adjacent_pixel_correlation']
        if corr is not None and not np.isfinite(corr):report['adjacent_pixel_correlation']=None
        cut=target_slice(hd['SCI'],hd)
        if cut is not None:
            a=sci[cut];cutouts.extend([a.copy(),wht[cut].copy()]);titles.extend(['Final SCI electrons/s','Final WHT'])
            report['target']=morphology(a);report['target']['input_fluxes']=fluxes
            center=wht[cut];report['target_weight_low_fraction']=float((center<0.5*np.median(wht[good])).mean())
            report['photometry_caveat']='Cutout sums are qualitative: native detector pixels have varying area, PSFs and coverage; no PAM or aperture correction applied'
            _,tsky,tnoise=sigma_clipped_stats(a)
            core=a>tsky+5*max(tnoise,1e-15)
            local=center[(~core)&(center>0)]
            ratio=float(np.median(center[core])/np.median(local)) if core.any() and local.size else None
            report['target_core_weight_ratio']=ratio
            report['warnings']=['Target-shaped WHT deficit; inspect rejected real emission'] if ratio is not None and ratio<0.8 else []
        image_plot([sci,wht,contributions],['Final SCI','WHT','Input contributions'],paths.file('diagnostics','drizzle','coverage','png') if paths else diagnostics/'coverage.png')
    if cutouts:image_plot(cutouts,titles,paths.file('diagnostics','drizzle','target_comparison','png') if paths else diagnostics/'target_comparison.png')
    report['validation_status']='review_required'
    report['limitations']=['CTE residuals and guiding smearing need visual review','Correlated noise measured qualitatively; no ERR or covariance map','Corner footprint approximation used; inspect detailed CTX coverage','No external photometric or official MAST comparison performed']
    return report
