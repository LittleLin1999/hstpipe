"""Detector-independent, mask-aware centroid QC. Never edits science pixels.

A centroid need not belong to its segmentation footprint. Missing essential
core measurements, rather than footprint membership, determine rejection.
"""
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import binary_dilation
from ..fits import image_data

QUALITY_VERSION='uvis-bitmask-3'
QUALITY_DEFAULTS=dict(core_radius_pix=2.5,quality_radius_pix=7,min_core_valid_fraction=.75,min_sector_valid_fraction=.25,max_blend_core_fraction=.3,min_sigma_pix=.65,max_sigma_pix=8.,max_axis_ratio=2.5,max_centroid_error_pix=1.,stability_floor_pix=.5,point_error_floor=.08,galaxy_error_floor=.2,detection_ignore_bits=4096,centroid_ignore_bits=4096,grow_dq_bits=8192,dq_grow_radius=1,max_diagnostic_cutouts=12)


def pixel_masks(hd,version,options,purpose='centroid'):
    from ..dq import masks
    return masks(hd,version,options,purpose)


def _number(row,name,default=np.nan):
    try:
        value=row[name]
        return default if np.ma.is_masked(value) else float(value)
    except (KeyError,ValueError,TypeError):return default


def _segment_patch(seg,sy,sx):
    if not isinstance(seg,dict):return seg[sy,sx]
    out=np.zeros((sy.stop-sy.start,sx.stop-sx.start),dtype=np.int32)
    data=seg['data'];x0=seg['x0'];y0=seg['y0']
    xa=max(sx.start,x0);xb=min(sx.stop,x0+data.shape[1]);ya=max(sy.start,y0);yb=min(sy.stop,y0+data.shape[0])
    if xa<xb and ya<yb:out[ya-sy.start:yb-sy.start,xa-sx.start:xb-sx.start]=data[ya-y0:yb-y0,xa-x0:xb-x0]
    return out


def validate_centroids(image,catalog,options,*,method='custom',segmentations=None,paths=None,label=None,return_all=False,contexts=None):
    """Validate any detector's zero-origin catalog; return accepted rows and audit.

    Optional segmentation maps describe morphology, never dictate centroid
    footprint membership. Detector DQ, coverage and dilation are recorded
    separately. Existing 4096 bits are recorded as input flags, not newly found CRs.
    No interpolation or automatic centroid replacement is performed.
    """
    options=dict(QUALITY_DEFAULTS,**options)
    audit=Table(catalog,copy=True);n=len(audit)
    for key in ('x','y'):
        if key not in audit.colnames:audit[key]=np.full(n,np.nan)
    if 'extver' not in audit.colnames:
        with fits.open(image) as hd:first=hd['SCI'].header.get('EXTVER',1)
        audit['extver']=np.full(n,first)
    metrics=[];cutouts=[]
    with fits.open(image) as hd:
        root=hd[0].header.get('ROOTNAME',Path(image).stem);contexts=dict(contexts or {})
        for v in np.unique(audit['extver']):
            v=int(v)
            if v in contexts:continue
            h=hd['SCI',v];a=np.asarray(h.data,float);dq,base,mask=pixel_masks(hd,v,options)
            _,sky,rms=sigma_clipped_stats(a,mask=mask)
            error=np.asarray(image_data(hd['ERR',v]),float) if ('ERR',v) in hd else np.full(a.shape,rms)
            contexts[v]=dict(data=a,dq=dq,base=base,mask=mask,sky=float(sky),rms=float(rms),error=error)
        for index,row in enumerate(audit):
            x=_number(row,'x');y=_number(row,'y');v=int(row['extver']);ctx=contexts[v];a=ctx['data']
            record=dict(source_id=str(row['source_id']) if 'source_id' in row.colnames else f'{root}-{method}-{v}-{index+1}',
                exposure_rootname=str(root),ccdchip=int(hd['SCI',v].header.get('CCDCHIP',-1)),detection_method=str(row['detection_method']) if 'detection_method' in row.colnames else method,
                snr=_number(row,'snr'),width_pix=_number(row,'width_pix'),minor_width_pix=_number(row,'minor_width_pix'),
                centroid_error_pix=_number(row,'centroid_error_pix'),ellipticity=np.nan,fwhm_pix=np.nan,
                local_masked_fraction=np.nan,core_valid_fraction=np.nan,input_core_valid_fraction=np.nan,min_sector_valid_fraction=np.nan,input_dq_fraction=np.nan,
                input_cr_fraction=np.nan,recovered_from_inherited_cr=False,processing_mask_fraction=np.nan,centroid_stability_pix=np.nan,centroid_in_segment=False,core_dq_bits=0)
            flags=[];seg=(segmentations or {}).get((v,int(_number(row,'segmentation_key',-1))),(segmentations or {}).get(v));label_id=int(_number(row,'segment_label',0))
            if not np.isfinite([x,y]).all():flags.append('nonfinite_centroid')
            elif x<0 or y<0 or x>=a.shape[1] or y>=a.shape[0]:flags.append('centroid_outside_image')
            else:
                ix=int(np.clip(round(x),0,a.shape[1]-1));iy=int(np.clip(round(y),0,a.shape[0]-1));radius=int(options.get('quality_radius_pix',7))
                sy=slice(max(0,iy-radius),min(a.shape[0],iy+radius+1));sx=slice(max(0,ix-radius),min(a.shape[1],ix+radius+1))
                data=a[sy,sx];dq=ctx['dq'][sy,sx];mask=ctx['mask'][sy,sx];base=ctx['base'][sy,sx];error=ctx['error'][sy,sx]
                yy,xx=np.indices(data.shape);xx=xx+sx.start-x;yy=yy+sy.start-y;r2=xx**2+yy**2
                core=r2<=options.get('core_radius_pix',2.5)**2;local=r2<=radius**2
                valid=~mask;core_valid=valid&core
                ideal_y,ideal_x=np.mgrid[-radius:radius+1,-radius:radius+1]
                expected_core=((ideal_x-(x-ix))**2+(ideal_y-(y-iy))**2<=options.get('core_radius_pix',2.5)**2).sum()
                coverage=float(core_valid.sum()/max(1,expected_core))
                sectors=[]
                for quadrant in range(4):
                    sector=core&(((xx>=0).astype(int)+2*(yy>=0).astype(int))==quadrant)
                    sectors.append(float((sector&valid).sum()/max(1,sector.sum())))
                record.update(local_masked_fraction=float(np.mean(mask[local])),core_valid_fraction=coverage,input_core_valid_fraction=float((core&~base).sum()/max(1,expected_core)),min_sector_valid_fraction=min(sectors),
                    input_dq_fraction=float(np.mean(dq[local]!=0)),input_cr_fraction=float(np.mean((dq[local]&4096)!=0)),
                    processing_mask_fraction=float(np.mean((mask&~(dq!=0))[local])),core_dq_bits=int(np.bitwise_or.reduce(dq[core],initial=0)))
                sky=_number(row,'local_sky',ctx['sky'])
                signal=np.where(valid&local,np.maximum(data-sky,0),0)
                if not np.isfinite(record['snr']):
                    good=valid&local;noise=np.sqrt(np.sum(error[good]**2));record['snr']=float(np.sum((data-sky)[good])/max(noise,1e-12))
                if not np.isfinite(record['width_pix']) or not np.isfinite(record['minor_width_pix']):
                    # Raw SCI core moments identify spikes; no smoothed image enters QC.
                    moment=signal*(r2<=3.5**2);total=moment.sum()
                    if total>0:
                        xc=(moment*xx).sum()/total;yc=(moment*yy).sum()/total
                        cov=np.array([[(moment*(xx-xc)**2).sum(),(moment*(xx-xc)*(yy-yc)).sum()],[(moment*(xx-xc)*(yy-yc)).sum(),(moment*(yy-yc)**2).sum()]])/total
                        eigen=np.maximum(np.linalg.eigvalsh(cov),0);record['minor_width_pix']=float(np.sqrt(eigen[0]));record['width_pix']=float(np.sqrt(eigen[1]))
                major=record['width_pix'];minor=record['minor_width_pix'];snr=record['snr']
                record['ellipticity']=float(1-minor/major) if major>0 else np.nan
                record['fwhm_pix']=float(2.355*np.sqrt(major*minor)) if major>=0 and minor>=0 else np.nan
                floor=options.get('galaxy_error_floor',.2) if major>2 else options.get('point_error_floor',.08)
                error_estimate=max(major/max(snr,1e-12),floor) if np.isfinite(major+snr) else np.inf
                record['centroid_error_pix']=max(record['centroid_error_pix'],error_estimate) if np.isfinite(record['centroid_error_pix']) else error_estimate
                centers=[]
                for aperture_radius in (3.5,5.5):
                    weights=signal*(r2<=aperture_radius**2);total=weights.sum()
                    centers.append(np.array([(weights*xx).sum(),(weights*yy).sum()])/total if total>0 else np.full(2,np.nan))
                record['centroid_stability_pix']=float(np.linalg.norm(centers[0]-centers[1]))
                center_masked=bool(ctx['mask'][iy,ix]);center_invalid=bool(ctx['base'][iy,ix]);bits=int(ctx['dq'][iy,ix])
                if bits&(256|2048):flags.append('saturated_or_invalid_core')
                elif center_invalid and not bits:flags.append('saturated_or_invalid_core')
                if center_masked:flags.append('masked_core_centroid')
                if coverage<options.get('min_core_valid_fraction',.75):flags.append('insufficient_core_coverage')
                if min(sectors)<options.get('min_sector_valid_fraction',.25) and coverage<.9:flags.append('asymmetric_core_mask')
                if seg is not None:
                    patch=_segment_patch(seg,sy,sx)
                    own=(patch==label_id) if label_id else patch>0
                    center_label=patch[iy-sy.start,ix-sx.start]
                    record['centroid_in_segment']=bool(center_label==label_id) if label_id else bool(center_label>0)
                    ring=(own&local&~core).sum()>=3 and not (own&core).any()
                    if (center_masked and (own&local).any()) or (ring and coverage<.9):flags.append('ring_halo_unreliable_center')
                    measurement=r2<=min(radius,max(3.5,2*major))**2
                    other=(patch>0)&~own&measurement
                    total=np.sum(signal[measurement]);contamination=float(np.sum(signal[other])/max(total,1e-12))
                    if contamination>options.get('max_blend_core_fraction',.3):flags.append('severe_blending')
                if 'blended' in row.colnames and bool(row['blended']):flags.append('severe_blending')
                if not np.isfinite(major+minor) or major>options.get('max_sigma_pix',8) or major/max(minor,1e-9)>options.get('max_axis_ratio',2.5):flags.append('unreliable_morphology')
                if minor<options.get('min_sigma_pix',.65):flags.append('cosmic_ray_like')
                if not np.isfinite(snr) or snr<options.get('snr',5):flags.append('low_snr')
                if not np.isfinite(record['centroid_error_pix']) or record['centroid_error_pix']>options.get('max_centroid_error_pix',1.):flags.append('poor_centroid_precision')
                if major<=2 and record['centroid_stability_pix']>max(options.get('stability_floor_pix',.5),3*record['centroid_error_pix']):flags.append('unstable_centroid')
                if any(f in flags for f in ('masked_core_centroid','ring_halo_unreliable_center','saturated_or_invalid_core')):
                    cutouts.append((index,data.copy(),dq.copy(),mask.copy(),_segment_patch(seg,sy,sx).copy() if seg is not None else np.zeros(data.shape,'i4'),x-sx.start,y-sy.start))
            if not flags and np.isfinite([x,y]).all() and not hd[0].header.get('HPCRRUN',False):
                previous=binary_dilation(dq!=0,iterations=max(1,int(options.get('edge_pixels',2))))
                other=binary_dilation((dq&np.uint16(65535^4096))!=0,iterations=max(1,int(options.get('edge_pixels',2))))
                yy0,xx0=np.indices(previous.shape);core0=(xx0-(x-sx.start))**2+(yy0-(y-sy.start))**2<=options['core_radius_pix']**2
                if core0.any():
                    old_bad=bool(previous[iy-sy.start,ix-sx.start]) or np.mean(~previous[core0])<options['min_core_valid_fraction']
                    other_bad=bool(other[iy-sy.start,ix-sx.start]) or np.mean(~other[core0])<options['min_core_valid_fraction']
                    record['recovered_from_inherited_cr']=bool(old_bad and not other_bad and np.any(dq&4096))
            record.update(centroid_valid=not flags,centroid_status='accepted' if not flags else 'rejected',quality_flags=';'.join(dict.fromkeys(flags)),rejection_reason=flags[0] if flags else '')
            metrics.append(record)
    fields=list(metrics[0]) if metrics else ['source_id','exposure_rootname','detection_method','snr','width_pix','minor_width_pix','centroid_error_pix','fwhm_pix','ellipticity','local_masked_fraction','core_valid_fraction','input_dq_fraction','input_cr_fraction','processing_mask_fraction','centroid_valid','centroid_status','quality_flags','rejection_reason']
    for key in fields:
        if metrics:audit[key]=[r[key] for r in metrics]
        else:audit[key]=np.array([],dtype=bool if key=='centroid_valid' else ('U64' if key in ('source_id','exposure_rootname','detection_method','centroid_status','quality_flags','rejection_reason') else float))
    selected=np.asarray(audit['centroid_valid'],bool).copy()
    for v in np.unique(audit['extver']):
        idx=np.flatnonzero(selected&(audit['extver']==v))
        if 'flux' in audit.colnames:idx=idx[np.argsort(np.nan_to_num(audit['flux'][idx],nan=-np.inf))[::-1]]
        selected[idx[int(options.get('max_sources',1000)):]]=False
    audit['selected_for_astrometry']=selected
    accepted=audit[selected]
    summary=dict(version=QUALITY_VERSION,detected=n,accepted=len(accepted),rejected=int((~np.asarray(audit['centroid_valid'],bool)).sum()),selection_limited=int((np.asarray(audit['centroid_valid'],bool)&~selected).sum()),rejection_counts={f:sum(f in r['quality_flags'].split(';') for r in metrics) for f in sorted({f for r in metrics for f in r['quality_flags'].split(';') if f})})
    audit.meta.update(image=Path(image).name,quality_parameters=dict(options),centroid_quality=summary,dq_provenance='Input-file DQ preserved; input 4096 is not classified as a newly introduced CR flag',core_policy='Mask topology/coverage; centroid membership in segmentation is not required')
    if paths:
        paths.area('intermediate','alignment',create=True);paths.area('diagnostics','alignment',create=True)
        label=label or Path(image).stem+'_'+method
        output=paths.file('intermediate','alignment',label+'_quality','ecsv');audit.write(output,format='ascii.ecsv',overwrite=True)
        accepted.meta.update(centroid_quality=summary,quality_file=str(output))
        _diagnostics(image,audit,cutouts,paths,label,options)
    else:
        accepted.meta.update(centroid_quality=summary,rejected_sources=[{k:r[k] for k in ('source_id','centroid_status','quality_flags','rejection_reason')} for r in metrics if not r['centroid_valid']])
    return (accepted,audit) if return_all else accepted


def _diagnostics(image,audit,cutouts,paths,label,options):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval
    if options.get('diagnostic_overlay',True):
        with fits.open(image) as hd:
            chips=[h for h in hd if h.name=='SCI'];fig,axes=plt.subplots(1,len(chips),figsize=(6*len(chips),5),squeeze=False)
            for h,ax in zip(chips,axes[0]):
                stride=max(1,h.shape[0]//512);sample=h.data[::stride,::stride];finite=sample[np.isfinite(sample)]
                lo,hi=ZScaleInterval().get_limits(finite) if finite.size else (0,1)
                ax.imshow(sample,origin='lower',extent=(0,h.shape[1],0,h.shape[0]),vmin=lo,vmax=hi,cmap='gray')
                rows=audit[audit['extver']==h.header.get('EXTVER',1)]
                for good,color,marker in ((True,'lime','o'),(False,'red','x')):
                    use=rows[np.asarray(rows['centroid_valid'])==good];ax.scatter(use['x'],use['y'],s=25,c=color,marker=marker,alpha=.65,label='accepted' if good else 'rejected')
                ax.legend();ax.set_title(f"{Path(image).name} SCI{h.header.get('EXTVER',1)}")
            fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment',label+'_quality','png'),dpi=110);plt.close(fig)
    if not cutouts:return
    stamps=[fits.PrimaryHDU()]
    for i,data,dq,mask,seg,x,y in cutouts:
        for name,array in [('SCI',data),('DQ',dq),('MASK',mask.astype('u1')),('SEG',seg.astype('i4'))]:
            h=fits.ImageHDU(array,name=name);h.header['EXTVER']=i+1;h.header['SOURCEID']=str(audit['source_id'][i]);h.header['XCENT']=x;h.header['YCENT']=y;h.header['REASON']=str(audit['rejection_reason'][i]);stamps.append(h)
    fits.HDUList(stamps).writeto(paths.file('intermediate','alignment',label+'_rejected_cutouts','fits'),overwrite=True,checksum=True)
    for start in range(0,min(len(cutouts),options.get('max_diagnostic_cutouts',12)),4):
        page=cutouts[start:start+4];fig,axes=plt.subplots(len(page),4,figsize=(12,3*len(page)),squeeze=False)
        for row,(i,data,dq,mask,seg,x,y) in zip(axes,page):
            finite=data[np.isfinite(data)];lo,hi=ZScaleInterval().get_limits(finite) if finite.size else (0,1)
            for ax,array,title in zip(row,(data,dq!=0,mask,seg),('SCI','Input DQ','Effective mask','Segmentation')):
                ax.imshow(array,origin='lower',cmap='gray',vmin=lo if title=='SCI' else None,vmax=hi if title=='SCI' else None)
                ax.axvline(x,color='red',alpha=.6);ax.axhline(y,color='red',alpha=.6);ax.set_title(title+' '+str(i+1))
        fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment',label+f'_masked_cores_{start+1:03d}','png'),dpi=100);plt.close(fig)


def ensure_catalog(image,catalog,options,*,paths=None,label=None):
    """Gate imported/native/SEP catalogs before an STScI matcher sees them.

    A table produced by this QC version can be reused across header-only WCS
    candidates because its detector pixels/masks have not changed. Imported
    tables without that certificate must undergo pixel-level validation.
    """
    certificate=catalog.meta.get('centroid_quality',{})
    if certificate.get('version')==QUALITY_VERSION and 'centroid_valid' in catalog.colnames:
        keep=np.asarray(catalog['centroid_valid'],bool)
        if 'selected_for_astrometry' in catalog.colnames:keep&=np.asarray(catalog['selected_for_astrometry'],bool)
        return catalog[keep]
    return validate_centroids(image,catalog,options,method='custom',paths=paths,label=label)
