"""Source measurements for relative registration, independent of WCS fitting."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from .wcs import detector_wcs
from ..fits import image_data
from scipy.ndimage import gaussian_filter,binary_dilation

def sources(image,options,extended=False,width=3.5,*,paths=None,label=None,return_all=False):
    """Return clean catalogs; every detector uses shared mandatory centroid QC."""
    from .centroids import pixel_masks,validate_centroids
    if not extended:
        from .matching import detect_sources
        return detect_sources(image,options,paths=paths,label=label,return_all=return_all)
    from photutils.segmentation import detect_sources,deblend_sources,SourceCatalog
    import inspect
    area_key='n_pixels' if 'n_pixels' in inspect.signature(detect_sources).parameters else 'npixels'
    level_key='n_levels' if 'n_levels' in inspect.signature(deblend_sources).parameters else 'nlevels'
    tables=[];segments={};contexts={}
    with fits.open(image) as hd:
        for h in hd:
            if h.name!='SCI':continue
            v=h.header.get('EXTVER',1);a=np.asarray(h.data,float);dq,base,centroid_mask=pixel_masks(hd,v,options,'centroid')
            _,_,mask=pixel_masks(hd,v,options,'detection')
            _,sky,rms=sigma_clipped_stats(a,mask=mask)
            if ('ERR',v) in hd:
                usable=np.asarray(image_data(hd['ERR',v]))[~mask];usable=usable[np.isfinite(usable)&(usable>0)]
                if usable.size:rms=max(float(rms),float(np.median(usable)))
            if not np.isfinite(rms) or rms<=0:continue
            error=np.asarray(image_data(hd['ERR',v]),float) if ('ERR',v) in hd else np.full(a.shape,rms)
            contexts[v]=dict(data=a,dq=dq,base=base,mask=centroid_mask,sky=float(sky),rms=float(rms),error=error)
            # Fill is for convolution/detection only; morphology and QC use original SCI.
            data=np.where(mask,0,a-sky);sigma=width/2.355;smooth=gaussian_filter(data,sigma)
            threshold=options.get('snr',5)*rms/max(1,2*np.sqrt(np.pi)*sigma)
            if 'WHT' in hd:
                # Diagnostic mosaics have correlated noise, unrejected CRs and
                # possible residual sky. The independent-pixel convolution
                # formula can otherwise merge entire CCD footprints into a
                # giant segment. Measure the actual filtered distribution;
                # this changes detection only, never SCI or centroid pixels.
                _,filtered_sky,filtered_rms=sigma_clipped_stats(smooth,mask=mask)
                if not np.isfinite(filtered_sky+filtered_rms) or filtered_rms<=0:continue
                threshold=max(threshold,float(filtered_sky+options.get('snr',5)*filtered_rms))
                contexts[v].update(filtered_background=float(filtered_sky),filtered_noise=float(filtered_rms),
                    detection_threshold=float(threshold),detection_noise_method='Measured filtered mosaic statistics; correlated noise and CR contamination included')
            seg=detect_sources(smooth,threshold,mask=mask,**{area_key:options.get('min_area',9)})
            if seg is None:continue
            seg=deblend_sources(smooth,seg,contrast=.005,progress_bar=False,**{area_key:options.get('min_area',9),level_key:32})
            segments[v]=seg.data
            cat=SourceCatalog(data,seg,error=error,mask=mask)
            x=np.asarray(cat.x_centroid if hasattr(cat,'x_centroid') else cat.xcentroid);y=np.asarray(cat.y_centroid if hasattr(cat,'y_centroid') else cat.ycentroid)
            major=(cat.semimajor_axis if hasattr(cat,'semimajor_axis') else cat.semimajor_sigma).value;minor=(cat.semiminor_axis if hasattr(cat,'semiminor_axis') else cat.semiminor_sigma).value
            flux=np.asarray(cat.segment_flux);snr=flux/np.maximum(np.asarray(cat.segment_flux_err if hasattr(cat,'segment_flux_err') else cat.segment_fluxerr),1e-12)
            world=detector_wcs(hd,h.header.get('EXTVER',1)).all_pix2world(np.column_stack([x,y]),0)
            tables.append(Table(dict(x=x,y=y,ra=world[:,0],dec=world[:,1],flux=flux,extver=np.full(len(x),v),width_pix=major,minor_width_pix=minor,snr=snr,segment_label=cat.labels,kind=np.where(major>2,'compact','point'))))
        raw=vstack(tables) if tables else Table(names=['x','y','ra','dec','flux','extver','width_pix','minor_width_pix','snr','segment_label','kind'],dtype=['f8']*5+['i8','f8','f8','f8','i8','U12'])
        result=validate_centroids(image,raw,options,method='segmentation',segmentations=segments,paths=paths,label=label,return_all=return_all,contexts=contexts)
        targets=result if not return_all else result[0]
        targets.meta['detection_chips']=[dict(extver=int(v),sky=c['sky'],noise=c['rms'],masked_fraction=float(c['mask'].mean()),threshold_snr=options.get('snr',5),error_floor=True,
            **{k:c[k] for k in ('filtered_background','filtered_noise','detection_threshold','detection_noise_method') if k in c}) for v,c in contexts.items()]
        return result


def project(catalog,image,tangent):
    """Apply full detector distortion, then the common linear TAN projection."""
    xy=np.empty((len(catalog),2))
    with fits.open(image) as hd:
        for v in np.unique(catalog['extver']):
            use=catalog['extver']==v;h=hd['SCI',int(v)]
            world=detector_wcs(hd,h.header.get('EXTVER',1)).all_pix2world(np.column_stack([catalog['x'][use],catalog['y'][use]]),0)
            xy[use]=tangent.wcs_world2pix(world,0)
    return xy


def plot_sources(image,catalog,destination):
    """Detector-source overlays for inspecting artifacts, chip edges and blends."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval
    with fits.open(image) as hd:
        chips=[h for h in hd if h.name=='SCI']
        fig,axes=plt.subplots(1,len(chips),figsize=(6*len(chips),5),squeeze=False)
        for h,ax in zip(chips,axes[0]):
            a=h.data;sample=a[::max(1,a.shape[0]//512),::max(1,a.shape[1]//512)]
            lo,hi=ZScaleInterval().get_limits(sample[np.isfinite(sample)])
            ax.imshow(sample,origin='lower',extent=(0,a.shape[1],0,a.shape[0]),vmin=lo,vmax=hi,cmap='gray')
            rows=catalog[catalog['extver']==h.header.get('EXTVER',1)]
            ax.scatter(rows['x'],rows['y'],s=40,facecolors='none',edgecolors='red',linewidths=.6)
            ax.set_title(f"{Path(image).name} SCI {h.header.get('EXTVER',1)}: {len(rows)} sources")
        fig.tight_layout();fig.savefig(destination,dpi=110);plt.close(fig)
