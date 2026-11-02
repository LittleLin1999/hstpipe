"""HST source detection and spherical match/validation statistics; no WCS fitting."""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.wcs import WCS
from .wcs import detector_wcs
from ..fits import image_data
from astropy.stats import sigma_clipped_stats,sigma_clip


def detect_sources(image,options,*,paths=None,label=None,return_all=False):
    """DAO detections go through the same pixel-level QC as other detectors."""
    from photutils.detection import DAOStarFinder
    from .centroids import pixel_masks,validate_centroids
    tables=[];contexts={}
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
            scale=np.sqrt(abs(np.linalg.det(detector_wcs(hd,h.header.get('EXTVER',1)).pixel_scale_matrix)))*3600
            fwhm=options['fwhm_arcsec']/scale if options.get('fwhm_arcsec') is not None else options.get('fwhm_pixels',2.5)
            cat=DAOStarFinder(options.get('snr',5)*rms,fwhm=fwhm)(np.where(mask,sky,a)-sky,mask=mask)
            if cat is None:continue
            xkey='x_centroid' if 'x_centroid' in cat.colnames else 'xcentroid';ykey='y_centroid' if 'y_centroid' in cat.colnames else 'ycentroid'
            x=np.asarray(cat[xkey]);y=np.asarray(cat[ykey]);world=detector_wcs(hd,h.header.get('EXTVER',1)).all_pix2world(np.column_stack([x,y]),0)
            tables.append(Table(dict(x=x,y=y,ra=world[:,0],dec=world[:,1],flux=np.asarray(cat['flux']),extver=np.full(len(x),v))))
        raw=vstack(tables) if tables else Table(names=['x','y','ra','dec','flux','extver'],dtype=['f8']*5+['i8'])
        result=validate_centroids(image,raw,options,method='DAOStarFinder',paths=paths,label=label,return_all=return_all,contexts=contexts)
        targets=result if not return_all else result[0]
        targets.meta['detection_chips']=[dict(extver=int(v),sky=c['sky'],noise=c['rms'],masked_fraction=float(c['mask'].mean()),threshold_snr=options.get('snr',5),error_floor=True) for v,c in contexts.items()]
        return result



def matches(sources,reference,options):
    names=['source_index','reference_index','x','y','ra','dec','ref_ra','ref_dec','dx_arcsec','dy_arcsec','accepted']
    if not len(sources) or not len(reference):return Table(names=names,dtype=['i8','i8']+['f8']*8+['bool'])
    sky=SkyCoord(np.asarray(sources['ra'])*u.deg,np.asarray(sources['dec'])*u.deg)
    ref=SkyCoord(np.asarray(reference['ra'])*u.deg,np.asarray(reference['dec'])*u.deg)
    index,distance,_=sky.match_to_catalog_sky(ref);reverse,_,_=ref.match_to_catalog_sky(sky)
    good=(distance.to_value(u.arcsec)<=options['search_radius_arcsec'])&(reverse[index]==np.arange(len(sky)))
    if len(ref)>1:
        _,second,_=sky.match_to_catalog_sky(ref,nthneighbor=2)
        good&=second.to_value(u.arcsec)>np.maximum(options.get('ambiguity_arcsec',.2),distance.to_value(u.arcsec)*options.get('ambiguity_ratio',2))
    selected=np.flatnonzero(good);ri=index[good]
    dx,dy=sky[selected].spherical_offsets_to(ref[ri]);dx=dx.to_value(u.arcsec);dy=dy.to_value(u.arcsec)
    if len(selected):
        # Prevent clipping spherical round-off when exact fixture/identity matches have zero scatter.
        def scatter(a,axis=None):return np.maximum(np.nanstd(a,axis=axis),options.get('numeric_floor_arcsec',1e-8))
        clipped_x=sigma_clip(dx,sigma=options.get('clip_sigma',3),maxiters=5,stdfunc=scatter);clipped_y=sigma_clip(dy,sigma=options.get('clip_sigma',3),maxiters=5,stdfunc=scatter)
        keep=~(np.ma.getmaskarray(clipped_x)|np.ma.getmaskarray(clipped_y))
    else:keep=np.array([],bool)
    return Table([selected,ri,np.asarray(sources['x'])[selected],np.asarray(sources['y'])[selected],
        np.asarray(sources['ra'])[selected],np.asarray(sources['dec'])[selected],np.asarray(reference['ra'])[ri],np.asarray(reference['dec'])[ri],dx,dy,keep],names=names)


def validate(sources,reference,options,shape=None):
    table=matches(sources,reference,options);kept=table[table['accepted']]
    result=dict(status='insufficient_matches',matched=len(kept),rejected=len(table)-len(kept),detected=len(sources),reference_count=len(reference),
                thresholds=dict(options),source_quality=sources.meta.get('centroid_quality'),rms_arcsec=None,robust_scatter_arcsec=None,offset_arcsec=None,uncertainty_arcsec=None)
    if len(kept)<options['min_matches']:return result,table
    dx=np.asarray(kept['dx_arcsec']);dy=np.asarray(kept['dy_arcsec']);offset=np.array([np.median(dx),np.median(dy)])
    scatter=float(1.4826*np.sqrt(np.median(abs(dx-offset[0]))**2+np.median(abs(dy-offset[1]))**2))
    rms=float(np.sqrt(np.mean(dx**2+dy**2)))
    span=min(np.ptp(kept['x'])/shape[1],np.ptp(kept['y'])/shape[0]) if shape is not None else None
    result.update(rms_arcsec=rms,robust_scatter_arcsec=scatter,offset_arcsec=offset.tolist(),span_fraction=span,
        uncertainty_arcsec=float(np.hypot(scatter/np.sqrt(len(kept)),options['systematic_floor_arcsec'])),
        uncertainty_caveat='Fit-star scatter is not independent external accuracy; catalog systematics and centroid/color effects remain')
    passed=(rms<=options['max_rms_arcsec'] and np.linalg.norm(offset)<=options['max_offset_arcsec'] and
            (span is None or span>=options.get('min_span_fraction',.1)))
    result['status']='validated' if passed else 'unreliable_residuals_or_distribution'
    return result,table


def science_header(path):
    with fits.open(path) as hd:return hd['SCI'].header.copy()


def science_shape(path):
    with fits.open(path) as hd:return hd['SCI'].shape
