"""Collapse acquisition target images through a nominal long-slit footprint."""
from pathlib import Path
import json
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import astropy.units as u
from .acquisition import target_anchored_wcs
from .io import digest


def polygon_area(polygon):
    p=np.asarray(polygon,float)
    if len(p)<3:return 0.
    return float(abs(np.dot(p[:,0],np.roll(p[:,1],-1))-np.dot(p[:,1],np.roll(p[:,0],-1)))/2)


def clip_rectangle(polygon,xlo,xhi,ylo,yhi):
    p=[np.asarray(v,float) for v in polygon]
    for axis,bound,sign in [(0,xlo,1),(0,xhi,-1),(1,ylo,1),(1,yhi,-1)]:
        if not p:break
        result=[];previous=p[-1];prev_in=(previous[axis]-bound)*sign>=0
        for current in p:
            inside=(current[axis]-bound)*sign>=0
            if inside!=prev_in:
                fraction=(bound-previous[axis])/(current[axis]-previous[axis]);result.append(previous+fraction*(current-previous))
            if inside:result.append(current)
            previous=current;prev_in=inside
        p=result
    return np.asarray(p)


def slit_bin_weights(corners,edges,cross_center,width):
    """Exact planar pixel-area fractions; no rotated/smoothed image is made."""
    corners=np.asarray(corners,float);weights=np.zeros((len(edges)-1,len(corners)))
    lo=cross_center-width/2;hi=cross_center+width/2
    for j,p in enumerate(corners):
        if p[:,0].max()<=lo or p[:,0].min()>=hi:continue
        area=polygon_area(p)
        if area<=0:continue
        bins=np.flatnonzero((edges[:-1]<p[:,1].max())&(edges[1:]>p[:,1].min()))
        for b in bins:weights[b,j]=polygon_area(clip_rectangle(p,lo,hi,edges[b],edges[b+1]))/area
    return weights


def collapse_acquisitions(metadata,output,spatial_limit=2.,bin_width=None,background_annulus=(1.5,2.3)):
    meta=json.loads(Path(metadata).read_text());results=[];provenance=[];centers=None
    for entry in meta['acquisitions']:
        path=Path(entry['raw']);before=digest(path);panel=next(p for p in entry['panels'] if p['sci_extver']==2)
        slit=entry['science'][0];pa=np.deg2rad(slit['pa_aper']);width=float(slit['aperture'].upper().split('X')[1])
        with fits.open(path) as hd:
            header=hd['SCI',2].header;data=hd['SCI',2].data.astype(float);ph=hd[0].header
            origin=SkyCoord(ph['RA_TARG']*u.deg,ph['DEC_TARG']*u.deg);w=target_anchored_wcs(header,panel['target_pixel_zero_based'],origin)
            if centers is None:
                actual_width=float(bin_width or panel['pixel_scale_arcsec']);n=int(np.floor(spatial_limit/actual_width));centers=np.arange(-n,n+1)*actual_width;edges=np.r_[centers-actual_width/2,centers[-1]+actual_width/2]
            ny,nx=data.shape;yy,xx=np.indices(data.shape)
            xc=(xx.ravel()[:,None]+np.array([-.5,.5,.5,-.5]));yc=(yy.ravel()[:,None]+np.array([-.5,-.5,.5,.5]))
            sky=w.pixel_to_world(xc,yc);east,north=origin.spherical_offsets_to(sky);east=east.arcsec;north=north.arcsec
            cross=east*np.cos(pa)-north*np.sin(pa);along=east*np.sin(pa)+north*np.cos(pa);corners=np.stack([cross,along],axis=-1)
            sc=SkyCoord(slit['ra_aper']*u.deg,slit['dec_aper']*u.deg);ee,nn=origin.spherical_offsets_to(sc);cross_center=float(ee.arcsec*np.cos(pa)-nn.arcsec*np.sin(pa))
            values=data.ravel();valid=np.isfinite(values)
            dq=hd['DQ',2].data
            if dq is not None:valid&=(np.asarray(dq).ravel()&int(header.get('SDQFLAGS',31743)))==0
            radius=np.hypot(east.mean(axis=1),north.mean(axis=1));sky_mask=valid&(radius>=background_annulus[0])&(radius<=background_annulus[1])
            if sky_mask.sum()<100:raise ValueError('Insufficient source-free sky pixels')
            sky_mean,background,sky_sigma=sigma_clipped_stats(values[sky_mask],sigma=3,maxiters=5)
            weights=slit_bin_weights(corners,edges,cross_center,width);weights[:,~valid]=0
            flux=weights@(np.where(valid,values,0)-background)
            area=np.array([polygon_area(p) for p in corners]);coverage=(weights@area)/(width*np.diff(edges))
            if np.min(coverage)<.98:raise ValueError('Incomplete acquisition/slit profile coverage')
            peak=float(flux.max())
            if peak<=0:raise ValueError('Acquisition target has no positive collapsed peak')
            normalized=flux/peak;results.append(normalized)
            # Diagnostic sky-placement alternatives are not confidence limits.
            sky_alt=[]
            for lo,hi in [(1.25,1.75),(1.75,2.3)]:
                mask=valid&(radius>=lo)&(radius<=hi);bg=float(sigma_clipped_stats(values[mask],sigma=3,maxiters=5)[1]);p=weights@(np.where(valid,values,0)-bg);sky_alt.append(dict(annulus_arcsec=[lo,hi],background_counts=bg,peak_normalized=(p/p.max()).tolist()))
            provenance.append(dict(root=ph['ROOTNAME'].strip(),path=str(path),sha256=before,filter=ph['APERTURE'],extver=2,target_pixel_zero_based=panel['target_pixel_zero_based'],slit_width_arcsec=width,slit_pa_deg=slit['pa_aper']%360,slit_cross_center_arcsec=cross_center,background_annulus_arcsec=list(background_annulus),background_counts=float(background),sky_sigma_counts=float(sky_sigma),sky_pixels=int(sky_mask.sum()),raw_counts_collapsed=(weights@np.where(valid,values,0)).tolist(),background_subtracted_counts=flux.tolist(),peak_normalized=normalized.tolist(),coverage=coverage.tolist(),background_alternatives=sky_alt))
        assert digest(path)==before
    curves=np.asarray(results);mean=curves.mean(axis=0);den=float(mean.max());mean/=den;curves/=den
    fields=dict(offset_along_slit_arcsec=centers,profile_peak_normalized=mean,visit_min=curves.min(axis=0),visit_max=curves.max(axis=0))
    for item,p in zip(provenance,curves):fields['profile_'+item['root']]=p
    table=Table(fields);table.meta.update(filter=provenance[0]['filter'],slit_width_arcsec=provenance[0]['slit_width_arcsec'],bin_width_arcsec=actual_width,n_visits=len(results),coordinate='Along slit relative to onboard optical acquisition centroid',normalization='Each visit peak normalized, then equal-visit mean and peak normalization',uncertainty='Visit min/max is descriptive repeat spread, not a confidence interval; acquisition ERR arrays are unavailable')
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True);table.write(out,format='ascii.ecsv',overwrite=True)
    out.with_suffix('.json').write_text(json.dumps(dict(method='SCI,2 target images; exact pixel-polygon overlap in slit coordinates; robust off-source constant sky; no smoothing',metadata_source=str(metadata),visits=provenance,bin_width_arcsec=actual_width,limitations='Optical acquisition profile includes CCD/filter PSF and host light; not a FUV PSF. Registration to UV trace is display-only; slit placement remains nominal.'),indent=2)+'\n')
    return table
