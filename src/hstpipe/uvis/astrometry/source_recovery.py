"""Recover missed measurements from a repeated HST object and a local HST guide."""
import warnings
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from photutils.centroids import centroid_2dg
from photutils.detection import find_peaks
from .wcs import detector_wcs
from .sources import project
from .centroids import pixel_masks,validate_centroids
from ..fits import image_data

def recover_missing(inputs,catalogs,common,detection,relative):
    if not len(common):return []
    groups={}
    for r in common:
        if r['status']=='common_candidate':groups.setdefault(str(r['source_id']),[]).append(r)
    byname={p.name:i for i,p in enumerate(inputs)}
    def node(row):
        i=byname[str(row['exposure'])];c=catalogs[i]
        k=int(np.flatnonzero(np.asarray(c['hst_source_id']).astype(str)==str(row['hst_source_id']))[0])
        return i,k
    guides=[]
    for sid,rows in groups.items():
        if len(rows)<2:continue
        if all(str(catalogs[i]['detection_method'][k])=='PSFPhotometry' for i,k in map(node,rows)):guides.append((sid,rows))
    if not guides:return []
    guide_id,guide_rows=max(guides,key=lambda item:len(item[1]));guide={byname[str(r['exposure'])]:node(r)[1] for r in guide_rows}
    with fits.open(inputs[0]) as hd:tangent=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
    xy=[project(c,p,tangent) for c,p in zip(catalogs,inputs)];records=[]
    for sid,rows in groups.items():
        if sid==guide_id or len(rows)<3:continue
        refs=[node(r) for r in rows if byname[str(r['exposure'])] in guide]
        if not refs:continue
        ref,k=min(refs,key=lambda n:(str(catalogs[n[0]]['detection_method'][n[1]])!='PSFPhotometry',float(catalogs[n[0]]['centroid_error_pix'][n[1]])))
        present={byname[str(r['exposure'])] for r in rows}
        for i,image in enumerate(inputs):
            if i in present or i not in guide:continue
            predicted=xy[ref][k]+xy[i][guide[i]]-xy[ref][guide[ref]]
            world=tangent.wcs_pix2world([predicted],0)[0]
            record=dict(source_id=sid,exposure=image.name,reference=inputs[ref].name,guide=guide_id,status='no_valid_measurement',external_coordinates_used_for_fit=False)
            with fits.open(image) as hd:
                for h in hd:
                    if h.name!='SCI':continue
                    v=int(h.header.get('EXTVER',1))
                    try:x,y=detector_wcs(hd,v).all_world2pix([world],0)[0]
                    except Exception as exc:record['reason']='WCS projection: '+str(exc);continue
                    if not np.isfinite([x,y]).all() or x<8 or y<8 or x>=h.shape[1]-8 or y>=h.shape[0]-8:continue
                    a=np.asarray(h.data,float);dq,base,mask=pixel_masks(hd,v,detection);err=np.asarray(image_data(hd['ERR',v]),float)
                    fits_xy=[];sky=0.;contexts=None
                    try:
                        ix=int(round(x));iy=int(round(y));sy=slice(iy-8,iy+9);sx=slice(ix-8,ix+9)
                        data=a[sy,sx];bad=mask[sy,sx];_,sky,noise=sigma_clipped_stats(data,mask=bad)
                        threshold=sky+detection.get('snr',5)*max(float(noise),float(np.median(err[sy,sx][~bad])))
                        peaks=find_peaks(data,threshold,mask=bad,box_size=3)
                        if peaks is None:raise ValueError('No independently detected peak in the HST-predicted region')
                        peaks=peaks[np.hypot(peaks['x_peak']+sx.start-x,peaks['y_peak']+sy.start-y)<=relative['match_tolerance_pix']]
                        if len(peaks)!=1:raise ValueError(f'{len(peaks)} plausible local peaks; no nearest/brightest choice')
                        ix=int(peaks['x_peak'][0])+sx.start;iy=int(peaks['y_peak'][0])+sy.start
                        record.update(detected_peak_native=[ix,iy],status='detected_but_rejected')
                        for radius in (3,5):
                            sy=slice(iy-radius,iy+radius+1);sx=slice(ix-radius,ix+radius+1)
                            data=a[sy,sx];bad=mask[sy,sx];edge=np.zeros(data.shape,bool);edge[[0,-1],:]=True;edge[:,[0,-1]]=True
                            if np.sum(edge&~bad)<10:raise ValueError('Insufficient uncontaminated local border')
                            _,sky,noise=sigma_clipped_stats(data[edge&~bad])
                            with warnings.catch_warnings():
                                warnings.simplefilter('error');cx,cy=centroid_2dg(data-sky,error=err[sy,sx],mask=bad)
                            fits_xy.append(np.array([cx+sx.start,cy+sy.start]))
                        center=fits_xy[0];stability=float(np.linalg.norm(fits_xy[0]-fits_xy[1]))
                        if not np.isfinite(center).all() or np.linalg.norm(center-[x,y])>relative['match_tolerance_pix']:raise ValueError('Local centroid did not remain within the HST-guided association region')
                        floor=detection.get('galaxy_error_floor',.2)
                        if stability>max(detection.get('stability_floor_pix',.5),3*floor):raise ValueError('Centroid depends on fitting window')
                        ctx=dict(data=a,dq=dq,base=base,mask=mask,error=err,sky=float(sky),rms=float(noise))
                        raw=Table(dict(x=[center[0]],y=[center[1]],extver=[v],local_sky=[float(sky)],centroid_error_pix=[max(floor,stability)]))
                        clean,audit=validate_centroids(image,raw,detection,method='HST_guided_centroid_2dg',contexts={v:ctx},return_all=True)
                        if not len(clean):raise ValueError('Pixel/centroid QC: '+str(audit['quality_flags'][0]))
                        yy,xx=np.indices(data.shape);aperture=(xx+sx.start-center[0])**2+(yy+sy.start-center[1])**2<=3.5**2
                        flux=float(np.sum((data-sky)[aperture&~bad]));template={n:catalogs[ref][n][k] for n in catalogs[ref].colnames}
                        for n in ('source_id','hst_source_id','detector_source_id','exposure_rootname'):
                            if n in catalogs[i].colnames:catalogs[i][n]=np.asarray(catalogs[i][n]).astype('U128')
                        for n in ('detection_method','quality_flags','rejection_reason'):
                            if n in catalogs[i].colnames:catalogs[i][n]=np.asarray(catalogs[i][n]).astype('U256')
                        for n in clean.colnames:
                            if n in template:template[n]=clean[n][0]
                        root=hd[0].header.get('ROOTNAME',image.stem);identity=root+'_recovered_'+sid
                        template.update(x=float(center[0]),y=float(center[1]),extver=v,flux=flux,hst_source_id=identity,source_id=identity,exposure_rootname=root,detection_method='HST_guided_centroid_2dg')
                        if 'detector_source_id' in template:template['detector_source_id']=identity
                        if 'ra' in template:
                            ra,dec=detector_wcs(hd,v).all_pix2world([center],0)[0];template.update(ra=ra,dec=dec)
                        if 'ccdchip' in template:template['ccdchip']=h.header.get('CCDCHIP',-1)
                        catalogs[i].add_row(template)
                        record.update(status='measured_candidate',extver=v,x=float(center[0]),y=float(center[1]),predicted_native=[float(x),float(y)],window_stability_pix=stability,centroid_error_pix=float(clean['centroid_error_pix'][0]),hst_source_id=identity)
                        break
                    except Exception as exc:record['reason']=str(exc)
            records.append(record)
    return records
