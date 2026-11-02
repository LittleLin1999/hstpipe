"""Recover repeated HST sources with bounded, mask-aware native centroids."""
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
        if all(str(catalogs[i]['detection_method'][k]) in ('PSFPhotometry','aperture_profile') for i,k in map(node,rows)):guides.append((sid,rows))
    with fits.open(inputs[0]) as hd:tangent=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
    xy=[project(c,p,tangent) for c,p in zip(catalogs,inputs)];records=[]
    for c in catalogs:
        if 'centroid_ignore_bits' not in c.colnames:c['centroid_ignore_bits']=np.full(len(c),int(detection.get('centroid_ignore_bits',4096)),dtype='i8')
    for sid,rows in groups.items():
        if len(rows)<3:continue
        source_nodes=[node(r) for r in rows];present={i for i,k in source_nodes}
        ref,k=min(source_nodes,key=lambda n:float(catalogs[n[0]]['centroid_error_pix'][n[1]]))
        for i,image in enumerate(inputs):
            if i in present:continue
            guide_id=None;predicted=None
            for gid,grows in sorted(guides,key=lambda item:-len(item[1])):
                if gid==sid:continue
                guide={byname[str(r['exposure'])]:node(r)[1] for r in grows}
                if i not in guide or ref not in guide:continue
                predicted=xy[ref][k]+xy[i][guide[i]]-xy[ref][guide[ref]];guide_id=gid;break
            prior='measured_HST_guide' if predicted is not None else 'existing_WCS_of_repeated_HST_source'
            if predicted is None:
                positions=np.asarray([xy[e][n] for e,n in source_nodes]);predicted=np.median(positions,axis=0)
                if np.max(np.linalg.norm(positions-predicted,axis=1))>relative['max_current_rms_pix']:
                    records.append(dict(source_id=sid,exposure=image.name,status='search_prior_unconstrained',search_prior=prior,relative_validation_required=True));continue
            # The existing WCS sets only a bounded search window. It never
            # supplies a measured centroid or an accepted transformation.
            world=tangent.wcs_pix2world([predicted],0)[0]
            record=dict(source_id=sid,exposure=image.name,reference=inputs[ref].name,guide=guide_id,status='no_valid_measurement',external_coordinates_used_for_fit=False,search_prior=prior,relative_validation_required=True)
            with fits.open(image) as hd:
                for h in hd:
                    if h.name!='SCI':continue
                    v=int(h.header.get('EXTVER',1))
                    try:x,y=detector_wcs(hd,v).all_world2pix([world],0)[0]
                    except Exception as exc:record['reason']='WCS projection: '+str(exc);continue
                    if not np.isfinite([x,y]).all() or x<8 or y<8 or x>=h.shape[1]-8 or y>=h.shape[0]-8:continue
                    a=np.asarray(h.data,float);dq,base,mask=pixel_masks(hd,v,detection);err=np.asarray(image_data(hd['ERR',v]),float)
                    record['mask_attempts']=[];measured=None
                    policies=[('configured_mask',dict(detection))]
                    ix,iy=int(round(x)),int(round(y))
                    inherited=bool(np.any(dq[iy-8:iy+9,ix-8:ix+9]&4096)) and not hd[0].header.get('HPCRRUN',False)
                    if inherited and int(detection.get('centroid_ignore_bits',4096))&4096:
                        policies.append(('exclude_inherited_4096',dict(detection,centroid_ignore_bits=int(detection.get('centroid_ignore_bits',4096))&~4096)))
                    for policy,measurement_options in policies:
                        attempt=dict(policy=policy,centroid_ignore_bits=measurement_options.get('centroid_ignore_bits',4096),status='rejected')
                        try:
                            dq,base,mask=pixel_masks(hd,v,measurement_options)
                            measured=_local_measurement(image,v,a,err,dq,base,mask,x,y,measurement_options,relative)
                            attempt.update(status='measured_candidate',window_stability_pix=measured['stability'])
                        except Exception as exc:attempt['reason']=str(exc)
                        record['mask_attempts'].append(attempt)
                        if measured is not None:break
                    if measured is None:
                        record['reason']=record['mask_attempts'][-1]['reason'];continue
                    try:
                        center=measured['center'];stability=measured['stability'];clean=measured['clean'];flux=measured['flux']
                        template={n:catalogs[ref][n][k] for n in catalogs[ref].colnames}
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
                        template['centroid_ignore_bits']=int(measurement_options.get('centroid_ignore_bits',4096))
                        near=np.flatnonzero((np.asarray(catalogs[i]['extver'])==v)&(np.hypot(catalogs[i]['x']-center[0],catalogs[i]['y']-center[1])<relative.get('dedup_radius_pix',2.)))
                        if len(near)>1:raise ValueError('Multiple existing detections near recovered centroid; association remains ambiguous')
                        if len(near)==1:
                            old_id=str(catalogs[i]['hst_source_id'][near[0]])
                            if any(str(q['exposure'])==image.name and str(q['hst_source_id'])==old_id for rs in groups.values() for q in rs):raise ValueError('Recovered centroid conflicts with another common-source identity')
                            record['replaced_detector_measurement']=dict(hst_source_id=str(catalogs[i]['hst_source_id'][near[0]]),x=float(catalogs[i]['x'][near[0]]),y=float(catalogs[i]['y'][near[0]]))
                            catalogs[i][int(near[0])]=template
                        else:catalogs[i].add_row(template)
                        record.update(status='measured_candidate',extver=v,x=float(center[0]),y=float(center[1]),predicted_native=[float(x),float(y)],window_stability_pix=stability,centroid_ignore_bits=template['centroid_ignore_bits'],core_valid_fraction=float(clean['core_valid_fraction'][0]),centroid_error_pix=float(clean['centroid_error_pix'][0]),hst_source_id=identity)
                        break
                    except Exception as exc:record['reason']=str(exc)
            records.append(record)
    return records


def _local_measurement(image,v,a,err,dq,base,mask,x,y,detection,relative):
    """Existing Photutils peak/centroid measurement under one explicit mask."""
    fits_xy=[]
    ix=int(round(x));iy=int(round(y));sy=slice(iy-8,iy+9);sx=slice(ix-8,ix+9)
    data=a[sy,sx];bad=mask[sy,sx];_,sky,noise=sigma_clipped_stats(data,mask=bad)
    threshold=sky+detection.get('snr',5)*max(float(noise),float(np.median(err[sy,sx][~bad])))
    peaks=find_peaks(data,threshold,mask=bad,box_size=3)
    if peaks is None:raise ValueError('No independently detected peak in the HST-predicted region')
    peaks=peaks[np.hypot(peaks['x_peak']+sx.start-x,peaks['y_peak']+sy.start-y)<=relative['match_tolerance_pix']]
    if len(peaks)!=1:raise ValueError(f'{len(peaks)} plausible local peaks; no nearest/brightest choice')
    ix=int(peaks['x_peak'][0])+sx.start;iy=int(peaks['y_peak'][0])+sy.start
    peak=[ix,iy]
    for radius in (3,5):
        sy=slice(iy-radius,iy+radius+1);sx=slice(ix-radius,ix+radius+1)
        data=a[sy,sx];bad=mask[sy,sx];edge=np.zeros(data.shape,bool);edge[[0,-1],:]=True;edge[:,[0,-1]]=True
        if np.sum(edge&~bad)<10:raise ValueError('Insufficient uncontaminated local border')
        _,sky,noise=sigma_clipped_stats(data[edge&~bad])
        with warnings.catch_warnings():
            warnings.simplefilter('error');cx,cy=centroid_2dg(data-sky,error=err[sy,sx],mask=bad)
        fits_xy.append(np.array([cx+sx.start,cy+sy.start]))
    if not np.isfinite(fits_xy).all():raise ValueError('Nonfinite local centroid')
    center=fits_xy[0];stability=float(np.linalg.norm(fits_xy[0]-fits_xy[1]))
    if not np.isfinite(center).all() or np.linalg.norm(center-[x,y])>relative['match_tolerance_pix']:raise ValueError('Local centroid did not remain within the HST-guided association region')
    floor=detection.get('galaxy_error_floor',.2)
    if stability>max(detection.get('stability_floor_pix',.5),3*floor):raise ValueError('Centroid depends on fitting window')
    ctx=dict(data=a,dq=dq,base=base,mask=mask,error=err,sky=float(sky),rms=float(noise))
    raw=Table(dict(x=[center[0]],y=[center[1]],extver=[v],local_sky=[float(sky)],centroid_error_pix=[max(floor,stability)]))
    clean,audit=validate_centroids(image,raw,detection,method='HST_guided_centroid_2dg',contexts={v:ctx},return_all=True)
    if not len(clean):raise ValueError('Pixel/centroid QC: '+str(audit['quality_flags'][0]))
    yy,xx=np.indices(data.shape);aperture=(xx+sx.start-center[0])**2+(yy+sy.start-center[1])**2<=3.5**2
    flux=float(np.sum((data-sky)[aperture&~bad]))
    return dict(center=center,stability=stability,clean=clean,flux=flux,peak=peak)


def recovery_figure(inputs,records,paths):
    """Explain a mask-dependent recovery without changing pixels or WCS."""
    accepted=[r for r in records if r.get('status')=='measured_candidate' and len(r.get('mask_attempts',[]))>1]
    if not accepted:return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import ImageNormalize,AsinhStretch
    with plt.rc_context({'font.family':'STIXGeneral','mathtext.fontset':'stix','xtick.direction':'in','ytick.direction':'in','xtick.minor.visible':True,'ytick.minor.visible':True,'xtick.top':True,'ytick.right':True,'ytick.left':True,'font.size':12}):
        fig,axes=plt.subplots(len(accepted),2,figsize=(9,4.4*len(accepted)),squeeze=False,layout='constrained')
        for row,r in enumerate(accepted):
            p=next(p for p in inputs if p.name==r['exposure']);v=r['extver'];x,y=r['x'],r['y'];ix,iy=round(x),round(y);radius=10
            with fits.open(p) as hd:
                data=hd['SCI',v].data[iy-radius:iy+radius+1,ix-radius:ix+radius+1].astype(float)
                dq=np.asarray(image_data(hd['DQ',v]),dtype='u2')[iy-radius:iy+radius+1,ix-radius:ix+radius+1]
            _,sky,noise=sigma_clipped_stats(data);norm=ImageNormalize(vmin=sky-noise,vmax=sky+25*noise,stretch=AsinhStretch(.1))
            extent=(ix-radius-.5,ix+radius+.5,iy-radius-.5,iy+radius+.5)
            for col in (0,1):
                ax=axes[row,col];shown=np.ma.masked_where((dq&4096)!=0,data) if col else data
                ax.set_facecolor('#e6a949');ax.imshow(shown,origin='lower',extent=extent,cmap='gray',norm=norm,interpolation='nearest')
                ax.plot(*r['predicted_native'],'+',color='#00cfff',ms=14,mew=1.5)
                if col:ax.plot(x,y,'o',mfc='none',mec='#20d382',ms=16,mew=1.5)
                ax.set_title(r['exposure'].removesuffix('_flc.fits')+'\n'+('Inherited 4096 included: measurement failed' if col==0 else 'Inherited 4096 excluded: centroid recovered'))
                ax.set_xlabel('Native detector X (pixels)');ax.set_ylabel('Native detector Y (pixels)')
            axes[row,1].text(.03,.03,f"Core valid: {r['core_valid_fraction']:.0%}\nWindow change: {r['window_stability_pix']:.4f} pixel",transform=axes[row,1].transAxes,color='white',fontsize=10,bbox=dict(facecolor='black',alpha=.7,pad=3))
        fig.suptitle('Conditional centroid recovery — same original SCI pixels\nGold: excluded 4096 pixels; cyan: search prior; green: independently measured centroid',fontsize=12)
        output=paths.file('diagnostics','alignment','centroid_recovery','jpg');output.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(output,dpi=135,pil_kwargs=dict(quality=86,optimize=True));plt.close(fig)
    return output
