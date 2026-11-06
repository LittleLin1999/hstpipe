"""Top-hat pixel overlap alignment with full redistribution covariance.

X1D supplies a wavelength/EXTRLOCY coordinate map, never the spatial PSF.
Input FLT diagonal variances are propagated through the overlap operator.
"""
from pathlib import Path
import logging,json
import numpy as np
from scipy.linalg import pinvh
from scipy.optimize import curve_fit
from scipy.ndimage import median_filter,maximum_filter1d
from astropy.io import fits
from .common import *


def overlap_profile(image,variance,valid,trace,columns,radius=190,aligned=True):
    n=2*radius+1; flux=np.zeros(n);diag=np.zeros(n);off=np.zeros(n-1);coverage=np.zeros(n)
    for col in columns:
        center=float(trace[col] if aligned else np.median(trace[columns]))
        low=max(0,int(np.floor(center))-radius-2);high=min(image.shape[0],int(np.ceil(center))+radius+3)
        yy=np.arange(low,high);t=yy-center+radius;k=np.floor(t).astype(int);f=t-k
        good=valid[yy,col];v=variance[yy,col];s=image[yy,col]
        for delta,w in [(0,1-f),(1,f)]:
            j=k+delta;m=good&(j>=0)&(j<n)
            np.add.at(flux,j[m],s[m]*w[m]);np.add.at(diag,j[m],v[m]*w[m]**2);np.add.at(coverage,j[m],w[m])
        m=good&(k>=0)&(k<n-1)
        np.add.at(off,k[m],v[m]*f[m]*(1-f[m]))
    cov=np.diag(diag)+np.diag(off,1)+np.diag(off,-1)
    return np.arange(-radius,radius+1,dtype=float),flux,cov,coverage/max(len(columns),1)


def subtract_background(y,p,cov,coverage,window,order=0):
    bg=(abs(y)>=window[0])&(abs(y)<=window[1])&(coverage>0.95)&(np.diag(cov)>0)
    if bg.sum()<20:raise ValueError('Fewer than 20 valid background rows')
    x=np.vander(y/100,order+1,increasing=True)
    # Iterated robust rejection only in the off-source background windows.
    for _ in range(4):
        ii=np.where(bg)[0];cb=cov[np.ix_(ii,ii)];w=pinvh(cb,rtol=1e-8)
        inv=pinvh(x[ii].T@w@x[ii]);a=inv@x[ii].T@w@p[ii]
        resid=p-x@a; sigma=np.sqrt(np.maximum(np.diag(cov),1e-30))
        nextbg=bg&(abs(resid)<4*sigma)
        if nextbg.sum()<20 or np.array_equal(nextbg,bg):break
        bg=nextbg
    # Linear estimator H includes covariance with the background samples.
    h=np.zeros((order+1,len(y)));h[:,ii]=inv@x[ii].T@w
    operator=np.eye(len(y))-x@h
    outcov=operator@cov@operator.T
    return p-x@a,outcov,x@a,bg,a,inv


def normalize(p,cov):
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(cov)):raise ValueError('Nonfinite profile or uncertainty')
    total=p.sum()
    if not np.isfinite(total) or total<=0:raise ValueError('Nonpositive finite-window flux')
    jac=np.eye(len(p))/total-p[:,None]/total**2
    return p/total,jac@cov@jac.T


def coordinate_map(hx,nx):
    spec=hx[1].data[0]
    wave=np.array(spec['WAVELENGTH'],float);trace=np.array(spec['EXTRLOCY'],float)-1
    if len(wave)!=nx or len(trace)!=nx:raise ValueError('X1D/FLT sampling mismatch; reference-table mapping required')
    if not np.all(np.diff(wave)>0) or not np.all(np.isfinite(trace)):raise ValueError('Invalid wavelength/trace mapping')
    return wave,trace


def refine_trace(s,v,valid,reference,wave,c):
    xs=[];dy=[];errors=[]
    for lo in range(0,s.shape[1],c['trace_chunk']):
        cols=np.arange(lo,min(lo+c['trace_chunk'],s.shape[1]));cols=cols[(wave[cols]>1165)&(wave[cols]<1680)]
        for a,b in c.get('trace_exclude_wavelengths',c['exclude_wavelengths']):cols=cols[~((wave[cols]>=a)&(wave[cols]<=b))]
        if len(cols)<12:continue
        y,p,cv,coverage=overlap_profile(s,v,valid,reference,cols,25)
        bg=np.median(p[(abs(y)>18)&(coverage>.9)])
        m=abs(y)<8;q=p[m]-bg;z=y[m];den=q.sum();vv=cv[np.ix_(m,m)]
        sn=den/np.sqrt(np.sum(vv))
        if sn<15 or den<=0:continue
        center=np.sum(z*q)/den;grad=(z-center)/den;err=np.sqrt(max(grad@vv@grad,1e-8))
        if abs(center)<5:xs.append(np.mean(cols));dy.append(center);errors.append(err)
    if len(xs)<4:return reference,dict(status='CALIBRATED_TRACE_ONLY',rms=None,nblocks=len(xs))
    xs=np.array(xs);dy=np.array(dy);errors=np.array(errors)
    coeff=np.polyfit(xs,dy,2,w=1/errors);res=dy-np.polyval(coeff,xs)
    good=abs(res)<np.maximum(.3,4*errors)
    if good.sum()>=4:coeff=np.polyfit(xs[good],dy[good],2,w=1/errors[good])
    rms=float(np.sqrt(np.mean((dy[good]-np.polyval(coeff,xs[good]))**2)))
    return reference+np.polyval(coeff,np.arange(len(reference))),dict(status='MEASURED_TRACE',rms=rms,nblocks=len(xs),x=xs.tolist(),offset=dy.tolist(),errors=errors.tolist(),coeff=coeff.tolist())


def fractional_aperture(y,height):
    width=float(np.median(np.diff(y)))
    return np.clip((np.minimum(y+width/2,height/2)-np.maximum(y-width/2,-height/2))/width,0,1)


def measures(y,p,cov,heights,wing_inner=10.):
    core=abs(y)<8;z=y[core];q=p[core]
    center=float(np.sum(z*q)/np.sum(q));g=np.zeros(len(y));g[core]=(z-center)/q.sum()
    centroid_error=np.sqrt(max(g@cov@g,0))
    peak=np.argmax(p);half=p[peak]/2
    left=np.where(p[:peak]<half)[0];right=np.where(p[peak+1:]<half)[0]
    if len(left) and len(right):
        l=left[-1];r=peak+1+right[0]
        xl=np.interp(half,[p[l],p[l+1]],[y[l],y[l+1]])
        xr=np.interp(half,[p[r],p[r-1]],[y[r],y[r-1]])
        fwhm=float(xr-xl)
    else:fwhm=np.nan
    total=p.sum();norm,ncov=normalize(p,cov)
    result=dict(centroid=center,centroid_error=float(centroid_error),fwhm=fwhm,flux=float(total),snr=float(total/np.sqrt(np.sum(cov))))
    for h in heights:
        a=fractional_aperture(y,h)
        result[f'ee_{h}']=float(a@norm);result[f'ee_{h}_error']=float(np.sqrt(max(a@ncov@a,0)))
    a=1-fractional_aperture(y,2*wing_inner);result['wing_fraction']=float(a@norm)
    result['wing_error']=float(np.sqrt(max(a@ncov@a,0)))
    result['wing_snr']=float((a@p)/np.sqrt(max(a@cov@a,np.finfo(float).tiny)))
    a=np.sign(y);result['asymmetry']=float(a@norm);result['asymmetry_error']=float(np.sqrt(max(a@ncov@a,0)))
    try:
        def gauss(x,a,mu,sig,b):return a*np.exp(-.5*((x-mu)/sig)**2)+b
        pars,_=curve_fit(gauss,y[core],p[core]/total,p0=[max(p)/total,center,1.5,0],sigma=np.sqrt(np.maximum(np.diag(cov)[core],max(np.max(np.diag(cov))*1e-12,np.finfo(float).tiny)))/total,absolute_sigma=True,maxfev=3000)
        result['gaussian_sigma']=float(abs(pars[2]))
    except Exception:result['gaussian_sigma']=None
    if np.isfinite(fwhm):
        def width(q):
            pk=np.argmax(q);hh=q[pk]/2;ll=np.where(q[:pk]<hh)[0];rr=np.where(q[pk+1:]<hh)[0]
            if not len(ll) or not len(rr):return np.nan
            l=ll[-1];r=pk+1+rr[0]
            return np.interp(hh,[q[r],q[r-1]],[y[r],y[r-1]])-np.interp(hh,[q[l],q[l+1]],[y[l],y[l+1]])
        grad=np.zeros(len(p))
        for i in range(len(p)):
            eps=max(abs(p[i])*1e-5,abs(total)*1e-8);up=p.copy();dn=p.copy();up[i]+=eps;dn[i]-=eps;grad[i]=(width(up)-width(dn))/(2*eps)
        result['fwhm_error']=float(np.sqrt(max(grad@cov@grad,0)))
    else:result['fwhm_error']=None
    return result,norm,ncov


def extract_one(r,c,mask_spectral_emission=True):
    paths=r['paths'];meta=[]
    with fits.open(paths['flt']) as hf, fits.open(paths['x1d']) as hx:
        hdr=hf[0].header
        if hdr.get('DETECTOR')!='FUV-MAMA' or hdr.get('OPT_ELEM')!='G140L':raise ValueError('Wrong FITS configuration')
        for k in ['DETECTOR','APERTURE','OPT_ELEM','CENWAVE','CAL_VER','DARKCORR','FLATCORR','WAVECORR','PCTAB','SPTRCTAB','DISPTAB','DARKFILE','PFLTFILE','POSTARG1','POSTARG2','OBSMODE']:
            r[k.lower()]=hdr.get(k)
        for sci in [h for h in hf if h.name=='SCI']:
            ver=sci.header.get('EXTVER',1);s=np.array(sci.data,float);err=np.array(hf['ERR',ver].data,float);dq=np.array(hf['DQ',ver].data,np.uint16)
            if s.shape!=err.shape or s.shape!=dq.shape:raise ValueError('SCI/ERR/DQ mismatch')
            unit=sci.header.get('BUNIT','').upper();r['flt_unit']=unit
            if unit not in ['COUNTS','COUNT']:raise ValueError('Expected native FLT counts, found '+unit)
            v=np.maximum(err**2,1.0) # Conservative floor for low-count MAMA background pixels.
            flags=int(sci.header.get('SDQFLAGS',hdr.get('SDQFLAGS',31743)))
            valid=((dq&flags)==0)&np.isfinite(s)&np.isfinite(err)&(err>=0)
            wave,ref=coordinate_map(hx,s.shape[1]);net=np.array(hx[1].data[0]['NET'],float);ne=np.array(hx[1].data[0]['NET_ERROR'],float);smooth=median_filter(net,size=31);emission=(net>smooth+5*np.maximum(ne,1e-12))&(net>2*np.maximum(smooth,1e-12));emission=maximum_filter1d(emission.astype(int),size=7)>0;trace,td=refine_trace(s,v,valid,ref,wave,c)
            r['trace']=td;r['trace_y']=float(np.median(trace));r['date']=sci.header.get('DATE-OBS',r.get('date'));r['exptime']=sci.header.get('EXPTIME');r['shifta1']=sci.header.get('SHIFTA1');r['shifta2']=sci.header.get('SHIFTA2');r['sdqflags']=flags
            np.savez_compressed(f'products/individual_profiles/{r["root"]}_trace.npz',wave=wave,reference_trace=ref,measured_trace=trace,valid=valid)
            rr=np.floor(trace).astype(int)[None,:]+np.arange(-80,81)[:,None];inside=(rr>=0)&(rr<s.shape[0]);rr=np.clip(rr,0,s.shape[0]-1);xx=np.arange(s.shape[1])[None,:]
            np.savez_compressed(f'products/individual_profiles/{r["root"]}_native_pixels.npz',wave=wave,trace=trace,detector_row=rr,offset_pixel=rr-trace[None,:],counts=s[rr,xx],variance=v[rr,xx],valid=valid[rr,xx]&inside,dq=dq[rr,xx])
            # Diagnostic source and background regions in detector coordinates.
            from .plotting import detector_plot
            detector_plot(r,s,valid,trace,wave,c)
            for band in c['bands']:
                a,b=band;cols=np.where((wave>=a)&(wave<=b))[0]
                for lo,hi in c['exclude_wavelengths']:cols=cols[~((wave[cols]>=lo)&(wave[cols]<=hi))]
                if mask_spectral_emission:cols=cols[~emission[cols]]
                if len(cols)<12:continue
                for mode in ['aligned','unaligned','x2d']:
                    if mode=='x2d':
                        if 'x2d' not in paths:continue
                        with fits.open(paths['x2d']) as h2:
                            im=np.array(h2['SCI',ver].data,float);e=np.array(h2['ERR',ver].data,float);d=np.array(h2['DQ',ver].data,np.uint16);hh=h2['SCI',ver].header
                            w=hh['CRVAL1']+(np.arange(im.shape[1])+1-hh['CRPIX1'])*hh['CD1_1']
                            cc=np.where((w>=a)&(w<=b))[0]
                            for lo,hi in c['exclude_wavelengths']:cc=cc[~((w[cc]>=lo)&(w[cc]<=hi))]
                            if mask_spectral_emission:
                                for badwave in wave[emission]:cc=cc[abs(w[cc]-badwave)>0.3]
                            # Rectification maps along-slit position onto output rows; locate it directly.
                            pp=np.nansum(im[:,cc],axis=1);pos=int(np.argmax(pp));sl=slice(max(0,pos-7),min(len(pp),pos+8));zz=np.arange(len(pp))[sl];qq=pp[sl]-np.median(pp)
                            center=float(np.sum(zz*qq)/np.sum(qq));tr=np.full(im.shape[1],center)
                            good=((d&flags)==0)&np.isfinite(im)&np.isfinite(e)&(e>=0)
                            vv=e**2;positive=vv[vv>0];vv=np.maximum(vv,np.median(positive)*1e-6 if len(positive) else 1e-60)
                            y,p,cv,coverage=overlap_profile(im,vv,good,tr,cc,190)
                    else:y,p,cv,coverage=overlap_profile(s,v,valid,trace,cols,190,mode=='aligned')
                    for bi,(window,order) in enumerate(zip(c['background_windows'],c['background_orders'])):
                        if mode=='unaligned' and bi>0:continue
                        try:
                            sub,scov,bg,bmask,coeff,bcov=subtract_background(y,p,cv,coverage,window,order)
                            keep=abs(y)<=c['profile_radius'];yy=y[keep];qq=sub[keep];ccov=scov[np.ix_(keep,keep)]
                            m,norm,ncov=measures(yy,qq,ccov,c['aperture_heights'],wing_inner=c.get('wing_inner_pixels',10))
                            quality=[]
                            if m['snr']<c['core_snr']:quality.append('LOW_CORE_SNR')
                            if min(coverage[abs(y)<=8])<c['min_good_fraction']:quality.append('CORE_DQ_COVERAGE')
                            if np.min(coverage[keep])<.95:quality.append('WING_DQ_COVERAGE')
                            if td['rms'] is not None and td['rms']>c['trace_rms_limit']:quality.append('UNRELIABLE_TRACE')
                            if td['status']!='MEASURED_TRACE':quality.append('UNVERIFIED_TRACE')
                            item=dict(root=r['root'],target=r['target'],program=r['program'],aperture=r['aperture'],date=r['date'],trace_y=r['trace_y'],band=f'{a}-{b}',lo=a,hi=b,mode=mode,background=bi,ncolumns=len(cols),quality=quality,wing_usable=(not quality and m['wing_snr']>=c['wing_snr']),**m)
                            file=f'products/individual_profiles/{r["root"]}_v{ver}_{a}_{b}_{mode}_bg{bi}.npz';item['profile_file']=file
                            np.savez_compressed(file,y=yy,profile=qq,covariance=ccov,normalized=norm,normalized_covariance=ncov,full_y=y,unsubtracted=p,background=bg,background_mask=bmask,coverage=coverage,background_coeff=coeff,background_cov=bcov)
                            meta.append(item)
                        except Exception as e:logging.warning('%s %s %s bg%d: %s',r['root'],band,mode,bi,e)
    return meta


def extract(c):
    rows=load_json('data/catalogs/observations.json');allm=[]
    for r in rows:
        if not all(x in r.get('paths',{}) for x in ['flt','x1d']):continue
        file=Path('products/measurements')/(r['root']+'.json')
        try:
            cachekey=sha(r['paths']['flt'])+sha(r['paths']['x1d'])+sha(__file__)+json.dumps(c,sort_keys=True)
            if r.get('extraction_cachekey')==cachekey and file.exists():
                cached=load_json(file)
                if all(Path(x['profile_file']).exists() for x in cached):
                    allm.extend(cached);logging.info('Reused extracted %s',r['root']);continue
            m=extract_one(r,c);r['extraction_cachekey']=cachekey;save_json(file,m);allm.extend(m)
            good=[x for x in m if x['mode']=='aligned' and x['background']==0 and not [q for q in x['quality'] if q!='WING_DQ_COVERAGE']]
            r['processing_status']='EXTRACTED';r['psf_status']='CORE_USABLE' if good else 'REJECT_QA';r['measured_continuum_snr']=max((x['snr'] for x in good),default=0);r['wing_usable_bands']=sum(x['wing_usable'] for x in good)
            r['point_source_prior']=r['program'] in c['priority_programs']+[17804,18113] or r['target'].upper().startswith(('B1.5IV','B2V','O7V','WD','GD','HZ','GRW','G191','FEIGE','AGK','BD+'))
            logging.info('Extracted %s core bands %s wing bands %s',r['root'],len(good),r['wing_usable_bands'])
        except Exception as e:r['processing_status']='ERROR';r['processing_error']=str(e);logging.exception('Extraction %s',r['root'])
        save_json('products/measurements/exposure_qa.json',[{k:v for k,v in row.items() if k not in ['paths','status','preview_done','selected_for_preview','reason']} for row in rows if row.get('processing_status')])
    save_json('products/measurements/profiles.json',allm)
    return allm


def refresh_core_quality(c):
    rows=load_json('data/catalogs/observations.json')
    qa=load_json('products/measurements/exposure_qa.json',[]);byroot={r['root']:r for r in qa}
    m=load_json('products/measurements/profiles.json',[])
    for r in rows:
        mm=[a for a in m if a['root']==r['root'] and a['mode']=='aligned' and a['background']==0]
        if not mm:continue
        core=[a for a in mm if not [q for q in a['quality'] if q!='WING_DQ_COVERAGE']]
        update=byroot.setdefault(r['root'],{'root':r['root']})
        update['psf_status']='CORE_USABLE' if core else 'REJECT_QA'
        update['measured_continuum_snr']=max((a['snr'] for a in core),default=0)
        update['point_source_prior']=r['program'] in c['priority_programs']+[17804,18113] or r['target'].upper().startswith(('B1.5IV','B2V','O7V','WD','GD','HZ','GRW','G191','FEIGE','AGK','BD+'))
    save_json('products/measurements/exposure_qa.json',list(byroot.values()))
