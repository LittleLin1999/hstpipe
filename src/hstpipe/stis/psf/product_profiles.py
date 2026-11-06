"""Measure native and drizzled PSFs from hstpipe reduction products only."""
from pathlib import Path
import logging
import numpy as np
from astropy.io import fits
from astropy.table import Table
from .common import *
from .extraction import overlap_profile,measures


def wavelength_weights(wavelength,lo,hi,excluded):
    """Fractional detector-pixel membership in a wavelength interval."""
    dispersion=np.abs(np.gradient(wavelength,axis=1))
    left=wavelength-dispersion/2;right=wavelength+dispersion/2
    weight=np.clip((np.minimum(right,hi)-np.maximum(left,lo))/dispersion,0,1)
    merged=[]
    for a,b in sorted(excluded):
        if merged and a<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],b)
        else:merged.append([a,b])
    for a,b in merged:
        aa=max(lo,a);bb=min(hi,b)
        if bb>aa:weight-=np.clip((np.minimum(right,bb)-np.maximum(left,aa))/dispersion,0,1)
    return np.clip(weight,0,1)


def feature_intervals(path,factor):
    from scipy.ndimage import median_filter,maximum_filter1d
    with fits.open(path) as h:
        a=h[1].data[0];w=np.asarray(a['WAVELENGTH'],float)
        if h[0].header.get('HELCORR')=='COMPLETE':w=w/factor
        net=np.asarray(a['NET'],float);err=np.asarray(a['NET_ERROR'],float);smooth=median_filter(net,size=31)
        flag=(err>0)&(net>smooth+5*err)&(net>2*np.maximum(smooth,1e-12));flag=maximum_filter1d(flag.astype(int),size=7)>0
    edges=np.diff(np.r_[False,flag,False].astype(int));starts=np.where(edges==1)[0];ends=np.where(edges==-1)[0]-1;dw=np.median(np.diff(w))
    return [[float(w[i]-dw/2),float(w[j]+dw/2)] for i,j in zip(starts,ends)]


def write_profile(item,y,p,cov,norm,ncov,coverage,background=None):
    stem=safe(f'{item["root"]}_{item["lo"]}_{item["hi"]}_{item["mode"]}')
    directory=Path('psf_library/native' if item['mode']=='hstpipe_native' else 'psf_library/drizzled');directory.mkdir(parents=True,exist_ok=True)
    path=directory/(stem+'.fits');npz=directory/(stem+'.npz')
    np.savez_compressed(npz,y=y,profile=p,covariance=cov,normalized=norm,normalized_covariance=ncov,coverage=coverage)
    header=fits.Header()
    for key,value in dict(ROOTNAME=item['root'],TARGNAME=item['target'],APERTURE=item['aperture'],PHYSAPER=item.get('physical_aperture',item['aperture']),WAVEMIN=item['lo'],WAVEMAX=item['hi'],SAMPLING=item['mode'],NORMRAD=item.get('normalization_radius',70),TRACEXY=item['trace_y'],DATEOBS=str(item['date']),GRADE='WING' if item['wing_usable'] else 'CORE',QFLAGS=','.join(item['quality']) or 'PASS').items():header[key]=value
    table=Table(dict(offset_pixel=y,offset_arcsec=y*.0246,flux=p,error=np.sqrt(np.maximum(np.diag(cov),0)),normalized=norm,normalized_error=np.sqrt(np.maximum(np.diag(ncov),0)),coverage=coverage))
    table['offset_pixel'].unit='pix';table['offset_arcsec'].unit='arcsec'
    header['BINWIDTH']=float(np.median(np.diff(y)));header['PIXSCAL']=.0246
    header['FLUXUNIT']='summed calibrated surface-brightness samples';header['WAVEFRAM']=item.get('wave_frame','OBSERVED')
    hd=[fits.PrimaryHDU(header=header),fits.BinTableHDU(table,name='PROFILE'),fits.ImageHDU(cov,name='COVAR'),fits.ImageHDU(ncov,name='NORMCOV')]
    if background is not None:hd.append(fits.ImageHDU(background,name='BKG'))
    fits.HDUList(hd).writeto(path,overwrite=True,checksum=True)
    item['profile_file']=str(npz);item['profile_fits']=str(path)
    return item


def measure_products(c):
    reductions=load_json('data/catalogs/hstpipe_reductions.json',[]);rows=catalog_apertures(c);byroot={r['root']:r for r in rows};out=[]
    for red in reductions:
        directory=Path(red.get('directory',''))
        if not (directory/'intermediate/background').is_dir():continue
        for root in red['roots']:
            r=byroot[root];bgfile=directory/'intermediate/background'/(root+'_flt_rate.fits')
            if not bgfile.exists():continue
            try:
                d=np.load(f'products/individual_profiles/{root}_trace.npz');wave=d['wave'];trace=d['measured_trace']
                with np.load(directory/'intermediate/flt/registration'/(root+'_flt.npz')) as op:
                    wavelength=op['observed_wavelength'];factor=op['factor'];mapping_valid=op['valid']
                control=Path('data/controls')/(root+'_flt_perform.fits');vhelio=fits.getheader(control,'SCI')['V_HELIO'];helio_factor=1-vhelio/299792.458
                excluded=c['exclude_wavelengths']+feature_intervals(r['paths']['x1d'],helio_factor)
                with fits.open(bgfile) as h:
                    t=float(h['SCI'].header['EXPTIME']);image=h['SCI'].data.astype(float)*t;var=(h['ERRBASE'].data.astype(float)*t)**2;bge=h['BKGERR'].data.astype(float)*t;dq=h['DQ'].data.astype(np.uint16);valid=(dq==0)&np.isfinite(image)&np.isfinite(var)&mapping_valid;bkg=h['BKG'].data.astype(float)*t
                image*=factor;var*=factor**2;bge*=factor;bkg*=factor
                for lo,hi in c['bands']:
                    spectral_weights=wavelength_weights(wavelength,lo,hi,excluded)
                    xx=np.arange(image.shape[1]);rr=np.floor(trace).astype(int)[None,:]+np.arange(-72,73)[:,None];rr=np.clip(rr,0,image.shape[0]-1)
                    cols=np.where(np.any(spectral_weights[rr,xx[None,:]]>0,axis=0))[0]
                    candidate_columns=len(cols);original_columns=cols.copy();radius=None
                    for candidate_radius in c.get('normalization_radii',[70,60,50,40,30,21,15,10,8]):
                        rr=np.floor(trace[original_columns]).astype(int)[None,:]+np.arange(-candidate_radius-2,candidate_radius+3)[:,None]
                        inside=(rr>=0)&(rr<valid.shape[0]);jj=np.clip(rr,0,valid.shape[0]-1)
                        complete=np.all(inside&valid[jj,original_columns[None,:]],axis=0)
                        if complete.sum()>=max(c.get('minimum_complete_columns',12),c.get('min_complete_column_fraction',.25)*candidate_columns):
                            cols=original_columns[complete];radius=candidate_radius;break
                    if radius is None:
                        logging.warning('%s %s-%s: no complete native PSF window',root,lo,hi);continue
                    y,p,cov,coverage=overlap_profile(image*spectral_weights,var*spectral_weights**2,valid,trace,cols,radius)
                    _,be,_,_=overlap_profile(bge*spectral_weights,np.zeros_like(var),valid,trace,cols,radius)
                    _,bp,_,_=overlap_profile(bkg*spectral_weights,np.zeros_like(var),valid,trace,cols,radius)
                    cov+=np.outer(be,be) # hstpipe BKGERR is a correlated bound.
                    try:m,norm,ncov=measures(y,p,cov,c['aperture_heights'],wing_inner=c.get('wing_inner_pixels',10))
                    except ValueError as error:
                        logging.warning('%s native %s-%s rejected band: %s',root,lo,hi,error);continue
                    quality=[]
                    for height in c['aperture_heights']:
                        if height>2*radius:
                            m[f'ee_{height}']=float('nan');m[f'ee_{height}_error']=float('nan')
                    if red.get('status')!='COMPLETE':quality.append('INCOMPLETE_REDUCTION')
                    if m['snr']<c['core_snr']:quality.append('LOW_CORE_SNR')
                    if min(coverage[abs(y)<=8])<c['min_good_fraction']:quality.append('CORE_DQ_COVERAGE')
                    if min(coverage)<.95:quality.append('WING_DQ_COVERAGE')
                    item=dict(root=root,roots=[root],target=r['target'],program=r['program'],aperture=r.get('configuration_aperture',r['aperture']),physical_aperture=r['aperture'],date=r['date'],trace_y=r['trace_y'],lo=lo,hi=hi,band=f'{lo}-{hi}',wing_inner_radius=c.get('wing_inner_pixels',10),mode='hstpipe_native',wave_frame='OBSERVED',background=0,normalization_radius=radius,ncolumns=len(cols),candidate_columns=candidate_columns,dq_column_fraction=1-len(cols)/candidate_columns,quality=quality,wing_usable=not quality and radius>=c.get('minimum_wing_radius',40) and m['wing_snr']>=c['wing_snr'],reduction_group=red['group'],input_background=str(bgfile),**m)
                    out.append(write_profile(item,y,p,cov,norm,ncov,coverage,bp))
            except Exception:logging.exception('hstpipe native PSF %s',root)
        if red.get('status')!='COMPLETE':continue
        product=next((p for p in red['products'] if p.endswith('_s2d_drz.fits')),None)
        if not product:continue
        try:
            r=byroot[red['roots'][0]]
            with fits.open(product) as h:
                hh=h['SCI'].header;im=h['SCI'].data.astype(float);err=h['ERR'].data.astype(float);dq=h['DQ'].data.astype(np.uint16);bge=h['BKGERR'].data.astype(float) if 'BKGERR' in h else np.zeros_like(im)
                wave=hh['CRVAL1']+(np.arange(im.shape[1])+1-hh['CRPIX1'])*hh['CD1_1'];y=(hh['CRVAL2']+(np.arange(im.shape[0])+1-hh['CRPIX2'])*hh['CD2_2'])/.0246
                helio_factors=[1-fits.getheader(Path('data/controls')/(x+'_flt_perform.fits'),'SCI')['V_HELIO']/299792.458 for x in red['roots']]
                helio_factor=float(np.average(helio_factors,weights=[byroot[x].get('exptime',1) or 1 for x in red['roots']]))
                wave=wave/helio_factor;feature_exclusions=[]
                for x,f in zip(red['roots'],helio_factors):feature_exclusions.extend(feature_intervals(byroot[x]['paths']['x1d'],f))
                spatial_step=hh['CD2_2'];factor=max(1.,c['hstpipe_pixfrac']*.0246/spatial_step)
                for lo,hi in c['bands']:
                    column_weights=wavelength_weights(wave[None,:],lo,hi,c['exclude_wavelengths']+feature_exclusions)[0]
                    cols=np.where(column_weights>0)[0];cutweights=column_weights[cols][None,:]
                    if len(cols)<12:continue
                    native_radii=[r['normalization_radius'] for r in out if r['mode']=='hstpipe_native' and r.get('reduction_group')==red['group'] and r['lo']==lo and r['hi']==hi]
                    radius=min(native_radii) if native_radii else 40
                    yy=abs(y)<=radius;yyi=np.where(yy)[0];good=(dq[np.ix_(yyi,cols)]==0)&np.isfinite(im[np.ix_(yyi,cols)]);exposure_coverage=h['NEXP'].data[np.ix_(yyi,cols)].astype(float)/len(red['roots']) if 'NEXP' in h else np.ones_like(good,dtype=float)
                    coverage=np.mean(np.where(good,np.clip(exposure_coverage,0,1),0),axis=1)
                    p=np.sum(np.where(good,im[np.ix_(yyi,cols)]*cutweights,0),axis=1)
                    # ERR includes background; separate that bound before applying a
                    # conservative spatial drizzle correlation factor to ERRBASE.
                    base=h['ERRBASE'].data[np.ix_(yyi,cols)].astype(float) if 'ERRBASE' in h else err[np.ix_(yyi,cols)]
                    vv=np.sum(np.where(good,base**2*cutweights**2,0),axis=1)*factor;be=np.sum(np.where(good,bge[np.ix_(yyi,cols)]*cutweights,0),axis=1);cov=np.diag(vv)+np.outer(be,be)
                    heights=[v*.0246/spatial_step for v in c['aperture_heights']]
                    # Measurement helper uses actual profile-bin coordinates for width;
                    # aperture measurements are recomputed with physical bin overlaps.
                    yp=y[yy]
                    try:m,norm,ncov=measures(yp,p,cov,c['aperture_heights'],wing_inner=c.get('wing_inner_pixels',10))
                    except ValueError as error:
                        logging.warning('%s drizzled %s-%s rejected band: %s',red['group'],lo,hi,error);continue
                    for height in c['aperture_heights']:
                        half=spatial_step/.0246/2;a=np.clip((np.minimum(yp+half,height/2)-np.maximum(yp-half,-height/2))/(2*half),0,1)
                        m[f'ee_{height}']=float(a@norm) if height<=2*radius else float('nan');m[f'ee_{height}_error']=float(np.sqrt(max(a@ncov@a,0))) if height<=2*radius else float('nan')
                    quality=[]
                    if not native_radii:quality.append('NATIVE_PROFILE_UNAVAILABLE')
                    if red.get('status')!='COMPLETE':quality.append('INCOMPLETE_REDUCTION')
                    if m['snr']<c['core_snr']:quality.append('LOW_CORE_SNR')
                    if min(coverage[abs(yp)<min(8,radius)])<c['min_good_fraction']:quality.append('CORE_DQ_COVERAGE')
                    item=dict(root=red['roots'][0],roots=red['roots'],target=r['target'],program=r['program'],aperture=r.get('configuration_aperture',r['aperture']),physical_aperture=r['aperture'],date=r['date'],trace_y=float(np.mean([byroot[x]['trace_y'] for x in red['roots']])),lo=lo,hi=hi,band=f'{lo}-{hi}',wing_inner_radius=c.get('wing_inner_pixels',10),mode='drizzled',normalization_radius=radius,background=0,quality=quality,wing_usable=not quality and radius>=c.get('minimum_wing_radius',40) and m['wing_snr']>=c['wing_snr'],reduction_group=red['group'],wave_frame='OBSERVED',drizzle_correlation_factor=factor,input_product=product,**m)
                    out.append(write_profile(item,yp,p,cov,norm,ncov,coverage))
        except Exception:logging.exception('Drizzled PSF %s',red['group'])
    old=load_json('products/measurements/profiles.json',[]);old=[r for r in old if r['mode'] not in ['hstpipe_native','drizzled']]
    save_json('products/measurements/profiles.json',old+out);save_json('products/measurements/hstpipe_psfs.json',out)
    return out
