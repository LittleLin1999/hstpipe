"""Persistent-source QA against independent, no-CR detector-image drizzles.

This is diagnostic only: STScI still owns production median/blot/CR rejection.
No detector pixels, reference WCS, or production CR parameters are fitted here.
"""
from pathlib import Path
import hashlib,json,warnings
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from ..common.target_runner import save
from ..common.provenance import identity,versions
from .astrometry.wcs import detector_wcs
from .drizzle import drizzle

VERSION='persistent-source-2'


def _sky(a,valid):
    return float(sigma_clipped_stats(a,mask=~valid)[1]) if valid.any() else np.nan


def aperture_assessment(stack,weights,science,weight,center,*,flux_tolerance=.02,width_ratio=1.05,min_weight_ratio=.8):
    """Compare like-grid aperture photometry, with empirical jackknife errors.

    Blank-field fluctuations and input jackknife scatter retain resampling
    covariance empirically; this is not an uncertainty map or precision claim.
    """
    yy,xx=np.indices(science.shape);rr=(xx-center[0])**2+(yy-center[1])**2
    ann=(rr>=12**2)&(rr<=17**2);core=rr<=2.5**2;valid=weights>0
    sky=np.array([_sky(a,ann&v) for a,v in zip(stack,valid)])
    samples=np.where(valid,stack-sky[:,None,None],np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning);reference=np.nanmedian(samples,axis=0)
    final_sky=_sky(science,ann&(weight>0));final=science-final_sky
    rows=[];failure=[];n=len(stack)
    for radius in (2.5,5.,10.):
        aperture=rr<=radius**2;area=int(aperture.sum());usable=(np.isfinite(reference)&(np.sum(valid,axis=0)>=max(2,n//2)))
        if not (usable&aperture).sum()>=area*.95 or np.mean(weight[aperture]>0)<.95:
            rows.append(dict(radius_pix=radius,status='insufficient_coverage'));failure.append('insufficient_source_coverage');continue
        fref=float(np.sum(reference[aperture]));fout=float(np.sum(final[aperture]));leave=[]
        for i in range(n):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning);med=np.nanmedian(np.delete(samples,i,axis=0),axis=0)
            leave.append(float(np.nansum(med[aperture])))
        jack=float(np.sqrt((n-1)/n*np.sum((np.asarray(leave)-np.mean(leave))**2))) if n>1 else np.inf
        # Empirical nearby same-area sums include pixel covariance. Restrict to
        # patches away from the actual source and reject contaminated high sums.
        blank=[]
        for angle in np.linspace(0,2*np.pi,16,endpoint=False):
            bx=center[0]+(18+radius)*np.cos(angle);by=center[1]+(18+radius)*np.sin(angle)
            q=(xx-bx)**2+(yy-by)**2<=radius**2
            if q.sum()>=area*.95 and np.mean(np.isfinite(reference[q]))>.95:blank.append(float(np.nansum(reference[q])))
        bnoise=float(sigma_clipped_stats(blank)[2]) if len(blank)>3 else 0.
        error=max(jack,bnoise,1e-12);snr=fref/error;delta=fref-fout;loss=delta/fref if fref>0 else np.nan
        failed=bool(fref>0 and delta>max(flux_tolerance*fref,3*error))
        if failed:failure.append('significant_persistent_source_flux_loss')
        # Positive excess can be residual CR, but median/mean and pixel-phase
        # differences require uncertainty-aware evidence, just as for deficits.
        excess=bool(fref>0 and -delta>max(flux_tolerance*fref,3*error))
        if excess:failure.append('significant_source_flux_excess')
        rows.append(dict(radius_pix=radius,status='review_required' if failed or excess else ('consistent' if snr>=3 else 'low_snr'),
            reference_flux=fref,final_flux=fout,loss_fraction=loss,empirical_uncertainty=error,reference_snr=snr,
            jackknife_uncertainty=jack,blank_aperture_uncertainty=bnoise,loss_significance=delta/error))
    outer=weight[(rr>=5**2)&(rr<=9**2)&(weight>0)]
    ratio=float(np.mean(weight[core])/np.median(outer)) if len(outer) else None
    if ratio is None or ratio<min_weight_ratio:failure.append('source_shaped_weight_deficit')
    def moments(a):
        signal=np.where((rr<=5**2)&np.isfinite(a),np.maximum(a,0),0);total=signal.sum()
        if total<=0:return None
        x=np.sum(xx*signal)/total;y=np.sum(yy*signal)/total
        cov=np.array([[np.sum((xx-x)**2*signal),np.sum((xx-x)*(yy-y)*signal)],
                      [np.sum((xx-x)*(yy-y)*signal),np.sum((yy-y)**2*signal)]])/total
        e=np.linalg.eigvalsh(cov)
        return dict(centroid=[float(x),float(y)],fwhm_pix=float(2.355*np.sqrt(np.mean(e))),ellipticity=float(1-np.sqrt(max(e[0],0)/max(e[1],1e-20))))
    before=moments(reference);after=moments(final)
    high_snr=any(x.get('reference_snr',0)>=20 and x.get('radius_pix')==5 for x in rows)
    shape_samples=[]
    for i in range(n):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning);sample=np.nanmedian(np.delete(samples,i,axis=0),axis=0)
        q=moments(sample)
        if q:shape_samples.append(q)
    width_error=float(np.std([q['fwhm_pix'] for q in shape_samples])) if shape_samples else None
    if before and after and high_snr:
        if after['fwhm_pix']-before['fwhm_pix']>max(before['fwhm_pix']*(width_ratio-1),3*(width_error or 0.)):failure.append('persistent_source_broadening')
        if np.linalg.norm(np.subtract(after['centroid'],before['centroid']))>.5:failure.append('source_centroid_inconsistency')
    return dict(status='validated' if not failure else 'review_required',apertures=rows,core_weight_ratio=ratio,
        before=before,after=after,fwhm_empirical_uncertainty_pix=width_error,failures=sorted(set(failure)),uncertainty_method='Input-image jackknife and empirical same-grid blank apertures; no pixelwise ERR or independent accuracy claim'),reference


def _common_sources(inputs,relative,output_wcs):
    grouped={}
    for entry in relative.get('detector_catalogs',[]):
        p=next((p for p in inputs if p.name==entry['exposure']),None)
        if p is None or not Path(entry['path']).exists():continue
        cat=Table.read(entry['path'],format='ascii.ecsv')
        if 'common_source_id' not in cat.colnames:continue
        with fits.open(p) as hd:
            for row in cat:
                if 'centroid_valid' in cat.colnames and not row['centroid_valid']:continue
                v=int(row['extver']);sky=detector_wcs(hd,v).all_pix2world([[row['x'],row['y']]],0)[0]
                xy=output_wcs.all_world2pix([sky],0)[0]
                grouped.setdefault(str(row['common_source_id']),[]).append(dict(exposure=p.name,extver=v,x=float(row['x']),y=float(row['y']),output_x=float(xy[0]),output_y=float(xy[1])))
    return {sid:rows for sid,rows in grouped.items() if len({r['exposure'] for r in rows})>=2}


def validate_source_rejection(product,inputs,paths,relative,drizzle_options,*,refresh=False):
    """Restartable QA of an existing final pass; does not redo that reduction."""
    product=Path(product);inputs=sorted(map(Path,inputs),key=lambda p:p.name);output=paths.file('intermediate','drizzle','source_validation')
    if len(inputs)<2:return dict(status='not_applicable_single_exposure',limitations=['No independent multi-exposure CR/source check'])
    with fits.open(product) as hd:
        w=WCS(hd['SCI'].header,hd);sources=_common_sources(inputs,relative,w)
        digest=hashlib.sha256(hd['SCI'].data.tobytes());digest.update(str(w.to_header(relax=True)).encode())
    if not sources:return dict(status='unvalidated',reason='No saved repeatable HST source identities for independent CR/source assessment')
    dependencies=[Path(__file__),Path(__file__).with_name('drizzle.py')]+[Path(c['path']) for c in relative.get('detector_catalogs',[]) if Path(c['path']).exists()]
    for p in inputs:
        root=p.stem.removesuffix('_flc')
        for dependency in (paths.exposure_file('intermediate','working',root,'flc'),paths.exposure_file('intermediate','alignment',root,'absolute_hdrlet')):
            if dependency.exists():dependencies.append(dependency)
    masks=[paths.exposure_file('intermediate','drizzle',p.stem.removesuffix('_flc'),'final_crmask') for p in inputs]
    cache=paths.cache();key=cache.key(dependencies+masks,dict(version=VERSION,product=digest.hexdigest(),sources=sources,parameters=drizzle_options,versions=versions()))
    if not refresh and cache.valid('uvis_source_validation',key):return dict(json.loads(output.read_text()),restored_from_cache=True)
    cutfile=paths.file('intermediate','drizzle','source_validation_cutouts','fits');stamps=[fits.PrimaryHDU()];measurements=[];panels=[]
    stamps[0].header['QAONLY']=True;stamps[0].header['QAVERS']=VERSION;stamps[0].header['REFIMAGE']=str(product)
    with paths.scratch('source-qa-') as tmp:
        temp=Path(tmp);copies=[]
        for p in inputs:
            dest=temp/p.name;paths.copy(p,dest);copies.append(dest)
        import os
        previous=Path.cwd()
        try:
            os.chdir(temp)
            drizzle(copies,temp/paths.group_id,dict(static=False,driz_separate=True,median=False,blot=False,driz_cr=False,driz_combine=False,
                driz_sep_wcs=True,driz_sep_refimage=f'{product}[SCI,1]',driz_sep_pixfrac=drizzle_options.get('final_pixfrac',1.),
                driz_sep_kernel=drizzle_options.get('final_kernel','square'),driz_sep_bits=drizzle_options.get('final_bits',80),runfile='source_qa_no_cr.log'))
        finally:os.chdir(previous)
        singles=[]
        for p in copies:
            root=p.stem.removesuffix('_flc');sci=temp/(root+'_single_sci.fits');wht=temp/(root+'_single_wht.fits')
            if not sci.exists():raise FileNotFoundError(f'Native separate drizzle did not produce {sci.name}')
            singles.append((sci,wht))
        with fits.open(product) as hd:
            for sid,observations in sources.items():
                x=float(np.median([q['output_x'] for q in observations]));y=float(np.median([q['output_y'] for q in observations]));ix,iy=round(x),round(y);radius=38
                from astropy.nddata import Cutout2D
                def cut(array,fill=np.nan):
                    return np.asarray(Cutout2D(np.asarray(array,dtype=float),(ix,iy),(2*radius+1,2*radius+1),mode='partial',fill_value=fill).data,float)
                stack=[];weights=[]
                for i,(sci,wht) in enumerate(singles):
                    with fits.open(sci) as sh,fits.open(wht) as wh:
                        a=cut(sh[0].data);weight=cut(wh[0].data,0)
                        unit=sh[0].header.get('BUNIT','').upper()
                        if unit in ('ELECTRONS','COUNTS'):a/=fits.getheader(inputs[i])['EXPTIME']
                        stack.append(a);weights.append(weight)
                final=cut(hd['SCI'].data);weight=cut(hd['WHT'].data,0)
                result,reference=aperture_assessment(np.array(stack),np.array(weights),final,weight,(x-ix+radius,y-iy+radius))
                result.update(source_id=sid,observations=observations,output_x=x,output_y=y);measurements.append(result)
                for j,a in enumerate(stack):
                    h=fits.ImageHDU(np.asarray(a,'f4'),name='INPUT');h.header['EXTVER']=len(stamps);h.header['SOURCEID']=sid;h.header['EXPOSURE']=inputs[j].name;h.header['BUNIT']='ELECTRONS/S';stamps.append(h)
                for name,a in (('MEDIAN',reference),('FINAL',final),('WHT',weight)):
                    h=fits.ImageHDU(np.asarray(a,'f4'),name=name);h.header['EXTVER']=len(stamps);h.header['SOURCEID']=sid;h.header['BUNIT']=hd['WHT'].header.get('BUNIT','WEIGHT') if name=='WHT' else 'ELECTRONS/S';stamps.append(h)
                panels.append((sid,reference,final,weight,result))
        paths.copy(temp/'source_qa_no_cr.log',paths.file('logs',None,'source_qa_no_cr','log'))
    valid=[r for r in measurements if r['status']=='validated'];failed=[r for r in measurements if r['status']=='review_required']
    result=dict(version=VERSION,status='validated' if valid and not failed else 'review_required',sources=measurements,
        validated_sources=len(valid),failed_sources=len(failed),independent_no_cr_inputs=[p.name for p in inputs],
        method='STScI separate square drizzle on the exact final WCS, followed by diagnostic pixel-median, aperture, PSF and WHT comparisons',
        original_pixels_modified=False,production_cr_parameters_modified=False,cutouts=str(cutfile),
        limitations=['Median/mean and detector pixel-phase differences enter empirical uncertainty; this is not calibrated aperture photometry or a covariance map'])
    fits.HDUList(stamps).writeto(cutfile,overwrite=True,checksum=True)
    _plot(panels,paths)
    save(output,result);cache.record('uvis_source_validation',key,[output,cutfile,paths.file('diagnostics','drizzle','source_validation','png')])
    return result


def _plot(panels,paths):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm
    if not panels:return
    paths.area('diagnostics',create=True)
    fig,axes=plt.subplots(len(panels),3,figsize=(10,3*len(panels)),squeeze=False)
    for row,(sid,reference,final,weight,result) in zip(axes,panels):
        for ax,a,title in zip(row,(reference,final,weight),('Independent no-CR median','Final SCI','Final WHT')):
            a=a[26:51,26:51]
            if np.isfinite(a).any():ax.imshow(a,origin='lower',cmap='gray',norm=simple_norm(a[np.isfinite(a)],'asinh',percent=99.5))
            else:ax.text(.5,.5,'No coverage',ha='center')
            ax.set_title(title+'\n'+sid,fontsize=8)
        row[0].set_ylabel(result['status'])
    fig.tight_layout();fig.savefig(paths.file('diagnostics','drizzle','source_validation','png'),dpi=105);plt.close(fig)


def publish_validated_image(product,paths,assessment):
    """Publish an assessed image without changing SCI/WHT/CTX or the WCS."""
    from .qa import validate_drizzle
    if assessment.get('status')!='validated':raise ValueError('Persistent-source validation must pass before publication')
    product=Path(product);validate_drizzle(product)
    with fits.open(product) as hd:
        if hd['SCI'].header.get('ASTSTAT')!='VALIDATED':raise ValueError('Absolute astrometry is not validated')
        if hd['SCI'].header.get('ASTREL') not in ('validated','single_exposure'):raise ValueError('Relative registration is not validated')
        before={h.name:hashlib.sha256(h.data.tobytes()).hexdigest() for h in hd if h.name in ('SCI','WHT','CTX')}
    destination=paths.product('drc')
    with paths.scratch('publish-') as tmp:
        candidate=Path(tmp)/destination.name;paths.copy(product,candidate)
        with fits.open(candidate,mode='update') as hd:
            hd[0].header['HPQASTAT']=('VALIDATED','Astrometry and persistent-source image QA passed')
            hd[0].header['HPSRCQA']=(assessment.get('version',VERSION),'Persistent-source CR/image assessment')
            hd[0].header['HPNQASRC']=(assessment.get('validated_sources',0),'Repeatable astronomical sources assessed')
            hd[0].header.add_history('Source/CR QA: independent no-CR input drizzles, empirical aperture/shape checks and WHT coverage; raw bright peaks are not classified as stars')
            for h in hd:
                if h.name in before and hashlib.sha256(h.data.tobytes()).hexdigest()!=before[h.name]:raise ValueError('Image arrays changed during QA annotation')
                h.add_checksum()
        validate_drizzle(candidate)
        if destination.exists():
            from datetime import datetime,timezone
            recovery=paths.recovery_directory('science-before-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))
            recovery.mkdir(parents=True);paths.copy(destination,recovery/destination.name)
        paths.copy(candidate,destination)
    return destination
