"""Resolve stellar identity using native HST pixels, never catalog-fit coordinates.

Aperture concentration is the default compact-source profile diagnostic.
Measured HST centroids retain shared pixel QC; concentration is not a stellar
classification by itself. Optional Gaussian PRF fitting remains available. Multiple credible stars remain ambiguous until an HST
reference profile discriminates them; distance/brightness never select a winner.
"""
from pathlib import Path
import json
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from .centroids import pixel_masks,validate_centroids
from ..fits import image_data
from .wcs import detector_wcs

DEFAULTS=dict(enabled=True,profile_method='aperture',aperture_radii_fwhm=[.5,2.],
    min_aperture_valid_fraction=.9,profile_concentration_fraction=.2,fit_shape=7,check_shape=9,xy_bounds=2.,fwhm_bounds_factors=[.55,1.8],
    max_qfit=.5,max_centroid_error_pix=.3,max_window_shift_pix=.2,duplicate_radius_pix=.6,
    group_separation_pix=5.,profile_flux_fraction=.2,profile_width_fraction=.2)

def stellar_options(options):return dict(DEFAULTS,**options.get('psf_resolution',{}))

def fit_candidates(image,catalog,version,detection,options,contexts=None):
    method=stellar_options(options)['profile_method']
    if method=='aperture':return aperture_candidates(image,catalog,version,detection,options,contexts)
    if method=='psf':return psf_candidates(image,catalog,version,detection,options,contexts)
    raise ValueError('profile_method must be aperture or psf')


def aperture_candidates(image,catalog,version,detection,options,contexts=None):
    """Native-pixel concentric aperture fluxes; never replace HST centroids.

    Fractional pixel weights come from Photutils. ERR propagation includes the
    covariance of nested apertures and their shared local background. Masked
    pixels are excluded, never filled or area-corrected to invent core flux.
    """
    from photutils.aperture import CircularAperture
    opts=stellar_options(options);version=int(version);rows=[]
    contexts=contexts if contexts is not None else {}
    with fits.open(image) as hd:
        if version not in contexts:
            data=np.asarray(hd['SCI',version].data,float);dq,base,mask=pixel_masks(hd,version,detection)
            contexts[version]=dict(data=data,dq=dq,base=base,mask=mask,error=np.asarray(image_data(hd['ERR',version]),float),sky=0.,rms=1.)
        ctx=contexts[version];ctx.setdefault('sky',0.);ctx.setdefault('rms',1.)
        a=ctx['data'];error=ctx['error'];mask=ctx['mask']
        scale=np.sqrt(abs(np.linalg.det(detector_wcs(hd,version).pixel_scale_matrix)))*3600
        fwhm=detection.get('fwhm_arcsec',.1)/scale if detection.get('fwhm_arcsec') else detection.get('fwhm_pixels',2.5)
        radii=np.asarray(opts['aperture_radii_fwhm'],float)*fwhm
        if len(radii)!=2 or not (0<radii[0]<radii[1]):raise ValueError('Require positive increasing inner/outer aperture radii')
        initial=np.column_stack([catalog['x'],catalog['y']]).astype(float)
        for index,row in enumerate(catalog):
            x,y=initial[index];entry=dict(index=index,source_id=str(row['hst_source_id']),extver=version,status='rejected')
            if not np.isfinite([x,y]).all():
                rows.append(dict(entry,reason='nonfinite_centroid'));continue
            radius=max(13,int(np.ceil(3*radii[1])));ix=int(round(x));iy=int(round(y))
            sx=slice(max(0,ix-radius),min(a.shape[1],ix+radius+1));sy=slice(max(0,iy-radius),min(a.shape[0],iy+radius+1))
            data=a[sy,sx];err=error[sy,sx];bad=mask[sy,sx]|~np.isfinite(err)|(err<=0)
            yy,xx=np.indices(data.shape);blank=~bad.copy()
            for bx,by in initial:
                if np.isfinite([bx,by]).all():blank&=np.hypot(xx+sx.start-bx,yy+sy.start-by)>max(7.,radii[1])
            if blank.sum()<30:
                rows.append(dict(entry,reason='insufficient_local_background'));continue
            _,sky,noise=sigma_clipped_stats(data[blank])
            # Conservative median-background error, including read noise from ERR.
            skyvar=np.pi/2*max(float(noise)**2,float(np.median(err[blank]**2)))/blank.sum()
            weights=[];coverage=[]
            for radius in radii:
                ap=CircularAperture((x-sx.start,y-sy.start),r=float(radius))
                w=ap.to_mask(method='exact').to_image(data.shape)
                if w is None:w=np.zeros(data.shape)
                w=np.where(bad,0.,w);coverage.append(float(w.sum()/ap.area));weights.append(w)
            wi,wo=weights;safe=np.where(bad,0.,data-float(sky));var=np.where(bad,0.,err**2)
            fi=float(np.sum(wi*safe));fo=float(np.sum(wo*safe))
            vi=float(np.sum(wi**2*var)+wi.sum()**2*skyvar);vo=float(np.sum(wo**2*var)+wo.sum()**2*skyvar)
            cov=float(np.sum(wi*wo*var)+wi.sum()*wo.sum()*skyvar)
            concentration=fi/fo if fo>0 else np.nan
            ratioerr=np.sqrt(max(0.,vi+concentration**2*vo-2*concentration*cov))/fo if fo>0 else np.nan
            snr=fo/np.sqrt(vo) if vo>0 else np.nan
            raw=Table(dict(x=[x],y=[y],extver=[version],flux=[fo],snr=[snr],local_sky=[float(sky)],
                centroid_error_pix=[float(row['centroid_error_pix']) if 'centroid_error_pix' in row.colnames else np.nan]))
            clean,audit=validate_centroids(image,raw,detection,method='aperture_profile',contexts={version:ctx},return_all=True)
            q=audit[0];reasons=[]
            if not len(clean):reasons.append('pixel_qc:'+str(q['quality_flags']))
            if min(coverage)<opts['min_aperture_valid_fraction']:reasons.append('insufficient_aperture_coverage')
            if not np.isfinite([fi,fo,vi,vo,concentration,ratioerr]).all() or fi<=0 or fo<=0 or vo<=0:reasons.append('invalid_aperture_flux')
            if float(q['centroid_error_pix'])>opts['max_centroid_error_pix']:reasons.append('poor_centroid_precision')
            if float(q['fwhm_pix'])>opts['fwhm_bounds_factors'][1]*fwhm:reasons.append('extended_profile')
            entry.update(initial_x=float(x),initial_y=float(y),x=float(x),y=float(y),fwhm_pix=float(q['fwhm_pix']),
                width_pix=float(q['width_pix']),ellipticity=float(q['ellipticity']),centroid_error_pix=float(q['centroid_error_pix']),
                snr=float(snr),flux=fo,flux_error=float(np.sqrt(vo)),flux_rate=fo/float(hd[0].header['EXPTIME']),
                inner_flux=fi,inner_flux_error=float(np.sqrt(vi)),concentration=float(concentration),concentration_error=float(ratioerr),
                aperture_radii_pix=radii.tolist(),aperture_valid_fraction=coverage,local_sky=float(sky),local_noise=float(noise),
                core_valid_fraction=float(q['core_valid_fraction']),window_stability_pix=float(q['centroid_stability_pix']),
                original_detection_method=str(row['detection_method']) if 'detection_method' in row.colnames else 'custom',
                detection_method='aperture_profile',profile_class='compact_candidate',reason=';'.join(reasons),
                status='rejected' if reasons else 'credible_stellar_fit')
            rows.append(entry)
    good=[]
    for r in sorted((r for r in rows if r['status']=='credible_stellar_fit'),key=lambda r:r['centroid_error_pix']):
        duplicate=next((q for q in good if np.hypot(q['x']-r['x'],q['y']-r['y'])<opts['duplicate_radius_pix']),None)
        if duplicate:r.update(status='duplicate_fit',reason='Duplicate native HST detection',duplicate_of=duplicate['source_id'])
        else:good.append(r)
    # Simple apertures do not deblend distinct neighboring compact sources.
    contaminated=set()
    for i,r in enumerate(good):
        for j,q in enumerate(good[:i]):
            if np.hypot(r['x']-q['x'],r['y']-q['y'])<radii[1]:contaminated.update((i,j))
    for i in contaminated:good[i].update(status='rejected',reason='neighbor_inside_outer_aperture; deblending requires optional PSF assessment')
    good=[r for r in good if r['status']=='credible_stellar_fit']
    for r in rows:
        for k,v in list(r.items()):
            if isinstance(v,(float,np.floating)) and not np.isfinite(v):r[k]=None
    return good,dict(status='unique_compact_profile' if len(good)==1 else ('ambiguous_stellar_profiles' if good else 'no_reliable_stellar_profile'),
        raw_candidates=len(catalog),duplicates=sum(r['status']=='duplicate_fit' for r in rows),rejected=sum(r['status']=='rejected' for r in rows),fits=rows,
        model='Photutils concentric aperture concentration; original HST centroids; no Gaussian profile residual gate',
        uncertainty_method='ERR weighted-sum covariance, nested apertures and shared local-background uncertainty; existing centroid QC',options=opts)


def psf_candidates(image,catalog,version,detection,options,contexts=None):
    """Fit supplied HST detections jointly when overlapping, then shared pixel QC.

    All coordinates are zero-origin native detector centroids. No external sky
    coordinate is accepted by this interface. Invalid source cores stay masked.
    """
    from photutils.psf import PSFPhotometry,CircularGaussianPRF,SourceGrouper
    from astropy.modeling.fitting import TRFLSQFitter
    opts=stellar_options(options);contexts=contexts if contexts is not None else {}
    if not len(catalog):return [],dict(status='no_candidates',raw_candidates=0,fits=[])
    version=int(version)
    with fits.open(image) as hd:
        if version not in contexts:
            a=np.asarray(hd['SCI',version].data,float);dq,base,mask=pixel_masks(hd,version,detection)
            err=np.asarray(image_data(hd['ERR',version]),float)
            contexts[version]=dict(data=a,dq=dq,base=base,mask=mask,error=err,sky=0.,rms=1.)
        ctx=contexts[version];ctx.setdefault('sky',0.);ctx.setdefault('rms',1.)
        a=ctx['data'];mask=ctx['mask'];error=ctx['error']
        initial=np.column_stack([catalog['x'],catalog['y']]).astype(float)
        groups=SourceGrouper(opts['group_separation_pix'])(initial[:,0],initial[:,1])
        rows=[]
        scale=np.sqrt(abs(np.linalg.det(detector_wcs(hd,version).pixel_scale_matrix)))*3600
        fwhm=detection.get('fwhm_arcsec',.1)/scale if detection.get('fwhm_arcsec') else detection.get('fwhm_pixels',2.5)
        bounds=np.asarray(opts['fwhm_bounds_factors'])*fwhm
        for group in np.unique(groups):
            indices=np.flatnonzero(groups==group);xy=initial[indices]
            x0=max(0,int(np.floor(xy[:,0].min()))-13);y0=max(0,int(np.floor(xy[:,1].min()))-13)
            x1=min(a.shape[1],int(np.ceil(xy[:,0].max()))+14);y1=min(a.shape[0],int(np.ceil(xy[:,1].max()))+14)
            sy=slice(y0,y1);sx=slice(x0,x1);data=a[sy,sx];bad=mask[sy,sx];err=error[sy,sx]
            yy,xx=np.indices(data.shape);blank=~bad.copy()
            for x,y in initial:blank&=np.hypot(xx+x0-x,yy+y0-y)>7
            if blank.sum()<30:
                rows.extend(dict(index=int(i),status='rejected',reason='insufficient_local_background') for i in indices);continue
            _,sky,noise=sigma_clipped_stats(data[blank]);passes=[]
            try:
                for shape in (opts['fit_shape'],opts['check_shape']):
                    model=CircularGaussianPRF(fwhm=fwhm);model.fwhm.fixed=False;model.fwhm.bounds=tuple(bounds);model.flux.bounds=(0,None)
                    photo=PSFPhotometry(model,int(shape),xy_bounds=opts['xy_bounds'],aperture_radius=3.,grouper=SourceGrouper(opts['group_separation_pix']),fitter=TRFLSQFitter(),fitter_maxiters=200)
                    passes.append(photo(data-float(sky),mask=bad,error=err,
                        init_params=Table(dict(x_init=xy[:,0]-x0,y_init=xy[:,1]-y0,group_id=np.ones(len(xy),int)))))
                for k,index in enumerate(indices):
                    r,s=passes[0][k],passes[1][k];x=float(r['x_fit'])+x0;y=float(r['y_fit'])+y0
                    width=float(r['fwhm_fit']);chi=float(r['reduced_chi2']);qfit=float(r['qfit'])
                    stability=float(np.hypot(r['x_fit']-s['x_fit'],r['y_fit']-s['y_fit']))
                    formal=float(np.hypot(r['x_err'],r['y_err']))
                    uncertainty=max(float(detection.get('point_error_floor',.08)),formal*np.sqrt(max(1.,chi)),stability)
                    reasons=[];flags=int(r['flags']);check_flags=int(s['flags'])
                    if not np.isfinite([x,y,width,chi,qfit,formal,uncertainty,r['flux_fit'],r['flux_err'],s['x_fit'],s['y_fit'],s['x_err'],s['y_err'],s['fwhm_fit'],s['flux_fit']]).all():reasons.append('nonfinite_or_unconstrained_fit')
                    if (flags|check_flags)&~1:reasons.append('native_psf_fit_flag')
                    if width<=bounds[0]*1.01 or width>=bounds[1]*.99:reasons.append('nonstellar_width_at_bound')
                    if qfit>opts['max_qfit']:reasons.append('nonstellar_profile_residual')
                    if stability>opts['max_window_shift_pix']:reasons.append('unstable_psf_centroid')
                    if uncertainty>opts['max_centroid_error_pix']:reasons.append('poor_psf_precision')
                    if float(r['flux_fit'])<=0 or float(r['flux_err'])<=0:reasons.append('nonpositive_flux_or_uncertainty')
                    snr=float(r['flux_fit'])/max(float(r['flux_err']),1e-12)
                    raw=Table(dict(x=[x],y=[y],extver=[version],flux=[float(r['flux_fit'])],snr=[snr],local_sky=[float(sky)]))
                    clean,audit=validate_centroids(image,raw,detection,method='PSFPhotometry',contexts={version:ctx},return_all=True)
                    if not len(clean):reasons.append('pixel_qc:'+str(audit['quality_flags'][0]))
                    entry=dict(index=int(index),source_id=str(catalog['hst_source_id'][index]),extver=version,
                        initial_x=float(initial[index,0]),initial_y=float(initial[index,1]),x=x,y=y,
                        fwhm_pix=width,width_pix=float(audit['width_pix'][0]),ellipticity=float(audit['ellipticity'][0]),
                        centroid_error_pix=uncertainty,formal_centroid_error_pix=formal,snr=snr,flux=float(r['flux_fit']),flux_error=float(r['flux_err']),
                        original_detection_method=str(catalog['detection_method'][index]) if 'detection_method' in catalog.colnames else 'custom',
                        original_segment_label=int(catalog['segment_label'][index]) if 'segment_label' in catalog.colnames and not np.ma.is_masked(catalog['segment_label'][index]) else None,
                        flux_rate=float(r['flux_fit'])/hd[0].header['EXPTIME'],local_sky=float(sky),local_noise=float(noise),
                        qfit=qfit,reduced_chi2=chi,window_stability_pix=stability,flags=flags,check_window_flags=check_flags,group_size=len(indices),
                        status='rejected' if reasons else 'credible_stellar_fit',reason=';'.join(reasons),
                        core_valid_fraction=float(audit['core_valid_fraction'][0]),detection_method='PSFPhotometry')
                    entry['model_warning']='Analytic PRF residuals exceed the ERR-only noise model; uncertainties inflated. Not calibrated PSF photometry.' if chi>1 else ''
                    rows.append(entry)
            except Exception as exc:
                rows.extend(dict(index=int(i),status='rejected',reason=f'psf_exception:{type(exc).__name__}: {exc}') for i in indices)
    good=[]
    for r in sorted((r for r in rows if r['status']=='credible_stellar_fit'),key=lambda r:r['centroid_error_pix']):
        same=next((q for q in good if np.hypot(q['x']-r['x'],q['y']-r['y'])<opts['duplicate_radius_pix']),None)
        if same:r.update(status='duplicate_fit',reason='Independent initial detections converged to the same stellar centroid',duplicate_of=same['source_id'])
        else:good.append(r)
    # Rejected native fits may have undefined covariance or morphology. Keep
    # that evidence as explicit JSON nulls, never as accepted measurements.
    for row in rows:
        missing=[key for key,value in row.items()
                 if isinstance(value,(float,np.floating)) and not np.isfinite(value)]
        if missing:
            row['nonfinite_fields']=missing
            for key in missing:row[key]=None
    return good,dict(status='unique_stellar_profile' if len(good)==1 else ('ambiguous_stellar_profiles' if good else 'no_reliable_stellar_profile'),
        raw_candidates=len(catalog),duplicates=sum(r['status']=='duplicate_fit' for r in rows),rejected=sum(r['status']=='rejected' for r in rows),fits=rows,
        model='Pixel-integrated CircularGaussianPRF centroid proxy; not calibrated PSF photometry',
        uncertainty_method='Maximum of pixel floor, formal radial error inflated by sqrt(reduced chi2), and fitting-window shift',options=opts)

def choose_profile(good,anchor,options):
    """HST-only profile/flux consistency; no coordinate zero-point requirement."""
    if len(good)==1:return good[0]
    if not anchor:return None
    opts=stellar_options(options);passing=[]
    for r in good:
        fluxerr=np.hypot(r['flux_error']/max(r['flux'],1e-12),anchor['flux_error']/max(anchor['flux'],1e-12))
        if 'concentration' in r and 'concentration' in anchor:
            if not np.allclose(r['aperture_radii_pix'],anchor['aperture_radii_pix']):continue
            delta=abs(r['concentration']-anchor['concentration'])
            allowed=max(opts['profile_concentration_fraction']*abs(anchor['concentration']),3*np.hypot(r['concentration_error'],anchor['concentration_error']))
            profile_ok=delta<=allowed
        else:profile_ok=abs(r['fwhm_pix']/anchor['fwhm_pix']-1)<=opts['profile_width_fraction']
        if abs(r['flux_rate']/anchor['flux_rate']-1)<=max(opts['profile_flux_fraction'],3*fluxerr) and profile_ok:passing.append(r)
    return passing[0] if len(passing)==1 else None

def install_measurement(catalog,result):
    idx=np.flatnonzero(np.asarray(catalog['hst_source_id']).astype(str)==result['hst_source_id'])
    if len(idx)!=1:raise ValueError('Profile assessment lacks its native HST detection identity')
    k=int(idx[0])
    for field in ('x','y','snr','width_pix','ellipticity','centroid_error_pix'):
        if field in catalog.colnames:catalog[field][k]=result[field]
    catalog['flux'][k]=result['psf_selected']['flux']
    if 'detection_method' in catalog.colnames:catalog['detection_method']=np.asarray(catalog['detection_method']).astype('U64');catalog['detection_method'][k]=result['psf_selected'].get('detection_method','PSFPhotometry')

def apply_selection(result,selected):
    if selected is None:return result
    result.update(status='valid_detection',reason=('Concentric-aperture compact profile and existing HST centroid QC' if selected.get('detection_method')=='aperture_profile' else 'Native-pixel PRF fit and shared centroid QC')+'; external position used for identification only',
        extver=selected['extver'],x=selected['x'],y=selected['y'],snr=selected['snr'],width_pix=selected['width_pix'],
        ellipticity=selected['ellipticity'],centroid_error_pix=selected['centroid_error_pix'],flux_rate=selected['flux_rate'],
        hst_source_id=selected['source_id'],psf_selected=selected)
    return result

def diagnostic(image,result,paths,label):
    if not paths or not result.get('psf_audit'):return
    from ...common.target_runner import save
    audit=result['psf_audit'];save(paths.file('intermediate','alignment',label+'_stellar_psf'),audit)
    # Individual fit records stay in the audit, not nested structures in ECSV.
    result['psf_audit_file']=str(paths.file('intermediate','alignment',label+'_stellar_psf'))

def profile_figure(inputs,records,paths,source_id):
    if not paths:return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import AsinhStretch,ImageNormalize,PercentileInterval
    present=[r for r in records if r.get('psf_audit')]
    if not present:return
    fig,axes=plt.subplots((len(present)+1)//2,2,figsize=(10,4*((len(present)+1)//2)),squeeze=False)
    for ax,r in zip(axes.flat,present):
        image=next(p for p in inputs if p.name==r['exposure']);v=int(r['extver']);x=r['predicted_x'];y=r['predicted_y'];radius=52
        with fits.open(image) as hd:
            a=hd['SCI',v].data;x0=max(0,int(x)-radius);y0=max(0,int(y)-radius)
            cut=a[y0:min(a.shape[0],int(y)+radius+1),x0:min(a.shape[1],int(x)+radius+1)]
            ax.imshow(cut,origin='lower',extent=(x0,x0+cut.shape[1],y0,y0+cut.shape[0]),cmap='gray',norm=ImageNormalize(cut,interval=PercentileInterval(99.5),stretch=AsinhStretch()))
        chosen=(r.get('psf_selected') or r.get('selected_hst') or {}).get('source_id')
        for q in r['psf_audit']['fits']:
            if q.get('x') is None or q.get('y') is None:continue
            color='lime' if q.get('source_id')==chosen else ('orange' if q['status']=='credible_stellar_fit' else 'red')
            ax.plot(q['x'],q['y'],'o',mfc='none',mec=color,ms=9);ax.annotate(str(q['index']+1),(q['x'],q['y']),color=color,fontsize=8)
        if r.get('selected_hst'):
            q=r['selected_hst'];ax.plot(q['x'],q['y'],'+',color='lime',ms=12,mew=2)
        ax.set_title(image.stem+' SCI'+str(v)+'\n'+r['status']);ax.set_xlabel('Native detector X');ax.set_ylabel('Native detector Y')
    for ax in list(axes.flat)[len(present):]:ax.axis('off')
    fig.suptitle(source_id+' — green selected; red rejected; orange unresolved\nHST profile assessment; catalog coordinates are identification priors only',fontsize=11)
    fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment','stellar_association_'+source_id,'png'),dpi=100);plt.close(fig)
