"""Conservative off-source diagnostics; AstroDrizzle owns sky subtraction.

Masks must be supplied and inspected by the caller. No object classification or
amplifier geometry is inferred from the science signal. Label maps use zero for
excluded pixels and positive integers for independently selected sky regions or
amplifiers. Statistics remain in each SCI extension's declared input units.
"""
import numpy as np
from astropy.io import fits

def assess_regions(exposure,mask):
    report=[]
    with fits.open(exposure) as science,fits.open(mask) as masks:
        for h in science:
            if h.name!='SCI':continue
            v=h.header.get('EXTVER',1);labels=masks['SCI',v].data
            if labels.shape!=h.data.shape:raise ValueError('Sky/amplifier mask shape differs from SCI')
            if not np.isfinite(labels).all() or np.any(labels<0) or not np.equal(labels,np.floor(labels)).all():
                raise ValueError('Masks must contain nonnegative integer region labels')
            dq=science['DQ',v].data
            for label in np.unique(labels[labels>0]):
                good=(labels==label)&(dq==0)&np.isfinite(h.data)
                values=h.data[good].astype(float)
                if len(values)<30:raise ValueError('Insufficient independent sky pixels for region assessment')
                report.append(dict(extver=v,region=int(label),npixel=len(values),
                    mean=float(values.mean()),rms=float(values.std(ddof=1)),
                    standard_error=float(values.std(ddof=1)/np.sqrt(len(values))),
                    unit=h.header.get('BUNIT','UNKNOWN'),
                    caveat='Conditional on supplied source-free masks; correlated noise not modeled'))
    return report

def estimate(data,err,dq,wcs,target,options):
    """Mask sources at native sampling; protect target independent of detection.

    Uncertainty includes spatial variation between sky tiles; it is deliberately
    larger than a naive pixel standard error. Auto never fits spatial surfaces.
    """
    from scipy.ndimage import gaussian_filter,binary_dilation
    from astropy.stats import sigma_clipped_stats
    from .dq import ignored_bits
    good=np.isfinite(data)&np.isfinite(err)&(err>0)&((dq&np.uint16(65535^ignored_bits(options,'background')))==0)
    if good.sum()<30:raise ValueError('No usable SCI/ERR/DQ sky pixels')
    _,med,std=sigma_clipped_stats(data,mask=~good,sigma=3,maxiters=5)
    noise=max(float(std),float(np.median(err[good])),np.finfo(float).eps)
    smoothed=gaussian_filter(np.where(good,data-med,0).astype(float),1)
    sources=(smoothed>options['threshold']*noise/3.0)|(data>med+5*noise)
    mask=(~good)|binary_dilation(sources,iterations=max(1,int(options['dilate'])))
    x,y=wcs.all_world2pix([target],0)[0]
    scale=np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))*3600
    yy,xx=np.ogrid[:data.shape[0],:data.shape[1]]
    mask|=(xx-x)**2+(yy-y)**2<(options['target_radius']/scale)**2
    blank=~mask;values=data[blank]
    n=int(blank.sum());fraction=n/data.size
    tiles=[];box=int(options['box_size'])
    for iy in range(0,data.shape[0],box):
        for ix in range(0,data.shape[1],box):
            tile=data[iy:iy+box,ix:ix+box];ok=blank[iy:iy+box,ix:ix+box]
            if ok.sum()>=max(30,0.25*tile.size):tiles.append(float(np.median(tile[ok])))
    reliable=n>=options['min_pixels'] and fraction>=options['min_blank'] and len(tiles)>=4
    value=float(np.median(values)) if n else 0.0
    error=max(noise/np.sqrt(max(n,1)),float(np.std(tiles))/np.sqrt(max(len(tiles),1))) if tiles else noise
    report=dict(measured=value,uncertainty=error,blank_fraction=fraction,npixels=n,
                tile_medians=tiles,reliable=reliable,mode=options['mode'],warnings=[])
    model=np.zeros(data.shape,dtype='f4');mode=options['mode']
    if mode=='user':raise ValueError('Supply user sky through prepare_background with per-extension values')
    if mode!='none' and reliable:
        if mode=='auto' and abs(value)<options['significance']*error:
            report['warnings'].append('Sky not significant; zero subtraction fallback')
        elif mode=='2d':
            from photutils.background import Background2D,MedianBackground
            model=Background2D(data,box_size=box,mask=mask,filter_size=3,
                bkg_estimator=MedianBackground(),exclude_percentile=80).background.astype('f4')
        else:model.fill(value)
    elif mode!='none':report['warnings'].append('Insufficient blank sky; zero subtraction fallback')
    report['applied']=float(np.median(model[blank])) if n else 0.0
    report['structure_peak_to_peak']=float(np.ptp(tiles)) if tiles else None
    return model,mask,report

def prepare_background(source,destination,directory,options,*,paths=None):
    """Subtract only on a working copy; record original sky for CR noise model."""
    import shutil
    from astropy.wcs import WCS
    from .qa import image_plot
    paths.copy(source,destination) if paths else shutil.copy2(source,destination)
    models=[fits.PrimaryHDU()];masks=[fits.PrimaryHDU()];reports=[]
    with fits.open(destination,mode='update') as hd:
        from .fits import image_data
        for h in hd:
            if h.name in ('ERR','DQ') and h.data is None:h.data=image_data(h)
        from .dq import native_compatible_dq
        native_compatible_dq(hd)
        if hd[0].header.get('HPSKYSUB',False):raise ValueError('Sky already subtracted by hstpipe')
        target=[hd[0].header['RA_TARG'],hd[0].header['DEC_TARG']]
        amp_report=None
        if options.get('amplifier_assess') or options.get('amplifier_correction'):
            from .amplifiers import equalize
            amp_models,amp_report=equalize(hd,options)
            if amp_models:
                ah=[fits.PrimaryHDU()]
                for av,ac in amp_models.items():
                    a=fits.ImageHDU(ac,name='SCI');a.header['EXTVER']=av;a.header['BUNIT']='ELECTRONS';ah.append(a)
                fits.HDUList(ah).writeto(paths.exposure_file('intermediate','background',source.stem.replace('_flc',''),'amp') if paths else directory/(source.stem.replace('_flc','')+'_amp.fits'),overwrite=True,checksum=True)
        for h in hd:
            if h.name!='SCI':continue
            v=h.header.get('EXTVER',1)
            opts=dict(options);mode=opts['mode'];filt=hd[0].header.get('FILTER','')
            pivot=h.header.get('PHOTPLAM',hd[0].header.get('PHOTPLAM',float('nan')))
            regime='narrowband' if filt.endswith('N') or filt.startswith('FQ') else ('low_background_UV' if pivot<3500 else 'broad_or_medium_optical')
            if options.get('source_fills_field'):regime='extended_source_fills_field'
            opts['mode']='none' if mode=='user' or (mode=='auto' and regime in ('narrowband','extended_source_fills_field')) else mode
            model,mask,report=estimate(h.data,hd['ERR',v].data,hd['DQ',v].data,WCS(h.header,hd),target,opts)
            report.update(requested_mode=mode,regime=regime,filter=filt)
            if mode=='auto' and opts['mode']=='none':report['warnings'].append('Automatic sky subtraction disabled for narrowband or field-filling emission; measured sky retained as information')
            if regime=='low_background_UV':report['warnings'].append('Low-background UV: source-masked significance/blank-sky checks required; residual CTE needs review')
            if mode=='user':
                user=options['user'].get(source.stem.replace('_flc',''),{})
                entry=user.get(str(v))
                if entry is None:raise ValueError(f'Missing user sky {source.name} SCI,{v}')
                value,error=entry
                if not np.isfinite([value,error]).all() or error<0:raise ValueError('Invalid user sky/uncertainty')
                model.fill(value);report.update(mode=mode,applied=float(value),uncertainty=float(error))
            h.data-=model
            h.header['HPSKY']=report['applied'];h.header['HPSKYERR']=report['uncertainty']
            h.header['HPSKYBG']=(report['measured'],'Measured pre-subtraction sky electrons; may be unreliable')
            h.header['HPSKYOK']=bool(report['reliable'])
            h.header['MDRIZSKY']=report['applied']
            report['sky_handling']='skysub=False, skyuser=HPSKY (already removed); MDRIZSKY records applied sky; HPSKYBG retains measured sky. Unsubtracted sky remains in SCI/reference for CR noise.'
            h.header.add_history('hstpipe: source-masked sky subtracted once; HPSKY is electrons')
            header=h.header.copy();models.append(fits.ImageHDU(model,header=header,name='SCI'))
            mask_header=header.copy();mask_header['BUNIT']='UNITLESS';mask_header['MASKONE']='EXCLUDED'
            masks.append(fits.ImageHDU(mask.astype('u1'),header=mask_header,name='SCI'))
            report['extver']=v
            if amp_report:report['amplifiers']=amp_report
            reports.append(report)
            image_plot([h.data,model,mask],['Sky-subtracted SCI','Applied sky (electrons)','Excluded sky pixels'],
                       paths.exposure_file('diagnostics','background',source.stem.replace('_flc',''),f'sci{v}','png',group=True) if paths else directory.parent.parent/'diagnostics'/'background'/f'{source.stem}_sci{v}.png')
        hd[0].header['HPSKYSUB']=True
        hd[0].header['HPWORK']=True
    root=source.stem.replace('_flc','')
    fits.HDUList(models).writeto(paths.exposure_file('intermediate','background',root,'sky') if paths else directory/f'{root}_sky.fits',overwrite=True,checksum=True)
    fits.HDUList(masks).writeto(paths.exposure_file('intermediate','background',root,'mask') if paths else directory/f'{root}_mask.fits',overwrite=True,checksum=True)
    return reports


def restore_background(source,destination,paths):
    """Restore a saved checkpoint by replaying its stored arrays, without estimation."""
    import shutil
    from .fits import image_data
    root=source.stem.replace('_flc','')
    sky=paths.exposure_file('intermediate','background',root,'sky')
    amp=paths.exposure_file('intermediate','background',root,'amp')
    paths.copy(source,destination)
    with fits.open(destination,mode='update') as hd,fits.open(sky) as models:
        from .dq import native_compatible_dq
        native_compatible_dq(hd)
        amps=fits.open(amp) if amp.exists() else None
        try:
            for h in hd:
                if h.name in ('ERR','DQ') and h.data is None:h.data=image_data(h)
                if h.name!='SCI':continue
                v=h.header.get('EXTVER',1)
                if amps is not None and ('SCI',v) in amps:h.data-=amps['SCI',v].data
                h.data-=models['SCI',v].data
                h.header['HPSKY']=models['SCI',v].header['HPSKY']
                h.header['HPSKYERR']=models['SCI',v].header['HPSKYERR']
                for key in ('HPSKYBG','HPSKYOK','MDRIZSKY'):
                    if key in models['SCI',v].header:h.header[key]=models['SCI',v].header[key]
                h.header.add_history('hstpipe: source-masked sky subtracted once; HPSKY is electrons')
            hd[0].header['HPSKYSUB']=True
            hd[0].header['HPWORK']=True
        finally:
            if amps is not None:amps.close()
