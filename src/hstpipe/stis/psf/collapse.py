"""Finite-bin projection and wavelength collapse of explicitly selected PSFs."""
from pathlib import Path
import json,hashlib
import numpy as np
from astropy.io import fits
from astropy.table import Table


def pixel_overlap_matrix(source_centers,target_centers,shift=0.):
    s=np.asarray(source_centers,float)+shift;t=np.asarray(target_centers,float)
    ds=float(np.median(np.diff(s)));dt=float(np.median(np.diff(t)))
    if ds<=0 or dt<=0 or not np.allclose(np.diff(s),ds) or not np.allclose(np.diff(t),dt):raise ValueError('Increasing uniform grids required')
    if t[0]-dt/2<s[0]-ds/2-1e-9 or t[-1]+dt/2>s[-1]+ds/2+1e-9:raise ValueError('Target rows exceed measured PSF support')
    return np.maximum(0,np.minimum(t[:,None]+dt/2,s[None,:]+ds/2)-np.maximum(t[:,None]-dt/2,s[None,:]-ds/2))/ds


def collapse_library_profile(library,request,output):
    req=json.loads(Path(request).read_text()) if not isinstance(request,dict) else request
    wave=np.asarray(req['wavelength_angstrom'],float);y=np.asarray(req['position_arcsec'],float);centers=np.asarray(req['trace_center_arcsec'],float);w=np.asarray(req.get('weights',np.ones(len(wave))),float)
    if len(wave)!=len(centers) or len(w)!=len(wave) or np.any(w<0) or w.sum()<=0:raise ValueError('Invalid collapse wavelengths/weights')
    w=w/w.sum();components=[];metadata=[]
    with fits.open(library,checksum=True) as h:
        entries=[h[name] for name in req['templates']];entries.sort(key=lambda x:(x.header['WAVEMIN']+x.header['WAVEMAX'])/2)
        if len(entries)!=2:raise ValueError('Two bracketing measured wavelength bands are required')
        a,b=entries;keys=['APERTURE','EPOCH','POSGROUP','SAMPLING','NORMRAD']
        if any(a.header[k]!=b.header[k] for k in keys):raise ValueError('Incompatible template groups')
        if a.header['APERTURE']!=req['science_aperture'] and not req.get('allow_aperture_proxy',False):raise ValueError('Unvalidated aperture transfer requires explicit proxy flag')
        if a.header['SAMPLING']!=req['sampling']:raise ValueError('Sampling mismatch')
        left=(a.header['WAVEMIN']+a.header['WAVEMAX'])/2;right=(b.header['WAVEMIN']+b.header['WAVEMAX'])/2
        f=(wave-left)/(right-left)
        if np.any(f<0) or np.any(f>1):raise ValueError('Wavelength extrapolation prohibited')
        for entry,alpha in [(a,1-f),(b,f)]:
            source=np.array(entry.data['offset_arcsec'],float);p=np.array(entry.data['profile'],float);cov=np.array(h['PCOV'+entry.name[3:]].data,float)
            operator=sum(float(ww*aa)*pixel_overlap_matrix(source,y,shift) for ww,aa,shift in zip(w,alpha,centers))
            components.append((operator@p,operator@cov@operator.T))
            metadata.append(dict(extension=entry.name,aperture=entry.header['APERTURE'],wave_min=entry.header['WAVEMIN'],wave_max=entry.header['WAVEMAX'],epoch=entry.header['EPOCH'],trace_y_min=entry.header['YMIN'],trace_y_max=entry.header['YMAX'],normalization_radius=entry.header['NORMRAD'],grade=entry.header['GRADE'],scatter_status=entry.header['SCATSTAT']))
    p=sum(v for v,cv in components);peak=int(np.argmax(p));peakval=p[peak]
    if peakval<=0:raise ValueError('No positive PSF peak')
    selector=np.eye(len(p))[peak];jac=np.eye(len(p))/peakval-p[:,None]*selector[None,:]/peakval**2
    # Cross-band covariance is unknown: sum marginal standard deviations as a
    # fully correlated bound rather than silently declaring bands independent.
    error=sum(np.sqrt(np.maximum(np.diag(jac@cov@jac.T),0)) for v,cov in components)
    table=Table(dict(position_arcsec=y,profile_peak_normalized=p/peakval,template_error_bound=error))
    table.meta.update(science_aperture=req['science_aperture'],template_aperture=metadata[0]['aperture'],sampling=req['sampling'],center_at_reference_arcsec=req['center_at_reference_arcsec'],proxy=req['science_aperture']!=metadata[0]['aperture'],weighting=req.get('weighting','explicit supplied wavelength weights'))
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True);table.write(out,format='ascii.ecsv',overwrite=True)
    meta=dict(library=str(library),library_sha256=hashlib.sha256(Path(library).read_bytes()).hexdigest(),templates=metadata,wavelength_min=float(wave.min()),wavelength_max=float(wave.max()),n_columns=len(wave),request=req,normalization='Peak normalized on the actual science spatial grid; no width fit',uncertainty='Conditional predictive covariance; fully correlated cross-band bound, fixed trace and peak index; single-star inter-star variance unknown',limitations='Regular aperture may be a proxy for D1; epochs/positions differ. Uniform spectral collapse is a visualization, not an extension significance test. Pixel overlap can alter undersampled cores.')
    out.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n');return table
