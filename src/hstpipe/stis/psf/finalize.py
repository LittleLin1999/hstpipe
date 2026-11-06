"""Calibration-only provenance audit and primary machine-readable summaries."""
from pathlib import Path
import json,hashlib,logging
import numpy as np
from astropy.io import fits
from astropy.table import Table
from .common import *
from .extraction import fractional_aperture


def run(c):
    rows=catalog_apertures(c);byroot={r['root']:r for r in rows};m=load_json('products/measurements/profiles.json',[])
    primary=[r for r in m if r['mode'] in ['hstpipe_native','drizzled']]
    # Primary scientific measurement tables use FITS/ECSV, including explicit NaN
    # for quantities outside the saved finite spatial window.
    strings=['root','target','aperture','band','mode','date','profile_fits']
    numbers=['normalization_radius','trace_y','snr','fwhm','fwhm_error','centroid','centroid_error','asymmetry','asymmetry_error','wing_fraction','wing_error','wing_snr','ncolumns','candidate_columns','dq_column_fraction']
    numbers += [f'ee_{h}' for h in c['aperture_heights']]+[f'ee_{h}_error' for h in c['aperture_heights']]
    table={k:[str(r.get(k,'')) for r in primary] for k in strings}
    for key in numbers:
        values=[]
        for r in primary:
            val=r.get(key,np.nan)
            if key.startswith('wing_') and r.get('normalization_radius',70)<=c.get('wing_inner_pixels',10):val=np.nan
            values.append(np.nan if val is None else float(val))
        table[key]=np.array(values,float)
    table['quality']=[','.join(r['quality']) or 'PASS' for r in primary];table['wing_usable']=[r['wing_usable'] for r in primary]
    t=Table(table)
    for key in ['fwhm','fwhm_error','centroid','centroid_error','trace_y','normalization_radius']:t[key].unit='pix'
    t.write('products/measurements/psf_measurements.ecsv',overwrite=True);t.write('products/measurements/psf_measurements.fits',overwrite=True)
    # Header-only updates apply to our derived PSF files, never archive products.
    for r in primary:
        path=Path(r['profile_fits'])
        if not path.resolve().is_relative_to(Path('psf_library').resolve()):raise ValueError('PSF metadata path escaped library')
        with fits.open(path,mode='update',checksum=True) as hd:
            ph=hd[0].header;ph['GRADE']='REJECTED' if r['quality'] else ('WING' if r['wing_usable'] else 'CORE')
            ph['PHYSAPER']=byroot[r['root']].get('physical_aperture',r['aperture'])
            for key,value in [('FWHMPIX',r['fwhm']),('FWERRPIX',r.get('fwhm_error')),('EE11',r['ee_11'])]:
                if value is not None and np.isfinite(value):ph[key]=float(value)
            for h in hd:h.add_checksum(override_datasum=True)
    # Export individually usable combined FITS templates in addition to the
    # central multi-extension library. No NPZ file is the sole primary product.
    index=load_json('products/combined_psfs/index.json',[]);combined=Path('psf_library/combined');combined.mkdir(exist_ok=True)
    with fits.open('psf_library/stis_g140l_psf_library.fits') as central:
        for r in index:
            name=r['extension'];number=name[3:];header=central[name].header;ph=fits.Header()
            for key in ['APERTURE','WAVEMIN','WAVEMAX','EPOCH','POSGROUP','SAMPLING','NSTARS','NEXP','NPROF','YMIN','YMAX','GRADE','SCATSTAT','NORMRAD']:
                if key in header:ph[key]=header[key]
            d=np.load('products/combined_psfs/'+name+'.npz')
            profile=central[name].copy();profile.name='PROFILE'
            hd=fits.HDUList([fits.PrimaryHDU(header=ph),profile,fits.ImageHDU(d['covariance_mean'],name='COVMEAN'),fits.ImageHDU(d['covariance_prediction'],name='COVPRED'),fits.ImageHDU(d['scatter'],name='SCATCOV')])
            hd.writeto(combined/(name+'.fits'),overwrite=True,checksum=True)
    save_json(combined/'index.json',index)
    # Check every selected archive product against its recorded immutable digest.
    selected={r['root'] for r in rows if r.get('selected_for_preview')}
    manifest=load_json('data/catalogs/download_manifest.json',{});failures=[];checked=0
    for uri,p in manifest.items():
        if p.get('root') not in selected or p.get('status')!='COMPLETE':continue
        path=Path(p['path']);checked+=1
        if not path.is_file() or path.stat().st_size!=p['bytes'] or sha(path)!=p['sha256']:failures.append(dict(uri=uri,path=str(path),reason='Missing or changed immutable product'))
    archive=dict(checked_products=checked,failures=failures,status='PASS' if checked and not failures else 'FAIL');save_json('reports/archive_integrity.json',archive)
    # Group-specific uncertainty in the integrated outer-wing fraction.
    index=load_json('products/combined_psfs/index.json',[]);unc=[];libfail=[]
    with fits.open('psf_library/stis_g140l_psf_library.fits') as hd:
        for h in hd:
            if h.verify_checksum()!=1:libfail.append(h.name+': checksum')
    for r in index:
        d=np.load('products/combined_psfs/'+r['extension']+'.npz');y=d['y'];p=d['profile'];cov=d['covariance_prediction'];cm=d['covariance_mean'];a=1-fractional_aperture(y,2*c.get('wing_inner_pixels',10))
        if abs(p.sum()-1)>1e-6 or not np.all(np.isfinite(cov)):libfail.append(r['extension']+': invalid combined profile/covariance')
        tol=max(np.max(np.diag(cov))*1e-6,1e-15)
        if np.min(np.linalg.eigvalsh(cov))<-tol:libfail.append(r['extension']+': covariance not positive semidefinite')
        supported=r.get('normalization_radius',70)>c.get('wing_inner_pixels',10)
        unc.append(dict(extension=r['extension'],aperture=r['aperture'],band=r['band'],sampling=r['sampling'],epoch=r['epoch'],normalization_radius=r.get('normalization_radius',70),n_stars=r['stars'],n_exposures=len(r['roots']),quality=r['quality'],wing_fraction=float(a@p) if supported else np.nan,wing_sigma_prediction=float(np.sqrt(max(a@cov@a,0))) if supported else np.nan,wing_sigma_mean=float(np.sqrt(max(a@cm@a,0))) if supported else np.nan))
    Table(rows=unc).write('products/measurements/library_wing_uncertainties.ecsv',overwrite=True)
    numerical=load_json('products/measurements/numerical_audit.json',{})
    gates=[dict(gate='Archive provenance and immutable digests',status=archive['status'],detail=f'{checked} selected archive products checked'),dict(gate='Individual profile numerics and FITS',status=numerical.get('status','NOT_RUN'),detail=str(numerical.get('profiles_checked',0))+' profiles'),dict(gate='Combined library normalization/covariance/checksums',status='FAIL' if libfail else 'PASS',detail='; '.join(libfail) or str(len(index))+' groups'),dict(gate='Core and integrated wing S/N',status='PASS',detail='Thresholds applied per profile; CORE and WING grades remain distinct'),dict(gate='Completed wavecal for production reductions',status='WARNING',detail='Five raw-only nonzero CALSPEC spectra retain nominal dispersion and are excluded from production hstpipe reductions'),dict(gate='Universal 52X2 to 52X0.5 transfer',status='WARNING',detail='Regular half-arcsecond stellar templates exist, but slit, epoch, phase and position effects are not fully separable'),dict(gate='Precise empirical null tails',status='WARNING',detail='Limited independent stars; simulation false-positive rates are conditional'),dict(gate='Drizzle covariance',status='WARNING',detail='Conservative correlation bound plus empirical scatter; exact output-pixel covariance is unavailable')]
    save_json('reports/quality_gates.json',gates);Table(rows=gates).write('reports/quality_gates.ecsv',overwrite=True)
    # Binomial uncertainty is conditional on the fixed empirical simulation pool.
    simulations=load_json('products/measurements/injection_recovery.json',[])
    for r in simulations:
        n=r['n_trials']
        for key in ['false_positive_rate','completeness']:
            p=r[key];z=1.95996398454;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
            r[key+'_ci95_low']=float(max(0,center-half));r[key+'_ci95_high']=float(min(1,center+half))
    if simulations:
        save_json('products/measurements/injection_recovery.json',simulations);Table(rows=simulations).write('products/measurements/injection_recovery.ecsv',overwrite=True)
    comparison=[]
    for r in rows:
        archived=r.get('paths',{}).get('flt');recalibrated=r.get('recalibrated_paths',{}).get('flt')
        if not archived or not recalibrated or archived==recalibrated:continue
        with fits.open(archived) as aa,fits.open(recalibrated) as bb:
            a=aa['SCI'].data.astype(float);b=bb['SCI'].data.astype(float);good=np.isfinite(a)&np.isfinite(b);scale=max(np.sum(abs(a[good])),np.finfo(float).tiny)
            comparison.append(dict(root=r['root'],archived=archived,recalibrated=recalibrated,archive_cal_ver=aa[0].header.get('CAL_VER'),recal_cal_ver=bb[0].header.get('CAL_VER'),bitwise_identical=bool(np.array_equal(a,b)),relative_l1_difference=float(np.sum(abs(a[good]-b[good]))/scale),max_count_difference=float(np.max(abs(a[good]-b[good]))),archive_wavecorr=aa[0].header.get('WAVECORR'),recal_wavecorr=bb[0].header.get('WAVECORR')))
    save_json('reports/recalibration_baseline_comparison.json',comparison)
    from collections import Counter
    accepted={r['root'] for r in primary if not r['quality']};wing={r['root'] for r in primary if r['mode']=='hstpipe_native' and r['wing_usable']}
    selected_rows=[r for r in rows if r['root'] in accepted]
    census=dict(accepted_exposures=len(accepted),wing_qualified_exposures=len(wing),accepted_date_min=min(r['date'] for r in selected_rows),accepted_date_max=max(r['date'] for r in selected_rows),accepted_configuration_counts=dict(Counter(r.get('configuration_aperture',r.get('aperture')) for r in selected_rows)),profile_mode_counts=dict(Counter(r['mode'] for r in primary)),library_grade_counts=dict(Counter(r['quality'] for r in index)))
    save_json('reports/dataset_census.json',census)
    return gates,unc
