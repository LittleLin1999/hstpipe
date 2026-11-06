"""PSF-only execution, calibration, and validation reports."""
from pathlib import Path
from collections import Counter,defaultdict
import json,sys,platform
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from .common import *
from .psf import eligible


def distribution(rows,key):
    v=np.array([r[key] for r in rows if r.get(key) is not None],float);v=v[np.isfinite(v)]
    if not len(v):return 'not measured'
    return f'n={len(v)}; median={np.median(v):.5g}; 16–84 percentile=[{np.quantile(v,.16):.5g}, {np.quantile(v,.84):.5g}]'


def numerical_audit(m):
    errors=[];checked=0
    for r in m:
        if r['mode'] not in ['hstpipe_native','drizzled']:continue
        d=np.load(r['profile_file']);cv=d['normalized_covariance'];n=d['normalized'];checked+=1
        if not np.all(np.isfinite(cv)) or not np.all(np.isfinite(n)):errors.append(r['profile_file']+': nonfinite')
        if abs(n.sum()-1)>1e-6:errors.append(r['profile_file']+': normalization')
        tol=max(np.max(np.diag(cv))*1e-6,1e-15)
        if np.max(abs(cv-cv.T))>tol or np.min(np.linalg.eigvalsh(cv)) < -tol:errors.append(r['profile_file']+': covariance')
        with fits.open(r['profile_fits']) as h:
            for hd in h:
                if hd.verify_checksum()!=1:errors.append(r['profile_fits']+': checksum')
    out=dict(profiles_checked=checked,failures=errors,status='PASS' if checked and not errors else 'FAIL');save_json('products/measurements/numerical_audit.json',out);return out


def primary_plots(c,rows,m):
    good=eligible(m,rows,modes=('hstpipe_native','drizzled'))
    for band in sorted(set(r['band'] for r in good)):
        fig,ax=plt.subplots(1,3,figsize=(15,4))
        for r in [x for x in good if x['band']==band]:
            d=np.load(r['profile_file']);y=d['y'];p=d['normalized'];dy=np.median(np.diff(y));a11=__import__('hstpipe.stis.psf.extraction',fromlist=['fractional_aperture']).fractional_aperture(y,11);den=float(a11@p);jac=np.eye(len(p))/den-p[:,None]*a11[None,:]/den**2;ratio_cov=jac@d['normalized_covariance']@jac.T;p=p/den;label=f'{r["root"]} {r["aperture"]} {r["mode"]} R={r.get("normalization_radius",70)}'
            ax[0].plot(y,p/dy,lw=.6,label=label);ax[1].plot(y,p/dy,lw=.6)
            heights=[3,5,7,9,11,15,21,41,81,140]
            from .extraction import fractional_aperture
            ee=[];ee_error=[]
            for height in heights:
                if height>2*r.get('normalization_radius',70):ee.append(np.nan);ee_error.append(np.nan);continue
                aa=fractional_aperture(y,height);ee.append(float(aa@p));ee_error.append(float(np.sqrt(max(aa@ratio_cov@aa,0))))
            ax[2].errorbar(np.array(heights)/2,ee,yerr=ee_error,lw=.5,markersize=1)
        ax[0].set(xlim=(-12,12),xlabel='Along-slit offset [native pixel]',ylabel='Profile / flux in 11-pixel aperture [pixel⁻¹]',title=band+' Å: hstpipe products')
        ax[1].set_yscale('symlog',linthresh=1e-5)
        ax[1].set(xlim=(-70,70),ylim=(-.001,.5),xlabel='Along-slit offset [native pixel]',ylabel='Profile / flux in 11-pixel aperture [pixel⁻¹]',title='Extended wings')
        ax[2].set(xlabel='Aperture half-height [native pixel]',ylabel='EE(height) / EE(11 pixels)',title='Common default-height normalization')
        ax[0].legend(fontsize=3,loc='upper left');fig.tight_layout();fig.savefig(f'figures/psf_{band}.png',dpi=150);plt.close(fig)
    for xkey,xlabel,name in [('trace_y','Detector trace position [native pixel]','detector_position'),('center_wave','Observed wavelength [Å]','wavelength')]:
        fig,ax=plt.subplots(figsize=(10,5))
        for ap in sorted(set(r['aperture'] for r in good)):
            g=[r for r in good if r['aperture']==ap and r['mode']=='hstpipe_native' and r['hi']-r['lo']<100]
            x=[(r['lo']+r['hi'])/2 if xkey=='center_wave' else r[xkey] for r in g]
            ax.errorbar(x,[r['fwhm']*.0246 for r in g],yerr=[(r.get('fwhm_error') or 0)*.0246 for r in g],fmt='.',ms=3,lw=.4,label=ap)
        ax.set(xlabel=xlabel,ylabel='Direct FWHM [arcsec]',title='Calibrator dependence; bands and epochs may be confounded');ax.legend();fig.tight_layout();fig.savefig('figures/'+name+'.png',dpi=140);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(12,4));all_dates=[r['mjd'] for r in rows];selected=[r['mjd'] for r in rows if 'flt' in r.get('paths',{})];ax[0].hist(all_dates,bins=35,label='All archive records',alpha=.5);ax[0].hist(selected,bins=35,label='Selected FLT sets',alpha=.7);ax[0].set(xlabel='Observation date [MJD]',ylabel='Records');ax[0].legend();aps=Counter(r.get('configuration_aperture',r.get('aperture','unknown')) for r in rows if 'flt' in r.get('paths',{}));ax[1].bar(list(aps),list(aps.values()));ax[1].set(xlabel='Exact FITS aperture',ylabel='Selected FLT exposure sets');fig.tight_layout();fig.savefig('figures/selection.png',dpi=140);plt.close(fig)


def pctab(c,rows,m):
    output=[]
    for r in m:
        if r['mode']!='hstpipe_native' or r['quality'] or r['aperture']!='52X2':continue
        p=Path('data/hstpipe_references/y2r16004o_pct.fits')
        if not p.exists():continue
        with fits.open(p) as h:
            t=h[1].data;ss=t[(t['CENWAVE']==1425)&(np.char.strip(t['APERTURE'])=='52X2')]
            for a in ss:
                height=int(a['EXTRHEIGHT']);key=f'ee_{height}'
                if key not in r:continue
                n=int(a['NELEM']);mid=(r['lo']+r['hi'])/2;wave=a['WAVELENGTH'][:n]
                if mid<min(wave) or mid>max(wave):continue
                output.append(dict(root=r['root'],band=r['band'],height=height,empirical_ratio_to_11=r[key]/r['ee_11'],pctab_ratio=float(np.interp(mid,wave,a['THROUGHPUT'][:n])),reference=str(p),sha256=sha(p)))
    save_json('products/measurements/pctab_comparison.json',output)
    if output:Table(rows=output).write('products/measurements/pctab_comparison.ecsv',overwrite=True)


def report(c):
    rows=catalog_apertures(c);m=load_json('products/measurements/profiles.json',[]);manifest=load_json('data/catalogs/download_manifest.json',{});reductions=load_json('data/catalogs/hstpipe_reductions.json',[]);library=load_json('products/combined_psfs/index.json',[])
    good=eligible(m,rows);accepted=set(r['root'] for r in good);audit=numerical_audit(m);primary_plots(c,rows,m);pctab(c,rows,m)
    from .finalize import run as finalize
    gates,wing_uncertainties=finalize(c)
    primary=[r for r in m if r['mode']=='hstpipe_native'];downloaded=[r for r in rows if 'flt' in r.get('paths',{})];screened=[r for r in rows if r.get('preview_done')];programs=[]
    for pid in [18169,17644,17649,18229,16966,18113,17804,17223]:
        rr=[r for r in rows if r['program']==pid];roots={r['root'] for r in rr};mm=[r for r in primary if r['root'] in roots];use=[r for r in mm if r['root'] in accepted and not r['quality']];ss=[r['snr'] for r in use];screen_snr=[r['snr'] for r in m if r['root'] in roots and r['mode']=='aligned' and r['background']==0]
        programs.append(dict(program=pid,matching_records=len(rr),previewed=sum(bool(r.get('preview_done')) for r in rr),flt_sets=sum('flt' in r.get('paths',{}) for r in rr),accepted_native=len(roots&accepted),wing_qualified=len({r['root'] for r in use if r['wing_usable']}),screen_rejected=sum(r['status'].startswith('REJECT') for r in rr),qa_rejected=sum(r.get('psf_status')=='REJECT_QA' for r in rr),downloaded_products=sum(p.get('status')=='COMPLETE' and p.get('root') in roots for p in manifest.values()),hstpipe_complete=sum(red['status']=='COMPLETE' and red.get('program')==pid for red in reductions),preview_snr_max=max((r.get('preview_snr',0) for r in rr),default=0),screen_profile_snr_min=min(screen_snr) if screen_snr else None,screen_profile_snr_max=max(screen_snr) if screen_snr else None,snr_min=min(ss) if ss else None,snr_max=max(ss) if ss else None))
    save_json('reports/program_summary.json',programs);pd.DataFrame(programs).to_csv('reports/program_summary.csv',index=False)
    keys=['root','obsid','program','target','date','aperture','physical_aperture','propaper','configuration_aperture','detector','opt_elem','cenwave','obsmode','exptime','preview_snr','measured_continuum_snr','status','reason','psf_status','processing_status','trace_y','postarg1','postarg2','shifta1','shifta2','cal_ver','sdqflags','wing_usable_bands']
    final_catalog=[]
    unsupported={root for red in reductions if red['status']=='UNSUPPORTED_WAVECORR' for root in red['roots']}
    for r in rows:
        item={k:r.get(k) for k in keys};item['accepted_hstpipe_native']=r['root'] in accepted
        if r['root'] in accepted:decision='ACCEPTED';reason='At least one native hstpipe profile passes configured QA'
        elif r['root'] in unsupported:decision='NOMINAL_ONLY';reason='WAVECORR is not COMPLETE; current hstpipe production gate is not satisfied'
        elif r['status']=='REJECT_ZERO_EXPTIME':decision='REJECTED';reason='Actual raw SCI EXPTIME and counts are zero'
        elif r.get('psf_status')=='REJECT_QA':
            decision='REJECTED';flags=sorted({q for a in m if a['root']==r['root'] and a['mode']=='aligned' and a['background']==0 for q in a['quality']});reason=', '.join(flags) or 'No supported profile'
        elif not r.get('selected_for_preview'):decision='NOT_SELECTED';reason='Outside the prioritized stellar/configuration selection'
        else:decision=r['status'];reason=r.get('reason','')
        item['final_selection_status']=decision;item['final_selection_reason']=reason;final_catalog.append(item)
    pd.DataFrame(final_catalog).to_csv('data/catalogs/exposure_catalog.csv',index=False)
    narrow=[r for r in good if r['aperture']=='52X0.5' and r['wing_usable']]
    transfer='NOT VALIDATED: no accepted regular 52X0.5 wing calibrator' if not narrow else 'CONDITIONAL: regular 52X0.5 templates are measured, but universal transfer from 52X2 is not established'
    slit=load_json('products/measurements/slit.json',[]);drizzle=load_json('products/measurements/drizzle_sensitivity.json',[]);epoch=load_json('products/measurements/epoch.json',[]);rect=load_json('products/measurements/rectification.json',[]);background=load_json('products/measurements/background.json',[]);null=load_json('products/measurements/leave_one_star_out.json',[]);sim=load_json('products/measurements/injection_recovery.json',[])
    archived_flts=sum(str(r['paths']['flt']).startswith('data/mast_calibrated/') for r in downloaded)
    summary=f'''# STIS G140L empirical PSF calibration

Working directory: `/home/linxiaojing/data/lrd-squared/hst/stis_psf`.
Pipeline repository: `/quasar/linxiaojing/lrd-squared/codes/hstpipe`.

- Public 2023–2026 G140L archive records: {len(rows)}; previewed: {len(screened)}.
- Selected FLT sets: {len(downloaded)} ({archived_flts} archived, {len(downloaded)-archived_flts} independently recalibrated); accepted hstpipe-native PSF exposures: {len(accepted)}.
- hstpipe reduction groups complete: {sum(r['status']=='COMPLETE' for r in reductions)} / {len(reductions)}.
- Primary native/drizzled PSF profiles: {sum(r['mode'] in ['hstpipe_native','drizzled'] for r in m)}; library groups: {len(library)}.
- Numerical/FITS checks: {audit['status']} ({audit['profiles_checked']} profiles).
- 52X2 → 52X0.5 transfer: **{transfer}**.

Observed narrow-minus-wide EE(11), at common supported normalization radius: {distribution(slit,'ee_11_difference')}. These correlated comparisons do not isolate slit width from epoch/position/centering.

Epoch FWHM differences [native pixel]: {distribution(epoch,'fwhm_difference')}.
Native minus drizzled FWHM [native pixel]: {distribution([r for r in rect if r['comparison']=='native_minus_drizzled'],'fwhm_difference')}.
Background wing-fraction differences: {distribution(background,'wing_fraction_difference')}.
Leave-one-star-out apparent broad-component fractions (0.2 arcsec diagnostic): {distribution([r for r in null if r.get('wing_usable',True)],'extended_fraction')}.

Drizzle pixfrac/scale FWHM shifts [native pixel]: {distribution(drizzle,'fwhm_change')}.
Drizzle pixfrac/scale wing-fraction shifts: {distribution(drizzle,'wing_fraction_change')}.

Library: `psf_library/stis_g140l_psf_library.fits`.
Per-program summary: `reports/program_summary.csv`.
Detailed methods, results and limits: `reports/validation_report.md`.

Reproduce with the existing HST environment:

```sh
hstpipe stis psf --config /home/linxiaojing/data/lrd-squared/hst/stis_psf/config/psf.yaml --stage all
```

This work contains calibration-star products only. No science-target analysis is included.
'''
    Path('reports/summary.md').write_text(summary)
    table=pd.DataFrame(programs).to_string(index=False)
    measures=[]
    grouping=defaultdict(list)
    for r in good:grouping[(r['aperture'],r['band'],r.get('normalization_radius',70))].append(r)
    for (ap,band,radius),g in sorted(grouping.items()):
        measures.append(dict(aperture=ap,band=band,normalization_radius=radius,n_exp=len(g),n_stars=len(set(r['target'] for r in g)),fwhm_median=np.median([r['fwhm'] for r in g]),ee11_median=np.nanmedian([r['ee_11'] for r in g]),wing_median=np.median([r['wing_fraction'] for r in g]),wing_exposure_scatter=np.std([r['wing_fraction'] for r in g],ddof=1) if len(g)>1 else None))
    pd.DataFrame(measures).to_csv('products/measurements/configuration_summary.csv',index=False)
    text=summary+'\n## Per-program execution\n\n```text\n'+table+'\n```\n\n## Measured configuration dependence\n\n```text\n'+pd.DataFrame(measures).to_string(index=False)+'\n```\n'+'''

## Architecture and reduction

The single maintained pipeline is hstpipe. Its existing STIS io, reference,
geometry, calibration, background, coadd, resampling, workflow, extraction and
pipeline modules perform the reduction. New hstpipe.stis.psf modules implement
archive selection, PSF measurements and validation. The scientific working
directory has no independent Python package, environment, or reduction framework.
Existing source modifications from other work were preserved. Every reduction
has a normal hstpipe configuration, stage cache, input checksums, reference
provenance, source hashes, logs and diagnostics. Archive products are immutable.

Reduction groups preserve target, exact aperture, visit, slit angle and coarse
detector position. Registration uses hstpipe's own coordinate operator and measured
continuum trace; recorded source offsets enter hstpipe source_offsets without
editing archive headers. Spectral registration uses the archive wavecal and
HSTCAL wavelength map; no unmeasured spectral shift is invented. hstpipe performs
individual 2D background subtraction, independent x2d controls, official flux
calibration checks, joint detector-to-output drizzle, x2d coaddition, and supported
boxcar/optimal 1D extraction. Existing optimal-extraction branch-flux caveats
remain documented by hstpipe; this PSF library is measured from the 2D products,
not from an optimal-extraction spatial model. DRZ/ERR/DQ/WHT/coverage and background arrays follow
its existing conventions. A completed run is not proof of PSF-wing validity.

Raw-only comparisons are recalibrated separately using stistools and preserved
CRDS provenance. Zero-duration/zero-count records are rejected explicitly. The
existing hstpipe gate requires WAVECORR=COMPLETE; nominal-dispersion raw-only
products without a real wavecal are retained as native diagnostic comparisons,
not falsely marked fully calibrated or forced through that gate.

## Native and drizzled PSF measurement

Primary native profiles use hstpipe's saved background-subtracted FLT-rate arrays,
original ERRBASE and DQ, the full hstpipe two-dimensional observed-wavelength
map, and the measured/calibrated trace coordinate map. Fractional spectral-pixel
band membership accounts for wavelength tilt across the wings. Native pixel
values and variances use hstpipe’s own validated flux-calibration factor, so
comparisons with drizzle share its surface-brightness convention. Strong narrow
stellar emission is excluded using archived X1D screening; the original counts
remain in the immutable FLTs and native pixel stamps.
Native-profile columns must have complete valid-pixel coverage across the full
PSF window; whole columns with a bad spatial sample are excluded rather than
allowing row-dependent missing pixels to bias a faint wing. The accepted and
candidate column counts are recorded. Pixel
overlap conserves flux while aligning subpixel trace offsets and propagates the
redistribution covariance. Original native pixel values, masks, wavelengths and
trace offsets are retained. Overlap alignment can broaden an undersampled core;
unaligned detector sums remain diagnostics, and no smoothing kernel is applied
to suppress wings. X1D supplies coordinates and spectral screening, never a PSF.

The hstpipe BKGERR is a correlated bound; the native profile adds its fully
correlated outer product to the measurement covariance. Drizzled measurements
use output ERRBASE plus the background bound, with a conservative spatial noise
correlation multiplier max(1,pixfrac*native_scale/output_scale). This is not an
exact drizzle covariance reconstruction. Cross-exposure/star scatter and
parameter variations must supplement these conditional errors.

PSF normalization radii are recorded per profile and library group. The largest
complete native window is selected from the configured radii (up to 70 pixels).
Windows smaller than 40 pixels are core-only and cannot receive a wing grade. This prevents the
repeller-wire shadow or a detector boundary from being treated as a zero-flux
wing. Slit comparisons use a common supported radius. Aperture heights exceeding
the measured window are explicitly unsupported (NaN), not assigned unit EE. Aperture heights 3/5/7/9/11/
15/21 and wider supported heights use fractional boundary overlap. PCTAB ratios
are compared in the common default-height normalization EE(h)/EE(11), rather
than equating finite-window fractions with infinite-aperture EE. Direct
half-maximum FWHM and Gaussian core widths are separate diagnostics. Centroid,
asymmetry, EE and wing fractions retain propagated uncertainties. Pixel values
in the FITS profiles are fractions per output bin; offsets and sampling identify
the bin size. Comparison plots convert these to density per native pixel.

Initial archive screening is followed by 2D core, trace, DQ and wing gates.
Core S/N ≥30 and wing S/N ≥5 are distinct requirements. Serious-DQ coverage,
trace residual limits and configurable band coverage govern acceptance. Native
screening uses conservative low-count variance floors. Strong narrow stellar
emission is excluded from continuum screening; stellar absorption/SED weighting
can remain. Automated screening does not prove absence of all field companions;
2D spectra and wing plots require inspection. Rejected profiles remain traceable.

## PSF-only validation and limits

Epoch differences compare the same GRW+70D5824 star at matched bands. Pairwise
comparisons share exposures and are correlated. Measurement-only z values are
labeled as such and do not by themselves establish physical focus/time changes.
Leave-one-star-out tests exclude the tested star, retain template prediction
covariance, and fit centroid/background nuisance terms. The small independent
star sample does not calibrate precise 3- or 5-sigma empirical tail probabilities.

Slit comparisons keep 52X0.5, 52X2, 52X0.2 and D1 configurations separate.
The FITS PROPAPER pseudo-aperture is recorded independently of the hardware
APERTURE; for example, a hardware 52X0.2 exposure requested at 52X0.2D1 is
classified as D1, not as regular 52X0.2. Detector position alone is not used to
infer a pseudo-aperture. Different stars, epochs, across-slit centering and detector locations remain
confounders. No universal slit correction is manufactured. Broad integrations are
source-spectrum and mask weighted; they are not universal 1330–1680 Angstrom
responses. Prefer matched narrow bands. A broad-band template needs compatible
spectral weights and validated sampling; no science-target reweighting or
matching is performed here. A high-S/N matched
52X0.5 calibrator is preferable to transferring a wide-slit PSF. The 18113 brown
dwarf is evaluated independently and cannot qualify from program membership.

Rectification comparisons require both profiles to pass their stated QA.
Comparisons to archived x2d can also include calibration-version differences;
the native-to-joint-drizzle comparison shares the current hstpipe calibration.
Old CAL_VER dates identify archive software, not the observation epoch.
Rectification comparisons quantify native versus standard x2d and hstpipe
joint-drizzle FWHM, EE, wings and asymmetry. Native covariance from explicit
pixel overlap is conditional on the fixed fitted trace. x2d ERR omits the full
rectification covariance. Drizzle-parameter sensitivity is reported separately
and is not replaced by the native-to-drizzled comparison.

Alternate off-source constant/linear backgrounds test wings with |y|=100–160
and 120–180 native pixels. These profiles are diagnostics rather than a second
reduction pipeline. Their spread supplements the hstpipe background uncertainty. Wing-variation
statistics are restricted to wing-qualified profiles and a common spatial window.
DARKCORR does not remove all sky/geocoronal emission or detector background
structure. Ly-alpha/O I airglow and configurable instrumental intervals are masked.

Library templates are stratified by wavelength, exact aperture, epoch, sampling,
and coarse detector position. Repeats first combine within each star, then stars
have equal weight. Mean and prediction covariance include measured repeat variability and
inter-star scatter. Raw sample scatter is not photon-noise-debiased and is
therefore conservative; it does not assign zero template variability merely
because a configuration group contains only one star. One-star groups explicitly lack an independent
inter-star scatter estimate. Coarse position groups are not proof of stability
at arbitrary detector positions. Interpolation is allowed only between bracketed
calibrated bands in a matching group; unconstrained extrapolation is prohibited.

Generic injection/recovery diagnostics, where present, convolve intrinsic
exponential profiles with measured PSFs. They test the extraction/null behavior
under Gaussian count-limit noise plus empirical held-out residuals. Their
conditional completeness and false-positive rates do not establish performance
for arbitrary sources or low-count event data. They are not target-specific.

## Recommended use

Use a high-grade template matched in aperture, wavelength, detector position and
epoch. Retain the individual native/drizzled profiles and predictive covariance.
Use empirical between-star scatter plus the measured background and drizzle
sensitivity when interpreting wings. A single-star group's photon errors cannot
stand in for template systematics. Where slit transfer remains conditional,
obtain repeated matched-aperture standards with controlled centering and sampled
subpixel phases. Focus and unknown across-slit offsets remain limitations.

## References consulted

- [Cycle 33 STIS Calibration Plan](https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/stis/calibration/_documents/calibration_files/Cycle_33/STIS_cyc33_annual_cal_approved_public.pdf)
- [STIS ISR 2026-01](https://arxiv.org/abs/2604.21966): sensitivity trends do not measure spatial PSF stability.
- [First-order spatial profiles](https://hst-docs.stsci.edu/stisihb/chapter-13-spectroscopic-reference-material/13-5-spatial-profiles/first-order-spatial-profiles)
- [Working with spectral images](https://hst-docs.stsci.edu/stisdhb/chapter-5-stis-data-analysis/5-4-working-with-spectral-images)
- [Calibration steps](https://hst-docs.stsci.edu/stisdhb/chapter-3-stis-calibration/3-4-descriptions-of-calibration-steps)
- [Point-source extraction/PCTAB](https://hst-docs.stsci.edu/stisdhb/chapter-5-stis-data-analysis/5-5-working-with-extracted-spectra)
'''
    text+='\n## Drizzle parameter sensitivity\n\n```text\n'+pd.DataFrame(drizzle).to_string(index=False)+'\n```\n'
    text+='\n## Accepted dataset dates and configurations\n\n```json\n'+json.dumps(load_json('reports/dataset_census.json',{}),indent=2)+'\n```\n'
    text+='\n## Explicit quality gates\n\n```text\n'+pd.DataFrame(gates).to_string(index=False)+'\n```\n'
    text+='\n## Integrated-wing predictive uncertainties\n\nThe ECSV table contains group-specific errors including repeat and inter-star scatter. CORE-only groups do not calibrate unsupported outer wings.\n\n```text\n'+pd.DataFrame(wing_uncertainties).to_string(index=False)+'\n```\n'
    text+='\n## Recalibration execution\n\n```json\n'+json.dumps([{k:v for k,v in r.items() if k!='references'} for r in load_json('data/catalogs/recalibration.json',[])],indent=2)+'\n```\n'
    text+='\n## Reduction failures or unsupported cases\n\n```json\n'+json.dumps([r for r in reductions if r['status']!='COMPLETE'],indent=2)+'\n```\n'
    text+='\n## Recommended library groups\n\n```json\n'+json.dumps(library,indent=2)+'\n```\n'
    import shutil
    for red in reductions:
        if red.get('status')!='COMPLETE':continue
        image=Path(red['directory'])/'diagnostics/background'/(red['roots'][0]+'_flt.png')
        if image.exists():shutil.copy2(image,'figures/hstpipe_background_example.png');break
    text+='\n## Diagnostic figures\n\n'
    figures=[('selection.png','Archive selection, dates, and exact requested aperture configurations'),('hstpipe_background_example.png','Representative hstpipe two-dimensional background reduction'),('psf_1420-1480.png','Native/drizzled core, logarithmic wings, and encircled energy at common measured 11-pixel normalization'),('wavelength.png','FWHM and wavelength dependence'),('detector_position.png','Detector-position dependence; other variables may be confounded'),('epoch.png','Matched-band temporal comparisons; pairs share exposures'),('slit.png','Wing-qualified slit comparisons at supported common normalization radii'),('rectification.png','Native versus rectified/drizzled profile differences'),('background.png','Alternative off-source background models'),('leave_one_star_out.png','Independent-star point-source null fits'),('drizzle_sensitivity.png','Measured pixfrac and output-scale sensitivity'),('injection_recovery.png','Generic calibration-algorithm injection/null diagnostics')]
    for filename,caption in figures:
        if Path('figures',filename).exists():text+=f'### {caption}\n\n![{caption}](../figures/{filename})\n\n'
    comparison=load_json('reports/recalibration_baseline_comparison.json',[])
    text+='\n## Independent native-count recalibration check\n\nOld CAL_VER strings identify archive software, not the observation epoch. The primary templates are observations from 2023–2026.\n\n```json\n'+json.dumps(comparison,indent=2)+'\n```\n'
    Path('reports/validation_report.md').write_text(text)
    save_json('reports/execution_status.json',dict(status='EXECUTED_WITH_LIMITATIONS',host=platform.node(),archive_records=len(rows),downloaded_flt_sets=archived_flts,available_flt_sets=len(downloaded),accepted_native_exposures=len(accepted),reduction_groups_complete=sum(r['status']=='COMPLETE' for r in reductions),library_groups=len(library),numerical_audit=audit,slit_transfer=transfer))
