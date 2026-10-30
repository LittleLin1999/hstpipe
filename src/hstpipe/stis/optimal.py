"""Empirical spatial-profile extraction without a Gaussian PSF assumption.

The profile is estimated from wavelength-averaged, background-subtracted data.
Only the profile model is constrained to be nonnegative; science pixels retain
their signed values. The estimator uses equal noise weights because zero or
underestimated low-count MAMA ERR values are unsafe as inverse-variance weights.
It is optimal for a known profile with spatially constant independent noise.
For drizzled data, a conservative positive covariance bound is reported.
"""

import numpy as np


def empirical_profile(science, valid, wavelength, aperture, return_metadata=False):
    """Return one normalized, empirical cross-dispersion template.

    Airglow wavelengths are omitted from template estimation. This has no
    effect on the science values extracted at those wavelengths.
    """
    wave = np.asarray(wavelength, float)
    quiet = np.isfinite(wave)
    for lo, hi in ((1210, 1223), (1298, 1310), (1348, 1364)):
        quiet &= ~((wave >= lo) & (wave <= hi))
    if quiet.sum() < 10:
        quiet = np.isfinite(wave)
    good = valid & quiet[None, :] & (aperture > 0)
    count = good.sum(axis=1)
    mean = np.divide(
        np.where(good, science, 0).sum(axis=1),
        count,
        out=np.zeros(science.shape[0], float),
        where=count > 0,
    )
    profile = np.maximum(mean, 0)
    profile[~np.any(aperture > 0, axis=1)] = 0
    fallback = not np.isfinite(profile).all() or profile.sum() <= 0
    if fallback:
        # No empirical source detection: equal spatial weights reduce to a
        # boxcar. The caller should flag this in its quality-control report.
        profile = np.any(aperture > 0, axis=1).astype(float)
    profile = profile / profile.sum()
    metadata = dict(fallback_to_flat=bool(fallback), template_columns=int(quiet.sum()),
                    method="EMPIRICAL_EQUALVAR", airglow_excluded=True)
    return (profile, metadata) if return_metadata else profile


def profile_extract(science, stat_error, bkg_error, valid, aperture, profile,
                    covariance_bound=False, min_profile_fraction=0.8,
                    min_information_fraction=0.5):
    """Fit full-aperture flux using only valid pixels of a fixed profile.

    As in unfold_jwst.extract_optimal, normalize the model BEFORE applying
    the pixel mask. Renormalizing on the surviving pixels would underestimate
    the full-aperture flux. Here q = P*A / sum(P*A), D = SCI*A, and
    F = sum_valid(q*D) / sum_valid(q**2). Equal-noise weights avoid near-zero
    MAMA pipeline variances. Negative data and zero ERR remain valid inputs.

    Missing flux is inferred from the profile, not filled in the image. Require
    at least 80% of profile light and 50% of equal-noise fitting information by
    default; these conservative thresholds are configurable and recorded.
    Statistical variance uses squared estimator coefficients, including the
    increased uncertainty after masking. Profile-shape uncertainty is separate
    and is NOT included in these conditional errors.
    """
    science, stat_error, bkg_error = map(
        lambda x: np.asarray(x, float), (science, stat_error, bkg_error))
    area = np.broadcast_to(np.asarray(aperture, float), science.shape)
    shape = np.asarray(profile, float)[:, None] * area
    norm = shape.sum(axis=0)  # full geometrical aperture, independent of DQ
    shape = np.divide(shape, norm, out=np.zeros_like(shape), where=norm > 0)
    good = (np.asarray(valid, bool) & np.isfinite(science)
            & np.isfinite(stat_error) & (stat_error >= 0)
            & np.isfinite(bkg_error) & (bkg_error >= 0))
    fitted_shape = np.where(good, shape, 0.)
    denominator = np.sum(fitted_shape**2, axis=0)
    profile_fraction = fitted_shape.sum(axis=0)
    information = np.divide(denominator, np.sum(shape**2, axis=0),
                            out=np.zeros_like(denominator), where=np.sum(shape**2, axis=0)>0)
    coefficient = np.divide(fitted_shape, denominator,
                           out=np.zeros_like(shape), where=denominator > 0)
    weighted = coefficient * area
    flux = np.sum(np.where(good, science, 0) * weighted, axis=0)
    if covariance_bound:
        stat = np.sum(np.where(good, stat_error, 0) * abs(weighted), axis=0)
    else:
        stat = np.sqrt(np.sum((np.where(good, stat_error, 0) * weighted)**2, axis=0))
    model = np.sum(np.where(good, bkg_error, 0) * abs(weighted), axis=0)
    accepted = ((norm > 0) & (denominator > 0)
                & (profile_fraction >= min_profile_fraction - 1e-12)
                & (information >= min_information_fraction - 1e-12))
    return {
        "FLUX_OPT": np.where(accepted, flux, np.nan),
        "ERR_OPT": np.where(accepted, np.hypot(stat, model), np.nan),
        "ERR_STAT_OPT": np.where(accepted, stat, np.nan),
        "BKGERR_OPT": np.where(accepted, model, np.nan),
        "DQ_OPT": np.where(accepted, 0, 4).astype("u2"),
        "PROFCOV_OPT": profile_fraction,
        "INFOFRAC_OPT": information,
        "OPT_MASKED": np.any(~good & (area > 0), axis=0),
    }


def profile_fit_diagnostics(science, stat_error, valid, aperture, profile, flux):
    """Exact forward model of profile_extract, including fractional edge areas.

    D = SCI * aperture; M = flux * normalized(profile * aperture).
    Report diagonal chi-square ONLY as a descriptive statistic. It is not the
    equal-noise fitting objective, ignores covariance, and treats the learned
    template as fixed. Nonpositive ERR pixels remain in extraction but cannot
    define standardized residuals. Their diagnostic mask is recorded explicitly.
    """
    area = np.broadcast_to(np.asarray(aperture, float), science.shape)
    active = valid & (area > 0)
    shape = np.asarray(profile)[:, None] * area
    norm = shape.sum(axis=0)
    shape = np.divide(shape, norm, out=np.zeros_like(shape), where=norm > 0)
    model_weighted = shape * np.asarray(flux)[None, :]
    model = np.divide(model_weighted, area, out=np.full_like(model_weighted, np.nan), where=area > 0)
    full_model = np.divide(np.asarray(profile)[:,None]*np.asarray(flux)[None,:], norm,
                           out=np.full_like(model,np.nan),where=norm>0)
    residual = np.where(active, science - model, np.nan)
    diagonal_ok = active & np.isfinite(stat_error) & (stat_error > 0) & np.isfinite(residual)
    chi = np.divide(residual, stat_error, out=np.full_like(residual, np.nan), where=diagonal_ok)
    n = diagonal_ok.sum(axis=0)
    chi2 = np.where(n > 0, np.nansum(chi**2, axis=0), np.nan)
    dof = n - 1  # nominal: one amplitude, conditioned on an estimated template
    reduced = np.divide(chi2, dof, out=np.full_like(chi2, np.nan), where=dof > 0)
    return dict(MODEL=model, MODEL_FULL=full_model, PROFILE_FULL=np.broadcast_to(np.asarray(profile)[:,None],science.shape), RESIDUAL=residual, CHI_DIAG=chi, PROFILE=shape,
                DIAGMASK=diagonal_ok, CHI2_DIAG=chi2, NDOF_NOMINAL=dof,
                CHI2_PER_DOF=reduced,
                ERROR_UNUSABLE_PIXELS=int(np.sum(active & ~diagonal_ok)),
                FULL_COLUMNS=np.all(valid | (area == 0), axis=0),
                ACCEPTED_COLUMNS=np.isfinite(flux))


def diagnostic_noise_scale(science, stat_error, valid, spatial, bin_columns=32):
    """Off-source noise floor for DISPLAY ONLY; never change extraction errors.

    Unclipped variance in |y|=100..150 native pixels, in fixed wavelength blocks,
    supplies a floor where per-pixel MAMA ERR is tiny. This can include residual
    sky structure and has no Gaussian/Poisson likelihood interpretation. No
    fitting residuals are used, so a poor source model cannot tune its own noise.
    """
    sky = valid & (abs(np.asarray(spatial)[:, None]) >= 100) & (abs(np.asarray(spatial)[:, None]) < 150)
    floor = np.full(science.shape[1], np.nan)
    for left in range(0,science.shape[1],bin_columns):
        right=min(left+bin_columns,science.shape[1]);m=sky[:,left:right]
        values=science[:,left:right][m]
        if len(values)<32: continue
        empirical=float(np.std(values,ddof=1))
        pipeline=float(np.sqrt(np.mean(stat_error[:,left:right][m]**2)))
        floor[left:right]=max(empirical,pipeline)
    return np.maximum(stat_error,floor[None,:]), floor


def continuum_center(science, valid, wavelength, spatial, guess=0., native_scale=.0246):
    """Signed broadband core centroid; no Gaussian PSF or faint-line registration.

    Exclude airglow and the strongest calibrated non-airglow feature (+/-10 A).
    Bootstrap 32-column blocks to screen centering precision. Only extraction
    apertures move; the image registration and dither coordinate grid do not.
    """
    wave=np.asarray(wavelength);y=(np.asarray(spatial)-guess)/native_scale
    quiet=np.isfinite(wave)&np.any(valid,axis=0)
    for lo,hi in ((1195,1235),(1285,1325),(1340,1375)):
        quiet &= ~((wave>=lo)&(wave<=hi))
    core=abs(y)<10
    signal=np.where(valid[core],science[core],0).sum(axis=0)
    candidates=np.flatnonzero(quiet)
    if len(candidates)<100:raise ValueError('Too little calibrated continuum for centering')
    peak=int(candidates[np.argmax(signal[candidates])]);quiet &= abs(wave-wave[peak])>10
    selected=valid&quiet[None,:]
    n=selected.sum(axis=1)
    collapsed=np.divide(np.where(selected,science,0).sum(axis=1),n,out=np.zeros(len(y)),where=n>0)
    center=float(y[np.flatnonzero(core)[np.argmax(collapsed[core])]])
    for _ in range(4):
        region=(abs(y-center)<6)&core
        denominator=collapsed[region].sum()
        if denominator<=0:raise ValueError('No positive broadband source for centering')
        center=float(np.sum(y[region]*collapsed[region])/denominator)
    blocks=[np.arange(i,min(i+32,len(wave))) for i in range(0,len(wave),32)]
    blocks=[b[quiet[b]] for b in blocks if quiet[b].any()]
    sums=np.array([np.where(valid[:,b],science[:,b],0).sum(axis=1) for b in blocks])
    counts=np.array([valid[:,b].sum(axis=1) for b in blocks])
    rng=np.random.default_rng(522);draws=[]
    for _ in range(128):
        pick=rng.integers(0,len(blocks),len(blocks));ss=sums[pick].sum(axis=0);nn=counts[pick].sum(axis=0)
        profile=np.divide(ss,nn,out=np.zeros_like(ss),where=nn>0)
        denom=profile[region].sum()
        if denom>0:draws.append(float(np.sum(y[region]*profile[region])/denom))
    error=float(np.std(draws,ddof=1)) if len(draws)>100 else np.inf
    if abs(center)>5 or error>.5:
        raise ValueError(f'Continuum center uncertain: shift={center:.3f}, block error={error:.3f} native pixels; configure header center explicitly')
    return dict(method='signed_continuum_core_centroid',center_arcsec=guess+center*native_scale,
                shift_native_pixels=center,block_error_native_pixels=error,
                excluded_feature_wavelength=float(wave[peak]),continuum_columns=int(quiet.sum()),
                limitation='Core centroid, not total-light centroid; block error omits full covariance')


def wide_profile(science, valid, wavelength, fit_region, spatial_step=1.):
    """Nonnegative nonparametric profile on a region independent of extraction.

    Select the weakest curvature penalty within 5% of the best two-fold
    wavelength-block prediction score. Signed collapsed data enter the fit;
    positivity applies only to the model. No symmetry or Gaussian is imposed.
    """
    from scipy.optimize import nnls
    wave=np.asarray(wavelength);quiet=np.isfinite(wave)
    for lo,hi in ((1210,1223),(1298,1310),(1348,1364)):
        quiet &= ~((wave>=lo)&(wave<=hi))
    region=np.broadcast_to(fit_region,science.shape)>0
    good=valid & quiet[None,:] & region
    rows=np.flatnonzero(np.any(region,axis=1))
    if len(rows)<5:raise ValueError('Profile fitting region too narrow')
    def collapse(mask):
        n=mask[rows].sum(axis=1)
        return np.divide(np.where(mask[rows],science[rows],0).sum(axis=1),n,
                         out=np.zeros(len(rows)),where=n>0), n>0
    mean,supported=collapse(good)
    scale=max(float(np.max(abs(mean))),np.finfo(float).tiny)
    curvature=np.diff(np.eye(len(rows)),n=2,axis=0)/spatial_step**2
    def fit(values,ok,penalty):
        design=np.vstack((np.eye(len(rows))[ok],np.sqrt(penalty)*curvature))
        target=np.r_[values[ok]/scale,np.zeros(len(curvature))]
        return nnls(design,target)[0]*scale
    candidates=(0.,.01,.1,1.)
    scores={p:0. for p in candidates}
    for fold in (0,1):
        train,ok=collapse(good & (((np.arange(len(wave))//32)%2)!=fold)[None,:])
        held,hk=collapse(good & (((np.arange(len(wave))//32)%2)==fold)[None,:])
        if np.sum(ok&hk)<5:raise ValueError('Insufficient independent profile-training blocks')
        for penalty in candidates:
            prediction=fit(train,ok,penalty)
            scores[penalty]+=float(np.mean(((prediction[hk]-held[hk])/scale)**2))
    best=min(scores.values());penalty=min(p for p in candidates if scores[p]<=1.05*best+1e-20)
    fitted=fit(mean,supported,penalty)
    fallback=not np.isfinite(fitted).all() or fitted.sum()<=0
    if fallback:raise ValueError('No constrained positive wide profile; use boxcar output or an explicit source model')
    profile=np.zeros(science.shape[0]);profile[rows]=fitted;profile/=profile.sum()
    return profile,dict(method='EMPIRICAL_EQUALVAR',fallback_to_flat=False,
        template_columns=int(np.sum(quiet&np.any(valid,axis=0))),airglow_excluded=True,
        fitting_rows=[int(rows[0]),int(rows[-1])],curvature_penalty=penalty,
        validation_scores=scores,profile_fitting_independent_of_boxcar=True,
        limitation='Empirical wavelength-averaged shape; unresolved wings remain uncertain')
