"""All residual-sky estimation, spatial model selection and uncertainty logic.

Signed means in native count-rate space; source mask, no ERR-based selection,
no clipping. Sparse airglow cells borrow the SAME wavelength from spatial
neighbors. This module never resamples SCI or changes calibration.
"""

from dataclasses import dataclass, replace
import numpy as np
from .config import BackgroundConfig

AIRGLOW = ((1195.0, 1235.0), (1285.0, 1325.0), (1340.0, 1375.0))


def wavelength_edges(wave, continuum_bin=20.0, dispersion=0.584):
    """Derive coverage from data; airglow wavelengths are instrument settings."""
    lo = np.floor(np.nanmin(wave) / continuum_bin) * continuum_bin
    hi = np.ceil(np.nanmax(wave) / continuum_bin) * continuum_bin
    edges = list(np.arange(lo, hi, continuum_bin)) + [hi]
    for left, right in AIRGLOW:
        if right <= lo or left >= hi:
            continue
        left = max(left, lo)
        right = min(right, hi)
        edges = [v for v in edges if not left < v < right]
        edges += list(np.arange(left, right, dispersion)) + [right]
    return np.unique(edges)


def spectra(counts, wave, sy, valid, edges, lo, hi):
    """Pixel sums/counts in wavelength bins and independent 8-row clusters."""
    mask = valid & (sy >= lo) & (sy < hi)
    bi = np.searchsorted(edges, wave, side="right") - 1
    region = mask & (bi >= 0) & (bi < len(edges) - 1)
    nbin = len(edges) - 1
    indices = bi[region]
    values = counts[region]
    n = np.bincount(indices, minlength=nbin)
    summ = np.bincount(indices, weights=values, minlength=nbin)
    squares = np.bincount(indices, weights=values * values, minlength=nbin)
    mean = np.divide(summ, n, out=np.full(nbin, np.nan), where=n > 0)
    pixelvar = np.divide(
        squares - np.divide(summ * summ, n, out=np.zeros(nbin), where=n > 0),
        n - 1,
        out=np.full(nbin, np.nan),
        where=n > 1,
    )
    # Cluster on detector Y, not interpolated image samples. DQ/coverage only.
    row_labels = np.nonzero(region)[0] // 8
    labels = np.unique(row_labels)
    sums = []
    nums = []
    for label in labels:
        selected = row_labels == label
        ii = indices[selected]
        sums.append(np.bincount(ii, weights=values[selected], minlength=nbin))
        nums.append(np.bincount(ii, minlength=nbin))
    sums = np.array(sums)
    nums = np.array(nums)
    return dict(
        mean=mean,
        n=n,
        pixel_se=np.sqrt(np.maximum(pixelvar, 0) / np.maximum(n, 1)),
        sums=sums,
        nums=nums,
        mask=mask,
    )


def draws(spec, rng, reps=160):
    k = len(spec["sums"])
    nb = len(spec["mean"])
    if k < 4:
        return np.full((reps, nb), np.nan)
    pick = rng.integers(0, k, size=(reps, k))
    nn = spec["nums"][pick].sum(1)
    ss = spec["sums"][pick].sum(1)
    return np.divide(ss, nn, out=np.full_like(ss, np.nan), where=nn > 0)


def fit(
    rate,
    variance,
    wave,
    sy,
    valid,
    inner=100.0,
    width=50.0,
    prediction_valid=None,
    options=None,
    edges=None,
):
    """Fit coarse signed means; return continuous model and correlation bound.

    Only DQ/coverage/finite SCI select pixels. ERR is an uncertainty floor, never
    a fit weight or selection criterion. 8-detector-row bootstrap uncertainty
    preserves native sampling. No clipping or positive floor is applied.
    """
    options = options or BackgroundConfig()
    edges = (
        wavelength_edges(wave[valid], options.continuum_bin) if edges is None else edges
    )
    centers = (edges[:-1] + edges[1:]) / 2
    intervals = []
    for sign in [-1, 1]:
        for distance in np.arange(inner, np.max(abs(sy[valid])) + width, width):
            limits = sorted([sign * distance, sign * (distance + width)])
            if np.count_nonzero(valid & (sy >= limits[0]) & (sy < limits[1])) > 100:
                intervals.append(limits)
    anchors, models, errors, records = [], [], [], []
    rng = np.random.default_rng(14023)
    for lo, hi in sorted(intervals):
        spec = spectra(rate, wave, sy, valid, edges, lo, hi)
        usable = (spec["n"] >= 12) & np.isfinite(spec["mean"])
        if usable.sum() < 10:
            continue
        mask = spec["mask"]
        anchor = float(np.mean(sy[mask]))
        boot = draws(spec, rng, reps=options.bootstrap_draws)
        error = np.fmax(np.nanstd(boot, axis=0, ddof=1), spec["pixel_se"])
        vb = spectra(variance, wave, sy, valid, edges, lo, hi)
        floor = np.sqrt(np.maximum(vb["mean"], 0) / np.maximum(spec["n"], 1))
        error = np.fmax(error, floor)
        usable &= np.isfinite(error)
        anchors.append(anchor)
        models.append(np.where(usable, spec["mean"], np.nan))
        errors.append(np.where(usable, error, np.nan))
        records.append(dict(rows=[lo, hi], anchor=anchor, valid_bins=int(usable.sum())))
    anchors = np.asarray(anchors)
    if not np.any(anchors < 0) or not np.any(anchors > 0):
        raise ValueError("No sky anchors on both sides of source")
    # Interpolation is applied to the fitted model only; SCI is never resampled.
    order = np.argsort(anchors)
    anchors = anchors[order]
    models, errors = np.array(models)[order], np.array(errors)[order]
    # Missing edge cells must borrow the SAME wavelength from spatial neighbors.
    # Interpolating their spectrum first can erase airglow when every fine line
    # cell fails the sample-count threshold but broad continuum cells survive.
    available = np.isfinite(models) & np.isfinite(errors)
    supported = available.any(axis=0)
    lower_support = np.full(len(centers), np.nan)
    upper_support = np.full(len(centers), np.nan)
    for k in np.flatnonzero(supported):
        ok = available[:, k]
        lower_support[k], upper_support[k] = anchors[ok][[0, -1]]
        models[:, k] = np.interp(anchors, anchors[ok], models[ok, k])
        errors[:, k] = np.interp(anchors, anchors[ok], errors[ok, k])
    # Entirely uncovered wavelength bins (including beyond the wavelength range)
    # have no spatial estimate. Fill for numerical evaluation only, flag below.
    if not supported.any():
        raise ValueError("No wavelength has a constrained sky estimate")
    for j in range(len(anchors)):
        models[j] = np.interp(centers, centers[supported], models[j, supported])
        errors[j] = np.interp(centers, centers[supported], errors[j, supported])
    lower_support = np.interp(centers, centers[supported], lower_support[supported])
    upper_support = np.interp(centers, centers[supported], upper_support[supported])
    idx = np.clip(np.searchsorted(anchors, sy) - 1, 0, len(anchors) - 2)
    fraction = np.clip((sy - anchors[idx]) / (anchors[idx + 1] - anchors[idx]), 0, 1)
    prediction_valid = valid if prediction_valid is None else prediction_valid
    model = np.zeros(wave.shape)
    bound = np.zeros(wave.shape)
    for j in range(len(anchors) - 1):
        use = prediction_valid & (idx == j)
        f = fraction[use]
        model[use] = (1 - f) * np.interp(wave[use], centers, models[j]) + f * np.interp(
            wave[use], centers, models[j + 1]
        )
        # Positive interpolation of standard errors is a correlation upper bound.
        bound[use] = (1 - f) * np.interp(wave[use], centers, errors[j]) + f * np.interp(
            wave[use], centers, errors[j + 1]
        )
    outside = prediction_valid & (
        (sy < np.interp(wave, centers, lower_support))
        | (sy > np.interp(wave, centers, upper_support))
        | (np.interp(wave, centers, supported.astype(float)) < 1)
        | (wave < centers[0])
        | (wave > centers[-1])
    )
    return model, bound, outside, records


@dataclass
class BackgroundResult:
    rate: np.ndarray
    error: np.ndarray
    flags: np.ndarray
    metadata: dict
    components: dict = None


def scoring_cells(rate, variance, wave, sy, mask):
    """Fixed target-independent cells for comparing held-out predictions."""
    if not mask.any():
        return []
    result = []
    for lo in np.arange(np.floor(wave[mask].min() / 50) * 50, wave[mask].max(), 50):
        band = mask & (wave >= lo) & (wave < lo + 50)
        values = rate[band]
        if len(values) < 40:
            continue
        pool = float(np.var(values, ddof=1))
        for bottom in np.arange(np.floor(sy[band].min() / 50) * 50, sy[band].max(), 50):
            use = band & (sy >= bottom) & (sy < bottom + 50)
            n = use.sum()
            if n < 40:
                continue
            v = max(
                float(np.var(rate[use], ddof=1)),
                float(np.mean(variance[use])),
                pool * 0.1,
                np.finfo(float).tiny,
            )
            result.append((use, np.sqrt(v / n)))
    return result


def select_width(rate, variance, wave, sy, valid, options, edges):
    """Two-fold native-row validation; choose coarsest model within 5% of best.

    These are prediction scores, not formal chi-square probabilities. The
    source and guard pixels participate in neither fitting nor selection.
    """
    widths = options.spatial_widths
    if len(widths) == 1:
        return widths[0], {}
    blocks = np.indices(rate.shape)[0] // 8
    scores = {w: [] for w in widths}
    for fold in (0, 1):
        train = valid & (blocks % 2 != fold)
        held = valid & (blocks % 2 == fold) & (abs(sy) >= options.sky_inner)
        cells = scoring_cells(rate, variance, wave, sy, held)
        if not cells:
            raise ValueError("Insufficient independent sky for model selection")
        for w in widths:
            model, _, _, _ = fit(
                rate,
                variance,
                wave,
                sy,
                train,
                inner=options.sky_inner,
                width=w,
                prediction_valid=valid,
                options=options,
                edges=edges,
            )
            scores[w].extend(
                (float(np.mean((rate - model)[mask])) / scale) ** 2
                for mask, scale in cells
            )
    scores = {w: float(np.mean(v)) for w, v in scores.items()}
    best = min(scores.values())
    return max(w for w in widths if scores[w] <= 1.05 * best), scores


def estimate(exposure, wave, sy, valid, options, dispersion=0.584):
    """Return a separate rate model for one original count-valued exposure."""
    rate = exposure.counts / exposure.time
    variance = exposure.variance / exposure.time**2
    if not options.enabled:
        return BackgroundResult(
            np.zeros_like(rate),
            np.zeros_like(rate),
            np.zeros(rate.shape, "u2"),
            {"enabled": False},
        )
    edges = wavelength_edges(wave[valid], options.continuum_bin, dispersion)
    width, scores = select_width(rate, variance, wave, sy, valid, options, edges)
    model, error, extrap, anchors = fit(
        rate,
        variance,
        wave,
        sy,
        valid,
        inner=options.sky_inner,
        width=width,
        options=options,
        edges=edges,
    )
    recovery = []
    for sigma in (20.0, 40.0, 80.0):
        injected = 1e-4 * np.exp(-0.5 * (sy / sigma) ** 2)
        recovered, _, _, _ = fit(
            rate + injected,
            variance,
            wave,
            sy,
            valid,
            inner=options.sky_inner,
            width=width,
            options=options,
            edges=edges,
        )
        roi = valid & (abs(sy) <= options.source_radius)
        recovery.append(
            {
                "sigma_native_pixels": sigma,
                "recovered_fraction": float(
                    np.sum((injected - recovered + model)[roi]) / np.sum(injected[roi])
                ),
            }
        )
    return BackgroundResult(
        model,
        error,
        np.where(valid, np.where(extrap, 2, 0), 1).astype("u2"),
        dict(
            enabled=True,
            selected_width=width,
            scores=scores,
            anchors=anchors,
            halo_recovery=recovery,
            wavelength_edges=edges.tolist(),
            uncertainty="Spatial block bootstrap plus sample/ERR floors; correlated interpolation upper bound",
            limitation="Source flux reaching sky strips is degenerate with background",
        ),
    )


def assess_x1d_residual(rate, dq, covered, trace, wavelength,
                        calstis_background, aperture_height, offset, size):
    """Audit CALSTIS's residual x1d sky against unclipped MAMA pixel means.

    CALSTIS cs6 performs iterative 3-sigma clipping independently in every
    wavelength column. For sparse dark-subtracted FUV-MAMA data this can
    reject the rare positive photon pixels and approach the negative dark
    pedestal. We therefore retain the official result as a diagnostic and
    accept its correction only if it agrees with an independent, unclipped
    mean of the *same* two off-source windows. Wavelength blocks, rather than
    individual pixels, set the uncertainty to acknowledge spectral covariance.
    """
    rate = np.asarray(rate, float)
    dq = np.asarray(dq)
    covered = np.asarray(covered, bool)
    trace = np.asarray(trace, float)
    wavelength = np.asarray(wavelength, float)
    if rate.ndim != 2 or rate.shape != dq.shape or rate.shape != covered.shape:
        raise ValueError("Residual-sky audit requires matching detector arrays")
    if trace.shape != (rate.shape[1],) or wavelength.shape != trace.shape:
        raise ValueError("Residual-sky audit requires one trace/wavelength per column")
    quiet = np.isfinite(wavelength)
    low, high = np.quantile(wavelength[quiet], [0.05, 0.95])
    quiet &= (wavelength >= low) & (wavelength <= high)
    for left, right in AIRGLOW:
        quiet &= ~((wavelength >= left) & (wavelength <= right))
    rows = np.arange(rate.shape[0])[:, None]
    sky = np.zeros(rate.shape, bool)
    for sign in (-1, 1):
        sky |= abs(rows - (trace[None, :] + sign * offset)) <= size / 2
    sky &= (dq == 0) & covered & np.isfinite(rate)
    columns = np.flatnonzero(quiet & (sky.sum(axis=0) >= size))
    if len(columns) < 100:
        raise ValueError("Too few independently covered continuum sky columns")
    selected = sky[:, columns]
    mean = float(np.mean(rate[:, columns][selected]))
    blocks = np.array_split(columns, max(8, int(np.ceil(len(columns) / 30))))
    block_means = np.array([
        np.mean(rate[:, block][sky[:, block]])
        for block in blocks if np.any(sky[:, block])
    ])
    sem = float(np.std(block_means, ddof=1) / np.sqrt(len(block_means)))
    hstcal = float(np.median(
        np.asarray(calstis_background, float)[columns] / aperture_height
    ))
    discrepancy = abs(hstcal - mean) / max(sem, np.finfo(float).eps)
    return {
        "unclipped_mean_rate": mean,
        "block_standard_error_rate": sem,
        "calstis_median_rate": hstcal,
        "difference_sigma": discrepancy,
        "continuum_columns": int(len(columns)),
        "reject_calstis_residual": bool(discrepancy > 5),
    }


def residual_measurement(rate, variance, primary_error, mask, row_block=4):
    """Unclipped rate mean with empirical 4-row x 32-column cluster error.

    Original ERR supplies only an uncertainty floor, never weights or selection.
    The primary-model error is a fully correlated bound: it does not decrease
    with sky-pixel count. Adding bounds linearly avoids assuming independence
    from a primary model inferred from the same exposure. These are screening
    scales, not calibrated Gaussian p-values.
    """
    yy, xx = np.nonzero(mask)
    values = rate[mask]
    if len(values) < 32:
        return dict(n=int(len(values)), mean=None, error=None, reason="insufficient pixels")
    mean = float(values.mean())
    groups, inverse = np.unique((yy // row_block) * (rate.shape[1] // 32 + 1) + xx // 32,
                                return_inverse=True)
    if len(groups) < 4:
        return dict(n=int(len(values)), mean=mean, error=None, reason="insufficient blocks")
    totals = np.bincount(inverse, weights=values - mean)
    cluster = np.sqrt(len(groups) / (len(groups) - 1) * np.sum(totals**2)) / len(values)
    sample = np.std(values, ddof=1) / np.sqrt(len(values))
    floor = np.sqrt(np.sum(np.nan_to_num(variance[mask]))) / len(values)
    bound = float(np.mean(primary_error[mask]))
    error = float(max(cluster, sample, floor) + bound)
    return dict(n=int(len(values)), blocks=int(len(groups)), mean=mean, error=error,
                primary_bound=bound, random_error=float(max(cluster, sample, floor)),
                significance=mean / error if error > 0 else None)


def residual_gate(measurements, options):
    """Require training detection, independent validation and two-sided agreement.

    Validation is a selection diagnostic, not an independent final goodness-of-
    fit test. A detected gradient is flagged; this routine never fits it away.
    """
    if any(m.get("error") is None or m["error"] <= 0 for m in measurements.values()):
        return False, "unconstrained"
    a, b = measurements["fit"], measurements["validation"]
    if abs(a["mean"]) < options.residual_detection_sigma * a["error"]:
        return False, "no significant training residual"
    if a["mean"] * b["mean"] <= 0 or abs(b["mean"]) < options.residual_validation_sigma * b["error"]:
        return False, "not independently confirmed"
    for left, right in [("fit", "validation"), ("lower", "upper")]:
        l, r = measurements[left], measurements[right]
        if abs(l["mean"] - r["mean"]) > 3 * (l["error"] + r["error"]):
            return False, "spatial inconsistency / possible gradient"
    if measurements["lower"]["mean"] * measurements["upper"]["mean"] <= 0:
        return False, "opposite-sided residuals"
    return True, "confirmed in disjoint validation strips"


def residual_correction(rate, variance, primary_error, wave, sy, valid, options,
                        dispersion=0.584, row_block=4):
    """Conservative second subtraction on the native, primary-subtracted image.

    Baseline continuum is order zero (a constant fitted over line-free sky).
    Broad 50-A bins test that simplification; inconsistent structure is flagged
    rather than removed with a flexible fit. Each airglow interval has a separate
    integrated detection gate and an unsmoothed native-dispersion mean profile.
    No correction is allowed to infer a cross-slit gradient automatically.
    """
    model, error = np.zeros_like(rate), np.zeros_like(rate)
    flags = np.zeros(rate.shape, "u2")
    if not options.residual_enabled:
        return model, error, flags, dict(enabled=False, applied=False, components=[])
    a = np.asarray(options.residual_fit_offsets) + options.source_radius
    b = np.asarray(options.residual_validation_offsets) + options.source_radius
    distance = abs(sy)
    fitmask = valid & (distance >= a[0]) & (distance < a[1])
    checkmask = valid & (distance >= b[0]) & (distance < b[1])
    # SCI validity never depends on ERR>0; valid is the unchanged geometry/DQ mask.
    masks = dict(fit=fitmask, validation=checkmask,
                 lower=fitmask & (sy < 0), upper=fitmask & (sy > 0))
    quiet = np.ones(rate.shape, bool)
    for lo, hi in AIRGLOW:
        quiet &= ~((wave >= lo) & (wave <= hi))
    lo, hi = np.quantile(wave[valid], [0.05, 0.95])
    quiet &= (wave >= lo) & (wave <= hi)
    records = []
    for label, region in [("continuum", quiet)] + [
        (f"airglow_{left:g}_{right:g}", (wave >= left) & (wave < right))
        for left, right in AIRGLOW
    ]:
        # Airglow is measured AFTER the accepted continuum correction, to avoid
        # counting a continuum offset twice. Its model uncertainty is included.
        residual = rate - model
        bound = primary_error + error
        measurements = {key: residual_measurement(residual, variance, bound, mask & region, row_block)
                        for key, mask in masks.items()}
        accepted, reason = residual_gate(measurements, options)
        record = dict(component=label, measurements=measurements, applied=False, reason=reason)
        if label == "continuum":
            bins = []
            for left in np.arange(np.floor(lo / 50) * 50, hi, 50):
                use = region & (wave >= left) & (wave < left + 50)
                item = {key: residual_measurement(residual, variance, bound, mask & use, row_block)
                        for key, mask in masks.items()}
                bins.append(dict(wavelength=[float(left), float(left + 50)], **item))
            record["bins"] = bins
            # A constant is only acceptable if broad independent bins agree.
            for row in bins:
                m = row["validation"]
                if (accepted and m.get("error") is not None and m["error"] > 0
                    and abs(m["mean"] - measurements["fit"]["mean"])
                        > 3 * (m["error"] + measurements["fit"]["error"])):
                    accepted, reason = False, "wavelength structure: constant model inadequate"
            candidate = np.full(rate.shape, measurements["fit"]["mean"] or 0.)
            candidate_error = np.full(rate.shape, measurements["fit"]["error"] or 0.)
        else:
            left, right = map(float, label.split("_")[1:])
            edges = np.r_[np.arange(left, right, dispersion), right]
            centers = (edges[:-1] + edges[1:]) / 2
            # No clipping, positivity floor or wavelength smoothing of airglow.
            specs = [spectra(residual, wave, sy, mask, edges, -np.inf, np.inf)
                     for mask in (fitmask, checkmask)]
            train, check = specs
            supported = (train["n"] >= 8) & (check["n"] >= 8)
            means = train["mean"]
            record["profile"] = dict(wavelength=centers.tolist(),
                fit=np.where(np.isfinite(means), means, 0).tolist(),
                validation=np.where(np.isfinite(check["mean"]), check["mean"], 0).tolist(),
                supported=supported.tolist())
            indices = np.flatnonzero(supported)
            # A fractional end bin can contain no detector centers. Permit one
            # such bin at each edge, but never bridge missing interior coverage.
            contiguous = (len(indices) >= 3 and indices[0] <= 1
                          and indices[-1] >= len(centers) - 2
                          and supported[indices[0]:indices[-1] + 1].all())
            if accepted and not contiguous:
                accepted, reason = False, "incomplete native-sampling airglow support"
            # Assess a resolved profile's consistency without forcing line noise
            # to zero. A failed shape test requires revisiting the primary model.
            se = np.hypot(train["pixel_se"], check["pixel_se"])
            positive = supported & np.isfinite(se) & (se > 0)
            shape_score = float(np.mean(((means[positive] - check["mean"][positive]) /
                                         se[positive])**2)) if positive.any() else None
            record["profile_discrepancy_score"] = shape_score
            if accepted and (shape_score is None or shape_score > 9):
                accepted, reason = False, "structured airglow discrepancy: revisit primary model"
            # Window edges are zero anchors, so there is no hard rectangular step.
            candidate = np.interp(wave, np.r_[left, centers[supported], right],
                                  np.r_[0., means[supported], 0.], left=0., right=0.)
            binned_bound = spectra(bound, wave, sy, fitmask, edges, -np.inf, np.inf)["mean"]
            binned_var = spectra(variance, wave, sy, fitmask, edges, -np.inf, np.inf)["mean"]
            se = np.fmax(train["pixel_se"], np.sqrt(np.maximum(binned_var, 0) /
                                                   np.maximum(train["n"], 1)))
            se = np.nan_to_num(se + binned_bound)
            candidate_error = np.interp(wave, np.r_[left, centers[supported], right],
                                        np.r_[0., se[supported], 0.], left=0., right=0.)
        record.update(applied=bool(accepted), reason=reason)
        if accepted:
            model[valid] += candidate[valid]
            error[valid] += candidate_error[valid]
        elif reason not in ("no significant training residual",):
            flags[valid & region] |= 1  # diagnostic only; does not change drizzle weights
        records.append(record)
    return model, error, flags, dict(enabled=True, applied=any(r["applied"] for r in records),
        fit_distances=a.tolist(), validation_distances=b.tolist(), components=records,
        continuum_order=0, gradient_fitted=False,
        uncertainty="Empirical clustered mean errors plus correlated primary-model bounds",
        limitation="Detection-gated estimates are conditional; very extended source flux in sky remains degenerate")


def add_residual(exposure, primary, wave, sy, valid, options, dispersion=0.584):
    """Append a validated residual model; preserve primary arrays separately."""
    if not primary.metadata.get("enabled", False):
        options = replace(options, residual_enabled=False)
    rate = exposure.counts / exposure.time - primary.rate
    variance = exposure.variance / exposure.time**2
    model, error, flags, metadata = residual_correction(
        rate, variance, primary.error, wave, sy, valid, options, dispersion)
    recovery = []
    if options.residual_enabled:
        for sigma in (20., 40., 80.):
            injected = 1e-4 * np.exp(-0.5 * (sy / sigma)**2)
            changed = residual_correction(rate + injected, variance, primary.error,
                         wave, sy, valid, options, dispersion)[0]
            roi = valid & (abs(sy) <= options.source_radius)
            recovery.append(dict(sigma_native_pixels=sigma, recovered_fraction=float(
                np.sum((injected - changed + model)[roi]) / np.sum(injected[roi]))))
    metadata["additional_stage_halo_recovery"] = recovery
    return BackgroundResult(primary.rate + model, primary.error + error, primary.flags,
        dict(primary.metadata, residual=metadata),
        dict(BKGPRIMARY=primary.rate, PRIMARYERR=primary.error,
             BKGRESID=model, RESERR=error, RESDQ=flags))


def estimate_rectified(science, variance, observed_wave, source_y, valid, options,
                       dispersion):
    """Fit one already calibrated x2d in its surface-brightness units.

    No exposure-time, sensitivity or pixel-area conversion is performed. The
    numerical normalization is reversed exactly. Spatial coordinates are in
    native-pixel equivalents so the same protected physical region is used.
    x2d pixels are correlated: block errors are diagnostic bounds, not a full
    reconstruction of the original detector covariance.
    """
    normalization = 1e-15
    values, var = science / normalization, variance / normalization**2
    if not options.enabled:
        return BackgroundResult(np.zeros_like(science), np.zeros_like(science),
                                np.zeros(science.shape, 'u2'), dict(enabled=False))
    edges = wavelength_edges(observed_wave[valid], options.continuum_bin, dispersion)
    width, scores = select_width(values, var, observed_wave, source_y, valid, options, edges)
    model, error, outside, anchors = fit(values, var, observed_wave, source_y, valid,
        inner=options.sky_inner, width=width, options=options, edges=edges)
    recovery = []
    amplitude = float(np.sqrt(np.median(var[valid])))
    amplitude = max(amplitude, np.finfo(float).eps)
    for sigma in (20., 40., 80.):
        injection = amplitude * np.exp(-.5 * (source_y / sigma)**2)
        changed = fit(values + injection, var, observed_wave, source_y, valid,
                      inner=options.sky_inner, width=width, options=options, edges=edges)[0]
        roi = valid & (abs(source_y) <= options.source_radius)
        recovery.append(dict(sigma_native_pixels=sigma, recovered_fraction=float(
            np.sum((injection - changed + model)[roi]) / np.sum(injection[roi]))))
    return BackgroundResult(model * normalization, error * normalization,
        np.where(valid, np.where(outside, 2, 0), 1).astype('u2'),
        dict(enabled=True, selected_width=width, scores=scores, anchors=anchors,
             halo_recovery=recovery, injection_peak_surface_brightness=amplitude*normalization,
             units="erg / (s cm2 Angstrom arcsec2)", stage="individual_MAST_x2d",
             wavelength_edges=edges.tolist(),
             uncertainty="Spatial blocks plus empirical/pipeline floors; pre-existing x2d covariance not fully known"))


def rectified_residual(science, base_error, primary_error, observed_wave, source_y,
                       valid, options, dispersion, spatial_step_native):
    """Assess ONLY the leftover background of an already primary-cleaned coadd.

    No image remapping occurs. The source mask remains in native-pixel units.
    Cluster height scales with output sampling so half-size pixels are not
    mistaken for twice as many independent spatial measurements.
    """
    options = options if options.enabled else replace(options, residual_enabled=False)
    normalization = 1e-15
    row_block = max(1, int(np.ceil(4 / spatial_step_native)))
    args = (base_error**2 / normalization**2, primary_error / normalization,
            observed_wave, source_y, valid, options, dispersion, row_block)
    model, error, flags, meta = residual_correction(science / normalization, *args)
    meta.update(units="erg / (s cm2 Angstrom arcsec2)",
                measurement_normalization=normalization, cluster_rows=row_block,
                stage="MAST_coadd_residual")
    # Measurement metadata remains normalized for stable diagnostics; record
    # this explicitly rather than presenting it as physical surface brightness.
    recovery = []
    amplitude = max(float(np.nanmedian(base_error[valid])), 1e-30) / normalization
    if options.residual_enabled:
        for sigma in (20., 40., 80.):
            injection = amplitude * np.exp(-.5 * (source_y / sigma)**2)
            changed = residual_correction(science / normalization + injection, *args)[0]
            roi = valid & (abs(source_y) <= options.source_radius)
            recovery.append(dict(sigma_native_pixels=sigma, recovered_fraction=float(
                np.sum((injection - changed + model)[roi]) / np.sum(injection[roi]))))
    meta['halo_recovery'] = recovery
    meta['injection_peak_surface_brightness'] = amplitude * normalization
    return BackgroundResult(model*normalization, error*normalization, flags, meta)
