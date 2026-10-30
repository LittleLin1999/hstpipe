"""Generic exposure and product QC; no target redshift or source-line windows."""

from pathlib import Path
import json
import numpy as np
from astropy.io import fits
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from .io import digest


def exposure_plot(path, e, result, valid):
    rate = e.counts / e.time
    limit = np.percentile(abs(rate[valid]), 99)
    norm = AsinhNorm(linear_width=max(limit / 30, 1e-8), vmin=-limit, vmax=limit)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for ax, data, label in zip(
        axes,
        [rate, result.rate, rate - result.rate],
        ["Original count rate", "Fitted 2D background", "Residual count rate"],
    ):
        ax.imshow(
            np.where(valid, data, np.nan),
            origin="lower",
            aspect="auto",
            norm=norm,
            cmap="RdBu_r",
            interpolation="none",
        )
        ax.set(title=label, xlabel="Detector X", ylabel="Detector Y")
    fig.suptitle(e.path.stem)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def validate_products(config, paths, grid, input_hashes, background_reports):
    """Checksums, matched photometry and independent source-adjacent residuals."""
    out = config.output / "diagnostics"
    out.mkdir(exist_ok=True)
    arrays = []
    for path in paths:
        with fits.open(path, checksum=True) as h:
            if not all(x.verify_checksum() == 1 and x.verify_datasum() == 1 for x in h):
                raise ValueError("FITS checksum failure")
            arrays.append(
                {
                    name: h[name].data.astype(float)
                    for name in ["SCI", "ERR", "EXPMAP", "NEXP"]
                }
            )
    with fits.open(paths[0]) as h:
        model_error = h["BKGERR"].data.astype(float)
        flags = h["BKGDQ"].data.copy()
    wave = grid.wave_start + np.arange(grid.nx) * grid.wave_step
    spatial = grid.spatial_start + np.arange(grid.ny) * grid.spatial_step
    base = np.load(config.output / "intermediate/flt/coadd/unsubtracted.npz")["SCI"]
    valid = (
        np.isfinite(base)
        & np.isfinite(arrays[1]["SCI"])
        & (abs(spatial[:, None]) < 0.5)
    )
    mast_baseline = config.output / "intermediate/x2d/unsubtracted.npz"
    mast_unclean = np.load(mast_baseline)["SCI"] if mast_baseline.exists() else arrays[1]["SCI"]
    # Photometric comparison includes sky in BOTH branches; background removal
    # is assessed separately. Wide-aperture sums reduce interpolation sensitivity.
    ratio = float(np.sum(base[valid]) / np.sum(mast_unclean[valid]))
    controls = []
    # Header scale is saved during prepare; no assumed source-line wavelength.
    stage = json.loads((config.output / "intermediate/flt/registration/inputs.json").read_text())
    native_scale = stage["native_spatial_scale"]
    radius = config.background.source_radius
    inner = config.background.sky_inner
    for lo in np.arange(np.ceil(wave.min() / 50) * 50, wave.max(), 50):
        for sign in (-1, 1):
            rows = (sign * spatial / native_scale >= radius + 4) & (
                sign * spatial / native_scale < inner - 2
            )
            mask = (
                rows[:, None]
                & (wave[None, :] >= lo)
                & (wave[None, :] < lo + 50)
                & np.isfinite(arrays[0]["SCI"])
                & ((flags & 3) == 0)
            )
            values = arrays[0]["SCI"][mask]
            if len(values) < 40:
                continue
            yy, xx = np.nonzero(mask)
            group = (yy // 16) * (grid.nx // 32 + 1) + xx // 32
            means = [float(values[group == k].mean()) for k in np.unique(group)]
            random = (
                float(np.std(means, ddof=1) / np.sqrt(len(means)))
                if len(means) > 3
                else np.inf
            )
            bound = float(np.mean(model_error[mask]))
            se = float(np.hypot(random, bound))
            controls.append(
                dict(
                    wavelength=[float(lo), float(lo + 50)],
                    side=int(sign),
                    mean=float(values.mean()),
                    uncertainty=se,
                    ratio=float(values.mean() / se),
                )
            )
    report = dict(
        input_hashes_preserved=all(
            digest(Path(p)) == sha for p, sha in input_hashes.items()
        ),
        products={p.name: digest(p) for p in paths},
        unsubtracted_joint_over_MAST=ratio,
        photometry_pass=abs(ratio - 1) < 0.03,
        independent_controls=controls,
        background=background_reports,
        limitations=[
            "Custom products; not official CALSTIS coadds.",
            "ERR omits full resampling covariance, including covariance already in MAST x2d.",
            "ERRBASE inherits low-count limitations of the pipeline ERR; zero errors are not grounds for pixel rejection.",
            "Background errors are correlated bounds, not independent pixel noise.",
            "Very broad source emission reaching sky strips can be subtracted.",
            "Header registration may require independently measured corrections for another target.",
        ],
    )
    if not report["input_hashes_preserved"]:
        raise ValueError("Input changed")
    (out / "report.json").write_text(json.dumps(report, indent=2))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    limit = np.nanpercentile(abs(base), 99.7)
    norm = AsinhNorm(linear_width=limit / 40, vmin=-limit, vmax=limit)
    extent = [wave[0], wave[-1], spatial[0], spatial[-1]]
    for ax, image, title in zip(
        axes,
        [arrays[1]["SCI"], base, arrays[0]["SCI"]],
        [
            "MAST x2d: per-exposure + residual sky",
            "Joint FLT drizzle before background",
            "Clean joint FLT drizzle",
        ],
    ):
        ax.imshow(
            image,
            origin="lower",
            aspect="auto",
            extent=extent,
            norm=norm,
            cmap="RdBu_r",
            interpolation="none",
        )
        ax.set(
            title=title,
            xlabel="Wavelength [A]",
            ylabel="Target-relative slit position [arcsec]",
        )
    fig.savefig(out / "comparison.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, name in zip(axes, ["EXPMAP", "NEXP"]):
        im = ax.imshow(arrays[0][name], origin="lower", aspect="auto", extent=extent)
        ax.set(title=name, xlabel="Wavelength [A]", ylabel="Slit position [arcsec]")
        fig.colorbar(im, ax=ax, label="seconds" if name == "EXPMAP" else "exposures")
    fig.savefig(out / "coverage.png", dpi=150)
    plt.close(fig)
    lines = [
        "# STIS reduction QC",
        "",
        f"Unsubtracted joint/MAST matched-aperture brightness ratio: {ratio:.7f}.",
        f'Photometric comparison within 3%: {report["photometry_pass"]}.',
        "",
        "| Input | Selected background strip width |",
        "| --- | ---: |",
    ]
    for name, r in background_reports.items():
        lines.append(f'| {name} | {r.get("selected_width",0):g} |')
    if controls:
        worst = max(controls, key=lambda r: abs(r["ratio"]))
        lines += [
            "",
            f'Maximum absolute independent-control residual/uncertainty ratio: {abs(worst["ratio"]):.3f}',
            f'(side {worst["side"]:+d}, {worst["wavelength"][0]:g}–{worst["wavelength"][1]:g} A).',
            "These are conservative screening ratios, not calibrated Gaussian detection significances.",
            "With native residual checks enabled, these guard regions participate in acceptance diagnostics and are not an independent final validation sample.",
        ]
    recoveries = {}
    for r in background_reports.values():
        for injection in r.get("halo_recovery", []):
            recoveries.setdefault(injection["sigma_native_pixels"], []).append(
                injection["recovered_fraction"]
            )
    if recoveries:
        lines += [
            "",
            "Fixed-model-width halo injection recovery inside the protected source region:",
            "",
            "| Gaussian sigma (native pixels) | Range across exposures |",
            "| ---: | ---: |",
        ]
        for sigma, values in sorted(recoveries.items()):
            lines.append(f"| {sigma:g} | {min(values):.4%}–{max(values):.4%} |")
        lines += [
            "",
            "This tests the selected background operator, not a new PSF model or complete injection/selection rerun.",
        ]
    lines += [
        "",
        "Final products:",
        *[f"- {p.name}" for p in paths],
        "",
        "Limitations:",
        *[f"- {s}" for s in report["limitations"]],
    ]
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")
    return report


def extraction_diagnostics(config, image_spectrum, individual_coadd,
                           selected_flt, selected_drz, residual_records):
    """Plot the selected final-product convention and a separate calibration check.

    ``FLUX_APER`` is always a diffuse-convention auxiliary column. It never
    changes meaning with ``flux_mode`` and is not the default final-product
    flux when ``flux_mode='point'``.
    """
    directory = config.output / "intermediate/x1d"
    wave = image_spectrum["WAVELENGTH"]
    good = (image_spectrum["DQ"] == 0) & (individual_coadd["DQ"] == 0)
    a, b = image_spectrum["FLUX"], individual_coadd["FLUX_APER"]
    bins = []
    for low in np.arange(np.ceil(wave.min() / 50) * 50, wave.max(), 50):
        use = good & (wave >= low) & (wave < low + 50)
        if use.sum() < 10:
            continue
        mean_a, mean_b = float(a[use].mean()), float(b[use].mean())
        bins.append(
            dict(
                wavelength=[float(low), float(low + 50)],
                pixels=int(use.sum()),
                from2d=mean_a,
                from_individual=mean_b,
                difference=mean_a - mean_b,
                ratio=mean_a / mean_b if mean_b != 0 else None,
            )
        )
    if not good.any():
        raise ValueError("No common fully covered extraction bins")
    ratio = float(np.sum(a[good]) / np.sum(b[good]))
    selected_good = (
        (selected_flt["DQ"] == 0) & (selected_drz["DQ"] == 0)
        & np.isfinite(selected_flt["FLUX_BOX"])
        & np.isfinite(selected_drz["FLUX_BOX"])
    )
    if not selected_good.any():
        raise ValueError("No common bins in the selected final 1D convention")
    selected_ratio = float(
        np.sum(selected_drz["FLUX_BOX"][selected_good])
        / np.sum(selected_flt["FLUX_BOX"][selected_good])
    )
    fig, axes = plt.subplots(
        2, 1, figsize=(12, 8), sharex=True, constrained_layout=True
    )
    mode = config.extraction.flux_mode.upper()
    rejected = [r["input"] for r in residual_records if not r["residual_applied"]]
    axes[0].plot(
        wave, selected_flt["FLUX_BOX"],
        label=f"FLT boxcar: selected {mode} calibration", lw=0.7,
    )
    axes[0].plot(
        wave, selected_drz["FLUX_BOX"],
        label=f"DRZ boxcar: selected {mode} calibration", lw=0.7, alpha=0.7,
    )
    axes[0].set_title("Final x1d product: FLT and DRZ in the same configured flux mode")
    axes[1].plot(
        wave, individual_coadd["FLUX"],
        label="FLT point-source response (shared native sky)", lw=0.7,
    )
    axes[1].plot(
        wave, b,
        label="FLT diffuse-calibration comparison (FLUX_APER; auxiliary)",
        lw=0.7, alpha=0.7,
    )
    axes[1].set_title("Calibration conventions on the same extracted FLT NET signal")
    for ax in axes:
        ax.set_ylabel("Flux density [erg s$^{-1}$ cm$^{-2}$ A$^{-1}$]")
        ax.legend(fontsize=8)
        ax.axhline(0, color="black", lw=0.5)
    axes[1].set_xlabel("Heliocentric wavelength [A]")
    fig.savefig(directory / "comparison.png", dpi=150)
    plt.close(fig)
    report = dict(
        matched_aperture_broadband_ratio=ratio,
        matched_aperture_agreement_within_10pct=bool(abs(ratio - 1) <= 0.10),
        selected_flux_mode=mode,
        selected_mode_drz_over_flt_broadband_ratio=selected_ratio,
        residual_background_not_applied_exposures=rejected,
        spectral_bin_comparison=bins,
        limitations=[
            "FLUX from 2D is finite-aperture diffuse flux; FLT FLUX uses the official point-source response.",
            "Compare 2D diffuse FLUX with individual-coadd FLUX_APER, not point-calibrated FLUX.",
            "CALSTIS BACKCORR=OMIT: any additional sky is already in the shared native model.",
            "Agreement of two branches sharing the same Stage-1 sky model is not proof that model is unbiased.",
            "x2d-native calibration and native/rectified aperture shapes differ; no forced flux matching.",
            "2D ERROR is a conservative covariance bound; ERR_DIAG ignores drizzle covariance.",
            "Native x1d ERR_STAT uses squared edge weights, with original CALSTIS ERROR saved separately per exposure.",
            "Low-count input ERR limitations remain; subtraction does not remove photon noise.",
            "Background bounds are correlated across wavelength; they cannot be quadrature-summed as independent bins.",
            "The native-coadd COVAR1_STAT_APER describes adjacent bins only, not the full covariance.",
        ],
    )
    lines = [
        "# 1D extraction QC",
        "",
        f"Aperture full height: {config.extraction.height_native_pixels:g} native pixels; center: {config.extraction.center_arcsec:g} arcsec.",
        "Official extraction: unweighted, maxsrch=0, BACKCORR=OMIT; original count ERR and DQ retained.",
        f"Shared native residual correction applied in {len(residual_records)-len(rejected)}/{len(residual_records)} exposures.",
        "Unclipped residual means require two-sided agreement and disjoint validation; no repeated x1d subtraction.",
        "",
        f"Matched finite-aperture broadband ratio (2D / individual x1d): {ratio:.7f}.",
        f"Selected {mode} mode DRZ/FLT broadband ratio: {selected_ratio:.7f}.",
        "The top plot panel uses the selected final-product calibration on both branches.",
        "The bottom panel isolates the point-versus-diffuse calibration difference on the selected FLT NET spectrum.",
        "This ratio is a comparison, not a calibration adjustment or a narrow-line validation.",
        ("QC: matched-aperture comparison differs by more than 10%; inspect the "
         "residual x1d background and spatial aperture before treating branches as interchangeable."
         if abs(ratio - 1) > 0.10 else
         "QC: matched-aperture broadband comparison agrees within 10%."),
        "",
        "| Wavelength [A] | 2D / individual aperture flux |",
        "| --- | ---: |",
    ]
    for row in bins:
        lines.append(
            f'| {row["wavelength"][0]:g}–{row["wavelength"][1]:g} | {row["ratio"]:.5f} |'
            if row["ratio"] is not None
            else "| zero denominator | undefined |"
        )
    lines += ["", "Limitations:", *["- " + x for x in report["limitations"]]]
    (directory / "REPORT.md").write_text("\n".join(lines) + "\n")
    return report


def residual_plot(path, exposure, result, unit_label="Residual [count/s/pixel]"):
    """Display decision measurements, including cases where nothing is removed."""
    report = result.metadata.get("residual", {})
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for ax, record in zip(axes.ravel(), report.get("components", [])):
        if record["component"] == "continuum":
            for key in ("fit", "validation", "lower", "upper"):
                bins = [b for b in record.get("bins", []) if b[key].get("error") is not None]
                ax.errorbar([np.mean(b["wavelength"]) for b in bins],
                            [b[key]["mean"] for b in bins],
                            yerr=[b[key]["error"] for b in bins], fmt=".-", ms=3,
                            lw=0.7, label=key)
        else:
            profile = record.get("profile", {})
            for key in ("fit", "validation"):
                if profile:
                    ax.plot(profile["wavelength"], profile[key], lw=.7, label=key)
        ax.axhline(0, color="black", lw=.6)
        ax.set(title=f'{record["component"]}: applied={record["applied"]}\n{record["reason"]}',
               xlabel="Observed wavelength [A]", ylabel=unit_label)
        ax.legend(fontsize=8)
    fig.suptitle(exposure.path.stem + " — additional residual check")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def residual_summary(config, path, grid, reports):
    """Compare primary-only and additional-residual coadds on identical grids."""
    primary = np.load(config.output / "intermediate/flt/coadd/primary_only.npz")["SCI"]
    with fits.open(path) as h:
        final = h["SCI"].data.astype(float)
        usable = (h["DQ"].data == 0) & np.isfinite(primary) & np.isfinite(final)
    wave = grid.wave_start + np.arange(grid.nx) * grid.wave_step
    spatial = grid.spatial_start + np.arange(grid.ny) * grid.spatial_step
    meta = json.loads((config.output / "intermediate/flt/registration/inputs.json").read_text())
    scale = meta["native_spatial_scale"]
    delta = primary - final
    regions = dict(aperture=abs(spatial / scale) <= config.extraction.height_native_pixels / 2,
                   wings=(abs(spatial / scale) >= 20) & (abs(spatial / scale) <= 80))
    rows = []
    for lo in np.arange(np.ceil(wave.min() / 50) * 50, wave.max(), 50):
        band = (wave >= lo) & (wave < lo + 50)
        for region, spatial_mask in regions.items():
            mask = usable & spatial_mask[:, None] & band[None, :]
            if not mask.any():
                continue
            before, after = float(primary[mask].sum()), float(final[mask].sum())
            rows.append(dict(wavelength=[float(lo), float(lo + 50)], region=region,
                primary_sum=before, final_sum=after, removed_sum=before-after,
                fractional_change=(after-before)/before if before != 0 else None))
    extent = [wave[0], wave[-1], spatial[0], spatial[-1]]
    limit = float(np.nanpercentile(abs(primary[usable]), 99.7))
    norm = AsinhNorm(linear_width=max(limit/40, 1e-30), vmin=-limit, vmax=limit)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    for ax, image, title in zip(axes, (primary, final, delta),
        ("Primary subtraction only", "Primary + validated residual", "Additional removed background")):
        ax.imshow(image, origin="lower", aspect="auto", extent=extent, norm=norm,
                  cmap="RdBu_r", interpolation="none")
        ax.set(title=title, xlabel="Wavelength [A]", ylabel="Slit position [arcsec]")
    out = config.output / "diagnostics"
    fig.savefig(out / "residual_comparison.png", dpi=150)
    plt.close(fig)
    bands = [(1205,1225), (1295,1315), (1345,1365), (1400,1450), (1450,1500), (1550,1600)]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for ax, (lo, hi) in zip(axes.ravel(), bands):
        mask = usable & (wave[None, :] >= lo) & (wave[None, :] < hi)
        for image, label in ((primary, "primary only"), (final, "with residual")):
            summed = np.where(mask, image, 0).sum(axis=1)
            n = mask.sum(axis=1)
            profile = np.divide(summed, n, out=np.full(len(n), np.nan), where=n > 0)
            ax.plot(spatial / scale, profile, lw=.8, label=label)
        ax.set(xlim=(-100,100), title=f"{lo}–{hi} A", xlabel="Native pixels from target",
               ylabel="Mean surface brightness")
        ax.axhline(0, color="black", lw=.5)
        ax.legend(fontsize=8)
    fig.savefig(out / "residual_profiles.png", dpi=150)
    plt.close(fig)
    decisions = {name: report.get("residual", {}) for name, report in reports.items()}
    summary = dict(exposures=decisions, flux_comparison=rows,
        max_abs_science_change=float(np.max(abs(delta[usable]))),
        note="Coadd comparisons use fixed calibration/mapping/weights. Float32 output rounding may remain.")
    (out / "residual_report.json").write_text(json.dumps(summary, indent=2))
    lines = ["# Additional residual-background assessment", "",
        "All estimates are signed count/s/native pixel. Errors include correlated primary-model bounds.",
        "The ratios below are conservative screening statistics, not Gaussian probabilities.", "",
        "| Exposure | Component | Fit mean ± bound | Fit ratio | Validation ratio | Applied / reason |",
        "| --- | --- | --- | ---: | ---: | --- |"]
    for name, report in decisions.items():
        for record in report.get("components", []):
            a, b = record["measurements"]["fit"], record["measurements"]["validation"]
            mean = f'{a["mean"]:.4g} ± {a["error"]:.3g}' if a.get("error") else "unconstrained"
            ratio_a = f'{a["significance"]:.2f}' if a.get("significance") is not None else "—"
            ratio_b = f'{b["significance"]:.2f}' if b.get("significance") is not None else "—"
            lines.append(f'| {name} | {record["component"]} | {mean} | {ratio_a} | {ratio_b} | {record["applied"]}: {record["reason"]} |')
    lines += ["", "## Additional-stage injection recovery", "",
        "Primary model held fixed; residual fitting AND decision gates are rerun. This is not an end-to-end primary-model test.", "",
        "| Exposure | Injected Gaussian sigma [native pixels] | Recovered flux inside ±80 |",
        "| --- | ---: | ---: |"]
    for name, report in decisions.items():
        for recovery in report.get("additional_stage_halo_recovery", []):
            lines.append(f'| {name} | {recovery["sigma_native_pixels"]:g} | {recovery["recovered_fraction"]:.3%} |')
    lines += ["", "## Limitations", "",
        "- Sky strips closer than the primary fit can contain extended emission. Inspect the masks and injection results.",
        "- A constant continuum is the only automatic smooth residual model. Structured continuum/gradients require investigation, not higher-order flattening.",
        "- Airglow is corrected only when independently confirmed with adequate native-wavelength coverage.",
        "- Validation strips are used for model acceptance; they are not a separate unbiased final test sample.",
        "- Background uncertainty is a correlated bound, not independent photon noise. Original low-count ERR limitations remain.",
        "- Failure to detect a correction does not establish absence of a residual background.",
        "- CALSTIS BACKCORR=OMIT prevents repeating this correction in the FLT branch."]
    (out / "RESIDUAL_REPORT.md").write_text("\n".join(lines) + "\n")
    return summary


def rectified_background_plot(path, science, model, valid, wave, spatial, label):
    """Calibrated image, model and residual with identical display scaling."""
    limit = max(float(np.nanpercentile(abs(science[valid]), 99)), 1e-30)
    norm = AsinhNorm(linear_width=limit/30, vmin=-limit, vmax=limit)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    extent = [wave[0], wave[-1], spatial[0], spatial[-1]]
    for ax, data, title in zip(axes, (science, model, science-model),
                                ('Input', 'Subtracted model', 'After subtraction')):
        im = ax.imshow(np.where(valid, data, np.nan), origin='lower', aspect='auto',
                       extent=extent, cmap='RdBu_r', norm=norm, interpolation='none')
        ax.set(title=title, xlabel='Wavelength [A]', ylabel='Slit position [arcsec]')
    fig.colorbar(im, ax=axes, label='erg s$^{-1}$ cm$^{-2}$ A$^{-1}$ arcsec$^{-2}$')
    fig.suptitle(label)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def mast_background_report(config, grid, unclean, before, after, reports):
    """Document the two distinct x2d background operations and their impact."""
    out = config.output/'diagnostics'
    out.mkdir(exist_ok=True)
    wave = grid.wave_start + np.arange(grid.nx)*grid.wave_step
    spatial = grid.spatial_start + np.arange(grid.ny)*grid.spatial_step
    valid = (after['DQ'] == 0) & np.isfinite(after['SCI'])
    rectified_background_plot(out/'mast_individual_background.png', unclean,
                              before['BKG'], valid, wave, spatial, 'MAST x2d: individual models mapped to coadd')
    rectified_background_plot(out/'mast_coadd_residual.png', before['SCI'],
                              after['BKGRESID'], valid, wave, spatial, 'MAST x2d: additional coadd residual only')
    meta = reports['coadd_residual']
    # Reuse the shared decision plot, with the correct physical-unit label.
    from types import SimpleNamespace
    residual_plot(out/'mast_residual_decisions.png', SimpleNamespace(path=Path('MAST_coadd')),
                  SimpleNamespace(metadata={'residual':meta}),
                  unit_label='Residual surface brightness / 1e-15')
    native_scale = json.loads((config.output/'intermediate/flt/registration/inputs.json').read_text())['native_spatial_scale']
    inner = abs(spatial) <= config.extraction.height_native_pixels*native_scale/2
    bands = []
    for lo in np.arange(np.ceil(wave.min()/50)*50, wave.max(), 50):
        use = valid & inner[:,None] & (wave[None,:] >= lo) & (wave[None,:] < lo+50)
        if use.any():
            bands.append(dict(wavelength=[float(lo),float(lo+50)],
                original_sum=float(unclean[use].sum()), primary_sum=float(before['SCI'][use].sum()),
                final_sum=float(after['SCI'][use].sum())))
    report = dict(background=reports, aperture_sums=bands,
        max_coadd_residual=float(np.max(abs(after['BKGRESID'][valid]))),
        units='erg / (s cm2 Angstrom arcsec2)',
        covariance='Pre-existing x2d and new coadd covariance not fully known; errors are bounds.')
    (out/'mast_background_report.json').write_text(json.dumps(report, indent=2))
    lines = ['# MAST x2d background QC', '',
        'Each original archive x2d has an independently fitted surface-brightness background removed before combination.',
        'The second model uses the primary-cleaned coadd and source-free guard regions unused by the primary fits.',
        'SCI = original coadd - mapped individual backgrounds - accepted coadd residual.',
        'Original calibration/WCS are retained; neither exposure time nor pixel area is divided out again.', '',
        '| Input | Selected strip width [native pixels] |', '| --- | ---: |']
    for name, r in reports.items():
        if name != 'coadd_residual': lines.append(f'| {name} | {r.get("selected_width",0):g} |')
    lines += ['', '| Coadd component | Applied | Reason |', '| --- | --- | --- |']
    for c in meta['components']:
        lines.append(f'| {c["component"]} | {c["applied"]} | {c["reason"]} |')
    lines += ['', 'Errors retain original photon noise. Background errors are separate correlated bounds.',
        'The rectified-pixel covariance is not reconstructed; nominal significance is a screening statistic.',
        'Extended source flux reaching the sky strips is degenerate with background; inspect halo-recovery metadata.',
        f'Subtraction linearity relative discrepancy: {meta["subtraction_linearity_relative_error"]:.3g}.']
    (out/'MAST_BACKGROUND_REPORT.md').write_text('\n'.join(lines)+'\n')


def extraction_fit_plot(directory, label, science, stat_error, bkg_error, valid,
                        aperture, profile, optimal_flux, wave, spatial,
                        profile_metadata, units, correlated=False):
    """Save the actual extraction model and diagnostic residuals, not a refit.

    Pattern follows unfold_jwst's data/model/chi/spectrum display and JWST's
    save_profile/save_scene_model/save_residual_image outputs. Our equal-noise
    estimator and signed MAMA pixels are unchanged. Diagonal chi2 is diagnostic,
    not a likelihood test of this estimated profile on correlated low-count data.
    """
    from .optimal import profile_fit_diagnostics, diagnostic_noise_scale
    from astropy.table import Table

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    area = np.broadcast_to(aperture, science.shape)
    diag = profile_fit_diagnostics(science, stat_error, valid, area, profile, optimal_flux)
    spatial = np.asarray(spatial)
    diagnostic_error, noise_floor = diagnostic_noise_scale(science, stat_error, valid, spatial)
    scaled = profile_fit_diagnostics(science, diagnostic_error, valid, area, profile, optimal_flux)
    crop = np.flatnonzero((abs(spatial) <= max(45,profile_metadata.get('fitting_height_native',60)/2+4)) | np.any(area > 0, axis=1))
    section = slice(crop.min(), crop.max()+1)
    yy = spatial[section]
    extent = [wave[0],wave[-1],yy[0],yy[-1]]
    active = valid & (area > 0)
    columns = diag['ACCEPTED_COLUMNS'] & np.isfinite(optimal_flux)
    diag_mask = diag['DIAGMASK'] & columns[None,:]
    # The independently trained wider template is shown beyond the amplitude
    # fitting aperture. Residual statistics use only the extraction footprint.
    data = np.where(valid, science, np.nan)[section]
    limit = max(float(np.nanpercentile(abs(data),99)), 1e-30)
    norm = AsinhNorm(linear_width=limit/30, vmin=-limit,vmax=limit)
    fig, axes = plt.subplots(5,1,figsize=(14,12),sharex=True,constrained_layout=True)
    panels = [('Data',data),('Full template prediction; amplitude fitted in the wider optimal aperture',diag['MODEL_FULL'][section]),
              ('Data − model (within fitted aperture)',diag['RESIDUAL'][section]),
              ('Raw ERRBASE standardized residual; near-zero errors can dominate',diag['CHI_DIAG'][section]),
              ('Off-source noise-floor display; not a Gaussian significance map',scaled['CHI_DIAG'][section])]
    for i,(ax,(title,image)) in enumerate(zip(axes,panels)):
        kw = dict(vmin=-5,vmax=5) if i>=3 else dict(norm=norm)
        im = ax.imshow(image,origin='lower',aspect='auto',extent=extent,cmap='RdBu_r',interpolation='none',**kw)
        ax.set(title=title,ylabel='Native pixels from target')
        fig.colorbar(im,ax=ax,label='Standardized residual' if i>=3 else units)
        # Actual fractional aperture footprint, including native trace variations.
        lower=np.full(len(wave),np.nan);upper=lower.copy()
        for j in range(len(wave)):
            rows=np.flatnonzero(area[:,j]>0)
            if len(rows): lower[j],upper[j]=spatial[rows[0]],spatial[rows[-1]]
        ax.plot(wave,lower,'k--',lw=.5);ax.plot(wave,upper,'k--',lw=.5)
        for edge in (-1,1):
            ax.axhline(edge*profile_metadata.get('boxcar_height_native',11)/2,color='tab:blue',ls=':',lw=.8)
    axes[-1].set_xlabel('Heliocentric wavelength [A]')
    fig.suptitle(label+' — empirical equal-noise profile extraction')
    fig.savefig(directory/(label+'_model.png'),dpi=140);plt.close(fig)

    quiet = np.isfinite(wave) & columns
    for lo,hi in ((1210,1223),(1298,1310),(1348,1364)):
        quiet &= ~((wave>=lo)&(wave<=hi))
    peak = None
    if quiet.any():
        candidates=np.flatnonzero(quiet)
        peak=float(wave[candidates[np.nanargmax(optimal_flux[candidates])]])
    bands=[('Template wavelengths',quiet)]
    if peak is not None: bands.append((f'Brightest non-airglow feature ±5 A ({peak:.1f} A)',columns&(abs(wave-peak)<=5)))
    fig,axes=plt.subplots(2,2,figsize=(13,8),constrained_layout=True)
    for ax,(title,band) in zip(axes[0],bands):
        for image,name in ((science,'Signed data'),(diag['MODEL_FULL'],'Full profile × fitted amplitude')):
            use=valid&band[None,:]&np.isfinite(image)
            n=use.sum(axis=1)
            mean=np.divide(np.where(use,image,0).sum(axis=1),n,out=np.full(len(n),np.nan),where=n>0)
            ax.plot(spatial,mean,lw=1,label=name)
        span=max(45,profile_metadata.get('fitting_height_native',151)/2+4)
        ax.set(xlim=(-span,span),title=title,xlabel='Native pixels from target',ylabel=units)
        ax.axhline(0,color='black',lw=.5)
        for sign in (-1,1):
            ax.axvline(sign*profile_metadata.get('boxcar_height_native',11)/2,color='tab:blue',ls=':',lw=.7)
            ax.axvline(sign*profile_metadata.get('optimal_height_native',140)/2,color='black',ls='--',lw=.7)
        ax.legend(fontsize=8)
    axes[1,0].plot(wave,np.where(columns,diag['CHI2_PER_DOF'],np.nan),lw=.7,label='Raw ERRBASE')
    axes[1,0].plot(wave,np.where(columns,scaled['CHI2_PER_DOF'],np.nan),lw=.7,label='Off-source noise floor')
    axes[1,0].legend(fontsize=8)
    axes[1,0].set(yscale='log',xlabel='Wavelength [A]',ylabel='χ² diagonal / nominal dof',
                  title='Descriptive only: covariance and learned profile omitted')
    values=diag['CHI_DIAG'][diag_mask]
    axes[1,1].hist(values,bins=np.linspace(-8,8,81),density=True,histtype='step',label='Raw ERRBASE')
    scaled_values=scaled['CHI_DIAG'][scaled['DIAGMASK']&columns[None,:]]
    axes[1,1].hist(scaled_values,bins=np.linspace(-8,8,81),density=True,histtype='step',label='Off-source floor')
    x=np.linspace(-8,8,200)
    axes[1,1].plot(x,np.exp(-x*x/2)/np.sqrt(2*np.pi),'k--',lw=.7,label='N(0,1), reference only')
    axes[1,1].set(xlabel='Diagonal standardized residual',ylabel='Density',
                  title=f'{np.mean(abs(values)>8):.1%} outside plotted range' if len(values) else 'No usable diagnostic errors')
    axes[1,1].legend(fontsize=8)
    fig.suptitle(label+f' — collapsed profiles; flat fallback={profile_metadata["fallback_to_flat"]}')
    fig.savefig(directory/(label+'_profile.png'),dpi=140);plt.close(fig)

    h=fits.Header({'OFFICIAL':False,'OPTMETH':'EMPIRICAL_EQUALVAR','DIAGONLY':True,
                   'COVKNOWN':False,'ROWSTART':int(crop.min()),'FLATFALL':profile_metadata['fallback_to_flat']})
    h.add_history('MODEL is the exact profile_extract forward model, including fractional aperture pixels.')
    h.add_history('Chi2 uses ERRBASE diagonal only, one nominal fitted amplitude per column.')
    h.add_history('No p-values: template learned from these data; low-count errors and covariance uncertain.')
    hd=[fits.PrimaryHDU(header=h)]
    for name,array in dict(DATA=science,ERRBASE=stat_error,BKGERR=bkg_error,
         MODEL=diag['MODEL'],RESIDUAL=diag['RESIDUAL'],CHI_DIAG=diag['CHI_DIAG'],
         PROFILE=diag['PROFILE'],PROFILE_FULL=diag['PROFILE_FULL'],MODEL_FULL=diag['MODEL_FULL'],APERTURE=area,VALID=valid,DIAGMASK=diag['DIAGMASK'],DIAGERR=diagnostic_error,CHI_FLOOR=scaled['CHI_DIAG']).items():
        hh=fits.Header({'BUNIT':units if name in ('DATA','ERRBASE','BKGERR','MODEL','MODEL_FULL','RESIDUAL','DIAGERR') else ''})
        hd.append(fits.ImageHDU(np.asarray(array[section],dtype='f4'),header=hh,name=name))
    hd.append(fits.ImageHDU(yy,name='SPATIAL'))
    hd.append(fits.BinTableHDU(Table(dict(WAVELENGTH=wave,FLUX_OPT=optimal_flux,
         CHI2_DIAG=diag['CHI2_DIAG'],NDOF_NOMINAL=diag['NDOF_NOMINAL'],
         CHI2_PER_DOF=diag['CHI2_PER_DOF'],CHI2_FLOOR=scaled['CHI2_PER_DOF'],NOISE_FLOOR=noise_floor,FULL_COVERAGE=diag['FULL_COLUMNS'],ACCEPTED=columns)),name='DIAGNOSTICS'))
    fits.HDUList(hd).writeto(directory/(label+'_fit.fits'),overwrite=False,checksum=True)
    summary=dict(profile_metadata,diagonal_chi2_only=True,correlated_input=bool(correlated),
        unusable_error_pixels=diag['ERROR_UNUSABLE_PIXELS'],complete_columns=int(np.sum(columns & diag['FULL_COLUMNS'])),
        accepted_columns=int(columns.sum()),
        accepted_masked_columns=int(np.sum(columns & ~diag['FULL_COLUMNS'])),
        median_nominal_reduced_chi2=float(np.nanmedian(diag['CHI2_PER_DOF'][columns])),
        fitted_feature_wavelength=peak,science_unit=units,
        profile_centroid_native_pixels=float(np.sum(profile*spatial)),
        median_noise_floor_chi2=float(np.nanmedian(scaled['CHI2_PER_DOF'][columns])),
        diagnostic_noise_floor='Unclipped off-source variance, |y|=100..150, 32-column blocks; DISPLAY ONLY',
        limitation='No chi2 p-value or formal fit acceptance: equal-noise fitting, learned profile, low-count errors and covariance.')
    (directory/(label+'_fit.json')).write_text(json.dumps(summary,indent=2))
    return summary


def extraction_branch_plot(config, branches):
    """Compare FLT/DRZ/X2D in the same final calibration mode, with uncertainties."""
    out=config.output/'intermediate/x1d'
    fig,axes=plt.subplots(3,1,figsize=(14,10),sharex=True,constrained_layout=True)
    for ax,(name,data) in zip(axes,branches.items()):
        wave=data['WAVELENGTH']
        for method,color in (('BOX','tab:blue'),('OPT','tab:orange')):
            good=data['DQ_'+method]==0
            f=np.where(good,data['FLUX_'+method],np.nan)
            e=np.where(good,data['ERR_'+method],np.nan)
            ax.plot(wave,f,color=color,lw=.7,label='Boxcar' if method=='BOX' else 'Empirical profile (equal-noise)')
            ax.fill_between(wave,f-e,f+e,color=color,alpha=.15,lw=0)
        ax.set(title=name+' — '+config.extraction.flux_mode.upper(),ylabel='Flux density [erg s⁻¹ cm⁻² A⁻¹]')
        ax.axhline(0,color='black',lw=.5)
        ax.legend(fontsize=8)
    axes[-1].set_xlabel('Heliocentric wavelength [A]')
    fig.savefig(out/'extraction_methods.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(14,8),sharex=True,constrained_layout=True)
    for ax,method in zip(axes,('BOX','OPT')):
        for name,data in branches.items():
            ax.plot(data['WAVELENGTH'],np.where(data['DQ_'+method]==0,data['FLUX_'+method],np.nan),lw=.7,label=name)
        ax.set(title=method+' extraction; all branches background-treated',ylabel='Flux density')
        ax.legend();ax.axhline(0,color='black',lw=.5)
    axes[-1].set_xlabel('Heliocentric wavelength [A]')
    fig.savefig(out/'all_branches.png',dpi=150);plt.close(fig)
    # Summarize actual saved fits; this is part of the repeatable extraction stage.
    fitsummary={p.stem:json.loads(p.read_text()) for p in sorted((out/'diagnostics').glob('*_fit.json'))}
    lines=['# Extraction diagnostics', '',
        'All three final branches use the selected '+config.extraction.flux_mode+' calibration.',
        'Boxcar and empirical-profile spectra have separate error arrays. The latter uses equal-noise weights, not inverse-variance MAMA weights.',
        'Profiles are learned only within calibrated wavelength support. Fallback status and extraction-center-relative centroids are recorded.', '',
        '| Dataset | Profile centroid [native pixels] | Raw diagonal chi2/dof median | Noise-floor display median | Flat fallback |',
        '| --- | ---: | ---: | ---: | --- |']
    for name,r in fitsummary.items():
        lines.append(f'| {name} | {r["profile_centroid_native_pixels"]:.3f} | {r["median_nominal_reduced_chi2"]:.4g} | {r["median_noise_floor_chi2"]:.4g} | {r["fallback_to_flat"]} |')
    ratios={}
    ref=branches['FLT']
    for name,data in branches.items():
        ratios[name]={}
        for method in ('BOX','OPT'):
            good=(data['DQ_'+method]==0)&(ref['DQ_'+method]==0)
            ratios[name][method]=float(data['FLUX_'+method][good].sum()/ref['FLUX_'+method][good].sum())
    lines += ['', '## Same-mode broadband comparison', '',
        '| Branch / FLT | Boxcar | Empirical profile |','| --- | ---: | ---: |']
    for name,r in ratios.items():lines.append(f'| {name} | {r["BOX"]:.6f} | {r["OPT"]:.6f} |')
    lines += ['', '## Reading these diagnostics', '',
        '- `all_branches.png`: FLT, DRZ, and background-treated X2D spectra.',
        '- `extraction_methods.png`: boxcar/profile extraction with separate uncertainty bands.',
        '- `diagnostics/*_model.png`: data, actual model, residual, raw chi and off-source-noise-scaled display.',
        '- `diagnostics/*_profile.png`: collapsed profile, brightest calibrated non-airglow feature, diagonal chi2 and histogram.',
        '- `diagnostics/*_fit.fits`: numerical model, residuals, profiles, masks and per-wavelength statistics.',
        '', 'Large raw chi2 is expected when pipeline errors approach zero; it cannot certify a poor or good spatial model by itself.',
        'The alternate off-source noise floor is a display aid and does not change science uncertainties, weights, or masks.',
        'Optimal normalization is fixed before masking; the fitted flux includes profile-inferred missing pixels. No spectral-gap interpolation is used.',
        'Use DQ_BOX and DQ_OPT separately. DQ is a compatibility alias for DQ_BOX.',
        'Neither statistic has formal p-values: the learned profile, non-Gaussian low counts and resampling covariance violate that interpretation.',
        'Full templates are fitted over 151 native pixels by default; optimal amplitudes use 140 pixels and boxcar uses 11. Blue dotted boundaries mark boxcar; black dashed boundaries mark optimal.',
        'A profile defined over a finite region is not an infinite-aperture instrumental PSF; unmeasured wings remain uncertain.',
        'The continuum-based center is recorded in centering.json and used by all apertures. Image registration is unchanged.',
        '', 'Diagnostic layout references: local unfold_jwst/extract1d.py (data, model, chi and boxcar/optimal panels),',
        'and https://jwst.readthedocs.io/en/latest/jwst/extract_1d/arguments.html (save_profile, save_scene_model, save_residual_image).']
    (out/'EXTRACTION_DIAGNOSTICS.md').write_text('\n'.join(lines)+'\n')
    extraction_coverage_plot(config,branches)
    (out/'diagnostic_summary.json').write_text(json.dumps(dict(fits=fitsummary,branch_ratios=ratios),indent=2))


def extraction_coverage_plot(config, branches):
    """Show method-specific validity and observed-profile support, not filled gaps."""
    out=config.output/'intermediate/x1d'
    fig,axes=plt.subplots(3,1,figsize=(13,8),sharex=True,constrained_layout=True)
    summary={}
    for ax,(name,data) in zip(axes,branches.items()):
        wave=data['WAVELENGTH'];summary[name]={}
        for method,color in (('BOX','tab:blue'),('OPT','tab:orange')):
            accepted=data['DQ_'+method]==0
            summary[name][method]=dict(total=len(wave),finite=int(np.isfinite(data['FLUX_'+method]).sum()),
                accepted=int(accepted.sum()),nan=int(np.sum(~np.isfinite(data['FLUX_'+method]))))
            ax.plot(wave,accepted.astype(float),color=color,lw=.7,alpha=.7,label=method+' accepted')
        if 'PROFCOV_OPT' in data:
            ax.plot(wave,data['PROFCOV_OPT'],color='green',lw=.6,label='Observed profile fraction')
            ax.plot(wave,data['INFOFRAC_OPT'],color='purple',lw=.6,label='Information fraction')
        ax.set(title=name,ylim=(-.05,1.05));ax.legend(fontsize=7,ncol=4)
    axes[-1].set_xlabel('Heliocentric wavelength [A]')
    fig.savefig(out/'extraction_coverage.png',dpi=150);plt.close(fig)
    (out/'coverage_summary.json').write_text(json.dumps(summary,indent=2))
