"""Six named workflow stages orchestrating reusable science modules.

1 prepare: validate FLTs, official controls, coordinates, diffuse calibration.
2 background: independent native-rate 2D models and model-selection diagnostics.
3 joint drizzle: subtract native model, calibrate, resample science once.
4 MAST coadd: combine original officially rectified x2d images on the same grid.
5 validate: photometry, coverage, integrity and source-independent diagnostics.
6 extract: consolidated FLT, joint-drizzle, and MAST-x2d 1D spectra.
"""

from pathlib import Path
import json
import sys
import numpy as np
from astropy.io import fits
from astropy import units as u
from .config import Config
from .io import read_inputs, validate, save_background, save_rectified_background, digest
from .reference import prepare, operator
from .geometry import make_grid, rectified_coordinates
from .background import estimate, add_residual, estimate_rectified, rectified_residual
from .coadd import Coadd
from .calibration import UNIT
from .diagnostics import exposure_plot, residual_plot, residual_summary, validate_products


def stage1_prepare(config, exposures):
    validate(exposures, config)
    prepare(exposures, config)
    directory = config.output / "intermediate/flt/registration"
    directory.mkdir(parents=True, exist_ok=True)
    headers = [
        fits.getheader(
            config.path(config.control_dir) / (e.path.stem + "_perform.fits"), "SCI"
        )
        for e in exposures
    ]
    offsets = [
        float(e.p["POSTARG2"]) + float(config.source_offsets.get(e.path.name, 0))
        for e in exposures
    ]
    grid = make_grid(headers, config, offsets)
    metadata = {}
    for e in exposures:
        print("  Prepare", e.path.name, flush=True)
        op = operator(e, config)
        op["mapping"] = grid.map(op["wavelength"], op["spatial"])
        np.savez_compressed(
            directory / (e.path.stem + ".npz"),
            **{k: v for k, v in op.items() if isinstance(v, np.ndarray)}
        )
        metadata[e.path.name] = dict(
            input_sha256=e.sha, postarg2=e.p["POSTARG2"], **op["meta"]
        )
    (directory / "inputs.json").write_text(
        json.dumps(
            dict(
                exposures=metadata,
                native_spatial_scale=headers[0]["CD2_2"] * 3600,
                grid=vars(grid),
            ),
            indent=2,
        )
    )
    return grid


def stage2_background(config, exposures):
    directory = config.output / "intermediate/background"
    directory.mkdir(parents=True, exist_ok=True)
    reports = {}
    for e in exposures:
        print("  Background", e.path.name, flush=True)
        with np.load(config.output / "intermediate/flt/registration" / (e.path.stem + ".npz")) as op:
            header = fits.getheader(
                config.path(config.control_dir) / (e.path.stem + "_perform.fits"), "SCI"
            )
            result = estimate(
                e,
                op["observed_wavelength"],
                op["source_y"],
                op["valid"],
                config.background,
                dispersion=header["CD1_1"],
            )
            result = add_residual(e, result, op["observed_wavelength"], op["source_y"],
                                  op["valid"], config.background, header["CD1_1"])
            save_background(directory / (e.path.stem + "_rate.fits"), e, result)
            exposure_plot(directory / (e.path.stem + ".png"), e, result, op["valid"])
            residual_plot(directory / (e.path.stem + "_residual.png"), e, result)
        reports[e.path.name] = result.metadata
    (directory / "models.json").write_text(json.dumps(reports, indent=2))
    return reports


def stage3_joint_drizzle(config, exposures, grid):
    directory = config.output / "intermediate/flt/coadd"
    directory.mkdir(parents=True, exist_ok=True)
    combined = Coadd(grid, config.pixfrac)
    baseline = Coadd(grid, config.pixfrac)
    primary_only = Coadd(grid, config.pixfrac)
    for e in exposures:
        print("  Joint drizzle", e.path.name, flush=True)
        with np.load(
            config.output / "intermediate/flt/registration" / (e.path.stem + ".npz")
        ) as op, fits.open(
            config.output / "intermediate/background" / (e.path.stem + "_rate.fits")
        ) as h:
            model = h["BKG"].data.astype(float) * e.time
            error = h["BKGERR"].data.astype(float) * e.time
            f = op["factor"]
            record = dict(
                input=e.path.name,
                sha256=e.sha,
                exptime=e.time,
                postarg2=float(e.p["POSTARG2"]),
                sky_width=float(h[0].header["SKYWIDTH"]),
                residual_applied=bool(h[0].header.get("RESSUB", False)),
            )
            kwargs = dict(
                mapping=op["mapping"],
                valid=op["valid"],
                covered=op["covered"],
                exptime=e.time,
                record=record,
            )
            combined.add(
                (e.counts - model) * f,
                e.variance * f**2,
                **kwargs,
                model=model * f,
                model_error=error * f,
                model_flags=h["BKGDQ"].data,
                input_dq=e.dq
            )
            baseline.add(e.counts * f, e.variance * f**2, **kwargs)
            primary = h["BKGPRIMARY"].data.astype(float) * e.time if "BKGPRIMARY" in h else model
            primary_only.add((e.counts - primary) * f, e.variance * f**2, **kwargs)
    np.savez_compressed(directory / "unsubtracted.npz", SCI=baseline.arrays()["SCI"])
    np.savez_compressed(directory / "primary_only.npz", SCI=primary_only.arrays()["SCI"])
    return combined.write(
        config.output / "products" / config.product_name("s2d_drz"),
        config,
        "CUSTOM_S2D_DRZ",
        background=True,
    )


def stage4_mast_coadd(config, exposures, grid):
    """Clean each MAST x2d, combine once, then assess the coadd residual."""
    from .diagnostics import rectified_background_plot, mast_background_report

    directory = config.output / "intermediate/x2d"
    directory.mkdir(parents=True, exist_ok=True)
    combined, unclean = Coadd(grid, 1.0), Coadd(grid, 1.0)
    options, reports, helio = config.mast_background, {}, []
    info = json.loads((config.output / "intermediate/flt/registration/inputs.json").read_text())
    native_scale = info['native_spatial_scale']
    for e in exposures:
        path = config.path(config.mast_dir) / e.path.name.replace("_flt", "_x2d")
        print("  MAST background + coadd", path.name, flush=True)
        with fits.open(path) as hd:
            if sum(x.name == "SCI" for x in hd) != 1:
                raise ValueError("MAST imsets must match single-imset FLTs")
            h = hd["SCI"].header
            science = hd["SCI"].data.astype(float)
            error = hd["ERR"].data.astype(float)
            dq = hd["DQ"].data.astype("u2")
            conversion = u.Unit(h["BUNIT"]).to(u.Unit(UNIT))
            if not np.isclose(conversion, 1):
                raise ValueError("Unexpected MAST SCI unit scaling")
            for key in ("OPT_ELEM", "DETECTOR", "CENWAVE", "PROPAPER"):
                if hd[0].header.get(key) != e.p.get(key):
                    raise ValueError("MAST/FLT configuration mismatch")
            for key in ("DARKCORR", "FLUXCORR"):
                if hd[0].header.get(key) != 'COMPLETE':
                    raise ValueError(f"MAST {key} must be COMPLETE")
            if not np.isclose(h["EXPTIME"], e.time):
                raise ValueError("MAST/FLT EXPTIME mismatch")
            wave, spatial = rectified_coordinates(h, science.shape, e.p["POSTARG2"],
                float(config.source_offsets.get(e.path.name, 0)))
            covered = (dq & 4) == 0  # archive no-data flag, not ERR-based selection
            valid = covered & (dq == 0) & np.isfinite(science)
            if np.any(valid & (~np.isfinite(error) | (error < 0))):
                raise ValueError("Unknown MAST uncertainty")
            hfactor = 1 - float(h['V_HELIO']) / 299792.458
            helio.append(hfactor)
            observed = wave / hfactor
            sy = spatial / native_scale
            result = estimate_rectified(science, error**2, observed, sy, valid,
                                        options, h['CD1_1']/hfactor)
            cleaned = directory / (path.stem + '_bkgsub.fits')
            save_rectified_background(cleaned, path, result)
            rectified_background_plot(directory / (path.stem + '.png'),
                science, result.rate, valid, wave[0], spatial[:,0], path.stem)
            reports[path.name] = dict(result.metadata, original_sha256=digest(path),
                                     output_sha256=digest(cleaned))
            mapping = grid.map(wave, spatial)
            record = dict(input=path.name, sha256=digest(path), exptime=e.time,
                          postarg2=float(e.p['POSTARG2']),
                          sky_width=result.metadata.get('selected_width', 0.),
                          residual_applied=False, bkgsub_sha256=digest(cleaned))
            kwargs = dict(mapping=mapping, valid=valid, covered=covered,
                          exptime=e.time, record=record, input_dq=dq)
            combined.add(science-result.rate, error**2, **kwargs,
                         model=result.rate, model_error=result.error, model_flags=result.flags)
            unclean.add(science, error**2, **kwargs)
    before = combined.arrays(background=True)
    unclean_arrays = unclean.arrays()
    np.savez_compressed(directory / 'unsubtracted.npz', SCI=unclean_arrays['SCI'])
    combined.write(directory / 'primary_only_coadd.fits', config,
                   'CUSTOM_MAST_X2D_COADD', background=True)
    wave = np.broadcast_to(grid.wave_start + np.arange(grid.nx)*grid.wave_step, grid.shape)
    sy = np.broadcast_to((grid.spatial_start + np.arange(grid.ny)[:,None]*grid.spatial_step)
                         / native_scale, grid.shape)
    # Common heliocentric grid has a tiny exposure-dependent airglow displacement.
    # Fail rather than conflate substantially different observed line geometries.
    hfactor = float(np.average(helio, weights=[e.time for e in exposures]))
    if np.ptp(helio) * np.max(wave) > .1 * grid.wave_step:
        raise ValueError('Airglow shifts require exposure-specific coadd residual geometry')
    valid = (before['DQ'] == 0) & np.isfinite(before['SCI'])
    residual = rectified_residual(before['SCI'], before['ERRBASE'], before['BKGERR'],
        wave/hfactor, sy, valid, options, grid.wave_step/hfactor,
        grid.spatial_step/native_scale)
    after = dict(before)
    after.update(SCI=before['SCI']-residual.rate,
        BKGPRIMARY=before['BKG'], PRIMARYERR=before['BKGERR'],
        BKGRESID=residual.rate, RESERR=residual.error, RESDQ=residual.flags,
        BKG=before['BKG']+residual.rate, BKGERR=before['BKGERR']+residual.error,
        ERR=np.hypot(before['ERRBASE'], before['BKGERR']+residual.error))
    # Same fixed weights must make model mapping commute with subtraction.
    domain = valid & np.isfinite(unclean_arrays['SCI'])
    discrepancy = np.max(abs((unclean_arrays['SCI']-after['BKG']-after['SCI'])[domain]))
    scale = max(float(np.max(abs(unclean_arrays['SCI'][domain]))), 1e-30)
    if discrepancy > 3e-6*scale:
        raise ValueError('MAST background/coadd linearity validation failed')
    reports['coadd_residual'] = dict(residual.metadata,
        subtraction_linearity_relative_error=float(discrepancy/scale),
        observed_frame_factor=hfactor)
    (directory/'models.json').write_text(json.dumps(reports, indent=2))
    final = combined.write(config.output/'products'/config.product_name('x2d_coadd'),
        config, 'CUSTOM_MAST_X2D_COADD', background=True, arrays=after,
        metadata=dict(COADDSUB=bool(residual.metadata['applied']),
                      RESSUB=bool(residual.metadata['applied']), BKGSTAGE='X2D_AND_COADD'))
    mast_background_report(config, grid, unclean_arrays['SCI'], before, after, reports)
    return final


def stage5_validate(config, products, grid, input_hashes, background_reports):
    """Validate both branches and save hashes for reproducibility."""
    report = validate_products(config, products, grid, input_hashes, background_reports)
    report["additional_residual"] = residual_summary(config, products[0], grid, background_reports)
    report["code_sha256"] = {
        p.name: digest(p) for p in Path(__file__).parent.glob("*.py")
    }
    report["reference_sha256"] = {
        p.name: digest(p)
        for p in config.path(config.reference_dir).glob("*")
        if p.is_file()
    }
    (config.output / "manifest.json").write_text(json.dumps(report, indent=2))
    return report


def stage6_extract(config):
    """Add 1D spectra without changing the established 2D calibration/drizzle."""
    from .extraction import run_extraction

    print("STAGE 6 — EXTRACT 1D", flush=True)
    return run_extraction(config)


def run(config):
    """Run the complete workflow; never overwrite final products implicitly."""
    if config.software_dir:
        software = config.path(config.software_dir)
        sys.path[:0] = [str(software / "drizzle-current"), str(software / "python")]
    exposures = read_inputs(sorted(Path(config.project_dir).glob(config.input_glob)))
    config.resolve_name(exposures)
    output = config.output
    products = output / "products"
    products.mkdir(parents=True, exist_ok=True)
    names = [config.product_name("s2d_drz"), config.product_name("x2d_coadd"), config.product_name("x1d")]
    if any((products / name).exists() for name in names):
        raise FileExistsError(
            "Final products exist; use a different output_dir for a new run"
        )
    hashes = {str(e.path): e.sha for e in exposures}
    hashes.update(
        {
            str(
                config.path(config.mast_dir) / e.path.name.replace("_flt", "_x2d")
            ): digest(
                config.path(config.mast_dir) / e.path.name.replace("_flt", "_x2d")
            )
            for e in exposures
        }
    )
    (output / "config.json").write_text(json.dumps(config.serializable(), indent=2))
    print("STAGE 1 — PREPARE", flush=True)
    grid = stage1_prepare(config, exposures)
    print("STAGE 2 — BACKGROUND", flush=True)
    reports = stage2_background(config, exposures)
    print("STAGE 3 — JOINT DRIZZLE", flush=True)
    clean = stage3_joint_drizzle(config, exposures, grid)
    print("STAGE 4 — MAST COADD", flush=True)
    mast = stage4_mast_coadd(config, exposures, grid)
    print("STAGE 5 — VALIDATE", flush=True)
    stage5_validate(config, [clean, mast], grid, hashes, reports)
    spectra = stage6_extract(config) if config.extraction.enabled else []
    print("Completed:", clean, mast, sep="\n", flush=True)
    return [clean, mast, *spectra]
