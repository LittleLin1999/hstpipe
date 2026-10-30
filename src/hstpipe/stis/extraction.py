"""Finite-aperture 2D extraction and official per-exposure CALSTIS x1d.

The 2D branch integrates diffuse surface brightness over slit area. The native
branch gives CALSTIS a COPY of each original FLT with SCI = COUNTS - t*BKG;
ERR, DQ, geometry, exposure time and calibration references are preserved.
The shared native-rate model includes any independently validated residual.
CALSTIS BACKCORR=OMIT prevents duplicate sky subtraction. Official point-source
FLUX is retained; FLUX_APER provides the diffuse convention for comparison.
"""

from pathlib import Path
from contextlib import contextmanager
import json
import os
import pkgutil
import sys
import tempfile
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from scipy.ndimage import map_coordinates
from .io import read_inputs, digest
from .calibration import UNIT, factor
from .optimal import wide_profile, continuum_center, profile_extract

FLUX_UNIT = "erg / (s cm2 Angstrom)"


def aperture_weights(centers, pixel_width, center, height):
    """Physical overlap length; divide by pixel_width for pixel fractions."""
    return np.maximum(
        0,
        np.minimum(np.asarray(centers) + pixel_width / 2, center + height / 2)
        - np.maximum(np.asarray(centers) - pixel_width / 2, center - height / 2),
    )


def bin_edges(wavelength):
    """Midpoint bin edges on a strictly increasing spectral grid."""
    wave = np.asarray(wavelength, float)
    if len(wave) < 2 or not np.isfinite(wave).all() or np.any(np.diff(wave) <= 0):
        raise ValueError("Wavelengths must be finite and strictly increasing")
    return np.r_[
        wave[0] - (wave[1] - wave[0]) / 2,
        (wave[:-1] + wave[1:]) / 2,
        wave[-1] + (wave[-1] - wave[-2]) / 2,
    ]


def overlap_matrix(input_wave, output_wave):
    """Flux-density redistribution: overlap length divided by OUTPUT bin width.

    No interpolation through gaps and no flux renormalization for missing data.
    Variances use the square of this matrix; correlated bounds use it linearly.
    """
    a, b = bin_edges(input_wave), bin_edges(output_wave)
    overlap = np.maximum(
        0, np.minimum(b[1:, None], a[None, 1:]) - np.maximum(b[:-1, None], a[None, :-1])
    )
    return overlap / np.diff(b)[:, None]


def sensitivity_domain(exposure, header, config):
    """Heliocentric range supported by PHOTTAB, excluding sensitivity extrapolation."""
    path = config.path(config.reference_dir) / exposure.p["PHOTTAB"].split("$")[-1]
    table = fits.getdata(path, 1)
    rows = table[
        (table["OPT_ELEM"] == exposure.p["OPT_ELEM"])
        & (table["CENWAVE"] == exposure.p["CENWAVE"])
        & (table["SPORDER"] == 1)
    ]
    if len(rows) != 1:
        raise ValueError("Ambiguous PHOTTAB selection")
    wave = rows[0]["WAVELENGTH"][: rows[0]["NELEM"]]
    hfactor = 1 - header["V_HELIO"] / 299792.458
    return float(np.min(wave) * hfactor), float(np.max(wave) * hfactor)


def calibrated_bins(wave, limits):
    """Accept a full spectral bin only when it lies inside reference support."""
    edges = bin_edges(wave)
    return (edges[:-1] >= limits[0]) & (edges[1:] <= limits[1])


def extract_image(
    science, error, model_error, valid, spatial, step, slit_width, center, height
):
    """Integrate surface brightness to F_lambda without multiplying by delta-lambda.

    ERR_DIAG ignores resampling covariance. ERR_STAT_BOUND uses the positive
    covariance upper bound sum(|weight| sigma); model uncertainty is likewise
    correlated, not added independently per spatial pixel. Incomplete apertures
    are flagged and not scaled up under an assumed unknown spatial profile.
    """
    weights = aperture_weights(spatial, step, center, height)
    if not np.isclose(weights.sum(), height, rtol=1e-6):
        raise ValueError("Extraction aperture falls outside the output spatial grid")
    area = weights[:, None] * slit_width
    coverage = np.sum(weights[:, None] * valid, axis=0) / height
    flux = np.sum(np.where(valid, science, 0) * area, axis=0)
    diag = np.sqrt(np.sum((np.where(valid, error, 0) * area) ** 2, axis=0))
    stat_bound = np.sum(np.where(valid, error, 0) * area, axis=0)
    bkg_bound = np.sum(np.where(valid, model_error, 0) * area, axis=0)
    good = coverage >= 1 - 1e-6
    result = dict(
        FLUX=flux,
        ERROR=np.hypot(stat_bound, bkg_bound),
        ERR_DIAG=diag,
        ERR_STAT_BOUND=stat_bound,
        BKGERR=bkg_bound,
        COVERAGE=coverage,
        DQ=np.where(good, 0, 4).astype("u2"),
    )
    for key in ("FLUX", "ERROR", "ERR_DIAG", "ERR_STAT_BOUND", "BKGERR"):
        result[key] = np.where(good, result[key], np.nan)
    return result


def write_spectrum(path, columns, metadata, provenance, descriptions=()):
    """Write a custom one-row-per-wavelength FITS table, with explicit conventions."""
    table = Table(columns)
    for name in table.colnames:
        if name == "WAVELENGTH":
            table[name].unit = "Angstrom"
        elif name.startswith("COVAR"):
            table[name].unit = u.Unit(FLUX_UNIT) ** 2
        elif name in (
            "FLUX",
            "FLUX_OPT",
            "FLUX_OPT_APER",
            "FLUX_APER",
            "ERROR",
            "ERR_OPT",
            "ERR_OPT_APER",
            "ERR_STAT_OPT",
            "ERR_STAT_OPT_APER",
            "BKGERR_OPT",
            "BKGERR_OPT_APER",
            "ERROR_APER",
            "ERR_DIAG",
            "ERR_STAT_BOUND",
            "ERR_STAT",
            "ERR_STAT_APER",
            "BKGERR",
            "BKGERR_APER",
            "ERROR_CALSTIS",
        ):
            table[name].unit = FLUX_UNIT
        elif name.startswith("EXPTIME"):
            table[name].unit = "s"
    header = fits.Header()
    header["OFFICIAL"] = False
    header["BACKCORR"] = "PERFORM" if metadata.get("BACKCORR") == "PERFORM" else "OMIT"
    for name, value in metadata.items():
        header[name] = value
    for line in descriptions:
        header.add_history(line)
    hd = fits.HDUList(
        [
            fits.PrimaryHDU(header=header),
            fits.BinTableHDU(table, name="SCI"),
            fits.BinTableHDU(
                Table(rows=[
                    {key: json.dumps(value, sort_keys=True) if isinstance(value, dict) else value
                     for key, value in record.items()}
                    for record in provenance
                ]),
                name="PROVENANCE",
            ),
        ]
    )
    hd.writeto(path, overwrite=False, checksum=True)


def extract_2d(config, path, native_scale, slit_width, calibration_limits=None):
    with fits.open(path) as h:
        sh = h["SCI"].header
        if not np.isclose(u.Unit(sh["BUNIT"]).to(u.Unit(UNIT)), 1):
            raise ValueError("2D input must contain diffuse surface brightness")
        sci = h["SCI"].data.astype(float)
        error = h["ERRBASE"].data.astype(float) if "ERRBASE" in h else h["ERR"].data.astype(float)
        bkgerr = h["BKGERR"].data.astype(float) if "BKGERR" in h else np.zeros_like(error)
        valid = (h["DQ"].data == 0) & np.isfinite(sci)
        if np.any(valid & (~np.isfinite(error) | ~np.isfinite(bkgerr))):
            raise ValueError("Unknown uncertainty in usable 2D pixels")
        spatial = (
            sh["CRVAL2"] + (np.arange(sci.shape[0]) + 1 - sh["CRPIX2"]) * sh["CD2_2"]
        )
        wave = sh["CRVAL1"] + (np.arange(sci.shape[1]) + 1 - sh["CRPIX1"]) * sh["CD1_1"]
        if calibration_limits is not None:
            valid &= calibrated_bins(wave, calibration_limits)[None, :]
        options = config.extraction
        height = options.height_native_pixels * native_scale
        result = extract_image(
            sci,
            error,
            bkgerr,
            valid,
            spatial,
            sh["CD2_2"],
            slit_width,
            options.center_arcsec,
            height,
        )
        fraction = aperture_weights(spatial, sh["CD2_2"], options.center_arcsec, options.optimal_height_native_pixels*native_scale) / sh["CD2_2"]
        aperture = np.broadcast_to(fraction[:, None], sci.shape)
        fit_region = aperture_weights(spatial, sh['CD2_2'], options.center_arcsec,
                                      options.profile_height_native_pixels*native_scale)[:,None]
        profile, profile_meta = wide_profile(sci, valid, wave, fit_region,
                                             sh['CD2_2']/native_scale)
        profile_meta.update(boxcar_height_native=options.height_native_pixels,
                            optimal_height_native=options.optimal_height_native_pixels,
                            fitting_height_native=options.profile_height_native_pixels,
                            min_profile_fraction=options.optimal_min_profile_fraction,
                            min_information_fraction=options.optimal_min_information_fraction)
        opt = profile_extract(
            sci * sh["CD2_2"] * slit_width,
            error * sh["CD2_2"] * slit_width,
            bkgerr * sh["CD2_2"] * slit_width,
            valid,
            aperture,
            profile,
            covariance_bound=True,
            min_profile_fraction=options.optimal_min_profile_fraction,
            min_information_fraction=options.optimal_min_information_fraction,
        )
        from .diagnostics import extraction_fit_plot
        extraction_fit_plot(config.output/'diagnostics/extraction',
            'DRZ' if path.name.endswith('_s2d_drz.fits') else 'X2D',
            sci * sh['CD2_2'] * slit_width, error * sh['CD2_2'] * slit_width,
            bkgerr * sh['CD2_2'] * slit_width, valid, aperture, profile, opt['FLUX_OPT'],
            wave, (spatial-options.center_arcsec)/native_scale, profile_meta, FLUX_UNIT,
            correlated=True)
        result["DQ_BOX"] = result["DQ"].copy()
        result.update(opt)
        box_rows = aperture_weights(spatial, sh['CD2_2'], options.center_arcsec,height)>0
        for method, rows in (("BOX", box_rows), ("OPT", fraction > 0)):
            used = valid[rows]
            # Conservative minimum among surviving pixels. These metadata do
            # not set flux weights; EXPMAP is not a new inverse-noise estimate.
            for key, extension in (("NEXP", "NEXP"), ("EXPTIME", "EXPMAP")):
                values=np.where(used,h[extension].data[rows],np.inf)
                minimum=np.min(values,axis=0)
                result[key+"_"+method]=np.where(np.isfinite(minimum),minimum,0)
            result['BKGDQ'+('_OPT' if method=='OPT' else '')] = (
                np.bitwise_or.reduce(np.where(used,h['BKGDQ'].data[rows],0),axis=0)
                if 'BKGDQ' in h else np.zeros(sci.shape[1],'u2'))
        result["NEXP"] = result["NEXP_BOX"].astype('u2')
        result["EXPTIME"] = result["EXPTIME_BOX"]
        result["WAVELENGTH"] = wave
    return result


@contextmanager
def official_environment(config):
    """Configure the existing isolated STIS installation; restore process env."""
    old = {key: os.environ.get(key) for key in ("oref", "PATH")}
    try:
        os.environ["oref"] = str(config.path(config.reference_dir)) + "/"
        if config.software_dir:
            software = config.path(config.software_dir)
            sys.path.insert(0, str(software / "python"))
            os.environ["PATH"] = (
                str(software / "hstcal-install/bin") + os.pathsep + os.environ["PATH"]
            )
            # Legacy setuptools namespace packages can hide the isolated stsci
            # packages. Use the standard namespace-path extension API.
            import stsci

            stsci.__path__ = pkgutil.extend_path(stsci.__path__, stsci.__name__)
        from stistools import x1d

        yield x1d.x1d
    finally:
        for key, value in old.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def count_copy(original, rate_model, destination):
    """Copy the ORIGINAL FLT layout; change SCI only, retain count-valued ERR.

    This is not the custom rate FITS relabelled as an FLT. Do not change EXPTIME,
    detector coordinates, wavecal offsets, ERR or DQ. Correlated model errors
    are propagated outside CALSTIS, whose image ERR expects diagonal variance.
    """
    with fits.open(original) as hd, fits.open(rate_model) as b:
        t = float(hd["SCI"].header["EXPTIME"])
        if hd["SCI"].header.get("BUNIT", "").upper() not in ("COUNT", "COUNTS"):
            raise ValueError("x1d copy requires original count-valued SCI")
        if (
            hd[0].header.get("ROOTNAME") != b[0].header.get("ROOTNAME")
            or hd["SCI"].data.shape != b["BKG"].data.shape
        ):
            raise ValueError("Background does not match the original exposure")
        if not np.isclose(u.Unit(b["BKG"].header["BUNIT"]).to(u.count / u.s), 1):
            raise ValueError("Native background must be a count rate")
        hd["SCI"].data = (hd["SCI"].data.astype(float) - t * b["BKG"].data).astype("f4")
        hd[0].header["BACKCORR"] = "OMIT"
        hd[0].header["RESSUB"] = bool(b[0].header.get("RESSUB", False))
        hd[0].header["BKG2DSUB"] = bool(b[0].header["BKG2DSUB"])
        hd[0].header["FILENAME"] = destination.name
        hd[0].header.add_history(
            "Custom SCI=original COUNTS-EXPTIME*BKG_RATE. No resampling."
        )
        hd[0].header.add_history(
            "Original ERR retained; correlated BKGERR propagated separately in 1D."
        )
        hd.writeto(destination, overwrite=False, checksum=True)


def pctab_height(path, aperture, cenwave, height, require_exact=True):
    """Audit HSTCAL GetPCT6's nearest-height choice before running x1d.

    HSTCAL interpolates in wavelength, NOT extraction height. Point products
    require a tabulated height rather than silently using another aperture.
    Diffuse/morph products cancel this point response, but record its selection.
    """
    with fits.open(path) as hd:
        table=hd[1].data
        match=(np.char.strip(table['APERTURE'])==aperture) & (
            (table['CENWAVE']==cenwave) | (table['CENWAVE']==-1))
        rows=table[match]
        if not len(rows):
            raise ValueError(f'No PCTAB rows for {aperture}, CENWAVE={cenwave}')
        heights=rows['EXTRHEIGHT'].astype(float)
        selected=int(np.argmin(abs(heights-int(height))))
        actual=float(heights[selected])
        if require_exact and not np.isclose(height,actual,rtol=0,atol=1e-8):
            raise ValueError(f'Point extraction height {height:g} has no exact PCTAB entry; '
                             f'HSTCAL would use {actual:g}. Choose one of {sorted(set(heights))}')
        if 'PEDIGREE' in table.names and str(rows[selected]['PEDIGREE']).startswith('DUMMY'):
            raise ValueError('DUMMY PCTAB row cannot calibrate a science extraction')
        return dict(requested_height=float(height),pctab_height=actual,
                    exact=bool(np.isclose(height,actual)),reference=Path(path).name,
                    reference_sha256=digest(path))


def official_exposure(config, e, directory, x1d):
    """Run CALSTIS-6 and validate its actual trace, net signal and calibration."""
    options = config.extraction
    pctab=config.path(config.reference_dir)/e.p['PCTAB'].split('$')[-1]
    height_audit={method:pctab_height(pctab,e.p['APERTURE'],e.p['CENWAVE'],height,
                                    require_exact=options.flux_mode=='point')
                  for method,height in [('BOX',options.height_native_pixels),
                                        ('OPT',options.optimal_height_native_pixels)]}
    rate_path = config.output / "intermediate/background" / (e.path.stem + "_rate.fits")
    native = directory / (e.path.stem + "_bkgsub_counts.fits")
    output = directory / (e.path.stem.replace("_flt", "") + "_bkgsub_x1d.fits")
    count_copy(e.path, rate_path, native)
    h = fits.getheader(
        config.path(config.control_dir) / (e.path.stem + "_perform.fits"), "SCI"
    )
    scale = h["CD2_2"] * 3600
    with np.load(config.output / "intermediate/flt/registration" / (e.path.stem + ".npz")) as op:
        sy = op["source_y"].copy()
        covered = op["covered"].copy()
    center = options.center_arcsec / scale
    guess = (
        np.interp(center, sy[:, e.counts.shape[1] // 2], np.arange(e.counts.shape[0]))
        + 1
    )
    kwargs = dict(
        backcorr="omit",
        ctecorr="omit",
        dispcorr="perform",
        helcorr="perform",
        fluxcorr="perform",
        maxsrch=0,
        extrsize=options.height_native_pixels,
        algorithm="unweighted",
    )

    def invoke(input_path, output_path, position, height=None):
        status = x1d(
            str(input_path),
            output=str(output_path),
            a2center=float(position),
            trailer=str(output_path.with_suffix(".log")),
            **dict(kwargs, extrsize=options.height_native_pixels if height is None else height),
        )
        if status != 0 or not output_path.exists():
            raise RuntimeError(f"CALSTIS x1d failed: {output_path}")
        with fits.open(output_path) as hd:
            if len(hd["SCI"].data) != 1 or hd[0].header.get("FLUXCORR") != "COMPLETE":
                raise ValueError("Expected one calibrated G140L spectral order")
            return {
                name: np.array(hd["SCI"].data[0][name], copy=True)
                for name in hd["SCI"].columns.names
            }

    # Positive unit-SCI control measures the official wavelength-dependent
    # calibration even where the real NET is zero/negative. It also verifies
    # that the CALSTIS extraction trace matches the established header mapping.
    with tempfile.TemporaryDirectory(prefix="x1d_control_", dir=directory) as tmp:
        tmp = Path(tmp)
        unit = tmp / "unit.fits"
        with fits.open(native) as hd:
            hd["SCI"].data[:] = 1.0
            hd.writeto(unit, checksum=True)
        for iteration in range(3):
            control = invoke(unit, tmp / f"control{iteration}.fits", guess)
            trace = control["EXTRLOCY"].astype(float) - 1
            columns = np.arange(len(trace))
            relative_y = map_coordinates(sy, [trace, columns], order=1, mode="nearest")
            trace_valid = (control["DQ"] == 0) & map_coordinates(
                covered.astype(float), [trace, columns], order=0
            ).astype(bool)
            if not trace_valid.any():
                raise ValueError("No valid trace control pixels")
            offset = float(np.median(relative_y[trace_valid] - center))
            if abs(offset) < 1e-3:
                break
            if iteration == 2:
                raise ValueError("CALSTIS trace registration did not converge")
            guess -= offset
        trace_error = float(np.max(abs(relative_y[trace_valid] - center)))
        if trace_error > options.trace_tolerance_pixels:
            raise ValueError(
                f"CALSTIS/x2d trace mismatch: {trace_error:g} native pixels"
            )
        response = control["FLUX"].astype(float) / control["NET"]
        if not np.all(np.isfinite(response[trace_valid]) & (response[trace_valid] > 0)):
            raise ValueError("Invalid official sensitivity control")
        wide_control = invoke(unit, tmp/'wide_control.fits', guess,
                              height=options.optimal_height_native_pixels)
        response_opt = wide_control['FLUX'].astype(float)/wide_control['NET']
        if not np.all(np.isfinite(response_opt[trace_valid]) & (response_opt[trace_valid]>0)):
            raise ValueError('Invalid official sensitivity at optimal extraction height')
        row = invoke(native, output, guess)

    np.testing.assert_allclose(row["EXTRLOCY"], control["EXTRLOCY"], rtol=0, atol=1e-5)
    np.testing.assert_array_equal(row["WAVELENGTH"], control["WAVELENGTH"])
    weights = aperture_weights(
        np.arange(e.counts.shape[0])[:, None],
        1,
        row["EXTRLOCY"][None, :] - 1,
        options.height_native_pixels,
    )
    with fits.open(native) as n, fits.open(rate_path) as b:
        native_counts = n["SCI"].data.astype(float)
        bad = (
            (e.dq != 0)
            | ~covered
            | ~np.isfinite(native_counts)
            | ~np.isfinite(e.variance)
        )
        full = (weights.sum(0) >= options.height_native_pixels - 1e-6) & ~np.any(
            (weights > 0) & bad, axis=0
        )
        calibration_limits = sensitivity_domain(e, h, config)
        cal_good = calibrated_bins(row["WAVELENGTH"], calibration_limits)
        good = full & (row["DQ"] == 0) & cal_good
        direct_net = np.sum(weights * np.nan_to_num(native_counts), axis=0) / e.time
        residual = row["BACKGROUND"].astype(float)
        net_difference = float(np.max(abs(direct_net[good] - residual[good] - row["NET"][good])))
        # EXTRLOCY is stored as float32, although CALSTIS sums using its
        # internal double-precision trace. Bound the resulting edge-area
        # difference using the actual boundary counts, not a loose flux rtol.
        edge = (weights > 0) & (weights < 1)
        quantization = (
            0.5
            * np.spacing(row["EXTRLOCY"])
            * np.sum(np.where(edge, abs(native_counts), 0), axis=0)
            / e.time
        )
        rounding = (
            8
            * np.finfo("f4").eps
            * np.sum(weights * abs(native_counts), axis=0)
            / e.time
        )
        if np.any(
            abs(row["NET"][good] + residual[good] - direct_net[good])
            > quantization[good] + rounding[good] + 1e-10
        ):
            raise ValueError(
                "Official NET plus residual background differs from native aperture sum"
            )
        np.testing.assert_allclose(
            row["FLUX"][good], row["NET"][good] * response[good], rtol=3e-6, atol=1e-30
        )
        # Recompute fractional-aperture variance with squared weights; retain
        # CALSTIS ERROR separately (cs6 accumulates area*ERR^2 at edge pixels).
        stat_rate = (
            np.sqrt(np.sum(weights**2 * np.nan_to_num(e.variance), axis=0)) / e.time
        )
        model_rate_bound = np.sum(weights * b["BKGERR"].data, axis=0)
        model_flags = np.bitwise_or.reduce(
            np.where(weights > 0, b["BKGDQ"].data, 0), axis=0
        )
        residual_metadata = json.loads((config.output / "intermediate/background/models.json").read_text())[e.path.name].get("residual", {})
        audit = dict(method="shared_native_unclipped", backcorr="OMIT",
                     decision=residual_metadata)
        opt_weights = aperture_weights(np.arange(e.counts.shape[0])[:,None],1,
            row['EXTRLOCY'][None,:]-1,options.optimal_height_native_pixels)
        fit_weights = aperture_weights(np.arange(e.counts.shape[0])[:,None],1,
            row['EXTRLOCY'][None,:]-1,options.profile_height_native_pixels)
        bkg_flags = b["BKGDQ"].data.copy()
        valid_pixels = ~bad & cal_good[None, :]
        profile, profile_meta = wide_profile(native_counts/e.time, valid_pixels,
            row['WAVELENGTH'], fit_weights)
        profile_meta.update(boxcar_height_native=options.height_native_pixels,
                            optimal_height_native=options.optimal_height_native_pixels,
                            fitting_height_native=options.profile_height_native_pixels,
                            min_profile_fraction=options.optimal_min_profile_fraction,
                            min_information_fraction=options.optimal_min_information_fraction)
        opt = profile_extract(
            native_counts / e.time,
            np.sqrt(e.variance) / e.time,
            b["BKGERR"].data.astype(float),
            valid_pixels,
            opt_weights,
            profile,
            min_profile_fraction=options.optimal_min_profile_fraction,
            min_information_fraction=options.optimal_min_information_fraction,
        )
        from .diagnostics import extraction_fit_plot
        extraction_fit_plot(config.output/'diagnostics/extraction',
            'FLT_' + e.path.stem.replace('_flt',''), native_counts/e.time,
            np.sqrt(e.variance)/e.time, b['BKGERR'].data.astype(float), ~bad & cal_good[None,:],
            opt_weights, profile, opt['FLUX_OPT'], row['WAVELENGTH'],
            np.arange(e.counts.shape[0])-np.median(row['EXTRLOCY']-1), profile_meta,
            'count / s', correlated=False)
    conversion, meta = factor(
        e,
        Path(config.project_dir),
        row["WAVELENGTH"],
        h,
        reference_dir=config.path(config.reference_dir),
    )
    aperture_response = (
        conversion * e.time * meta["slit_width"] * meta["native_spatial_scale"]
    )
    official_stat_rate = np.divide(
        row["ERROR"].astype(float), response,
        out=np.full_like(response, np.nan), where=response > 0,
    )
    # This column adopts the SAME native-sampling diffuse convention as x2d.
    # The ratio includes aperture/throughput and dispersion-convention effects;
    # do not call it a pure point-source aperture correction.
    result = dict(
        WAVELENGTH=row["WAVELENGTH"],
        FLUX=row["FLUX"],
        FLUX_APER=row["NET"] * aperture_response,
        ERR_STAT=np.maximum(stat_rate * response, row["ERROR"].astype(float)),
        BKGERR=model_rate_bound * response,
        ERR_STAT_APER=np.maximum(stat_rate, official_stat_rate) * aperture_response,
        BKGERR_APER=model_rate_bound * aperture_response,
        ERROR_CALSTIS=row["ERROR"],
        CAL_RATIO=np.divide(
            aperture_response,
            response,
            out=np.full_like(response, np.nan),
            where=response > 0,
        ),
        DQ=np.where(good, 0, np.asarray(row["DQ"], "u2") | 4).astype("u2"),
        BKGDQ=model_flags,
        CALDQ=(~cal_good).astype("u2"),
    )
    result["ERROR"] = np.hypot(result["ERR_STAT"], result["BKGERR"])
    result["ERROR_APER"] = np.hypot(result["ERR_STAT_APER"], result["BKGERR_APER"])
    # The profile-weighted estimate has its own fractional-pixel variance.
    # Retain at least the official boxcar statistical scale until a calibrated
    # per-pixel MAMA variance model is available.
    # Both science branches have already subtracted exactly the same total
    # native model. CALSTIS's BACKGROUND must therefore remain zero.
    np.testing.assert_allclose(residual, 0, atol=0)
    selected_net = row["NET"]
    result["FLUX"] = selected_net * response
    result["FLUX_APER"] = selected_net * aperture_response
    opt_net = opt["FLUX_OPT"]
    result["FLUX_OPT"] = opt_net * response_opt
    result["FLUX_OPT_APER"] = opt_net * aperture_response
    result["ERR_STAT_OPT"] = np.maximum.reduce((opt["ERR_STAT_OPT"], stat_rate, official_stat_rate)) * response_opt
    result["BKGERR_OPT"] = opt["BKGERR_OPT"] * response_opt
    result["ERR_OPT"] = np.hypot(result["ERR_STAT_OPT"], result["BKGERR_OPT"])
    result["CAL_RATIO_OPT"] = np.divide(
        aperture_response, response_opt,
        out=np.full_like(response_opt, np.nan), where=response_opt > 0)
    result["ERR_STAT_OPT_APER"] = result["ERR_STAT_OPT"] * result["CAL_RATIO_OPT"]
    result["BKGERR_OPT_APER"] = result["BKGERR_OPT"] * result["CAL_RATIO_OPT"]
    result["ERR_OPT_APER"] = np.hypot(result["ERR_STAT_OPT_APER"], result["BKGERR_OPT_APER"])
    result["DQ_BOX"] = result["DQ"].copy()
    result["DQ_OPT"] = opt["DQ_OPT"].copy()
    for key in ("PROFCOV_OPT", "INFOFRAC_OPT", "OPT_MASKED"):
        result[key] = opt[key]
    result['BKGDQ_OPT'] = np.bitwise_or.reduce(
        np.where(valid_pixels & (opt_weights>0), bkg_flags, 0),axis=0)
    record = dict(
        input=e.path.name,
        input_sha256=e.sha,
        exptime=e.time,
        a2center=float(guess),
        boxcar_height=options.height_native_pixels,
        optimal_height=options.optimal_height_native_pixels,
        profile_height=options.profile_height_native_pixels,
        optimal_point_calibration="official unit-SCI control at optimal height",
        pctab_height_audit=height_audit,
        trace_max_difference_pixels=trace_error,
        net_max_difference=net_difference,
        residual_background_median=float(np.nanmedian(residual[good])),
        residual_background_max=float(np.nanmax(abs(residual[good]))),
        residual_audit=audit,
        residual_applied=bool(residual_metadata.get("applied", False)),
        residual_method="SHARED_NATIVE_UNCLIPPED",
        background_sha256=digest(rate_path),
        count_input_sha256=digest(native),
        official_x1d_sha256=digest(output),
        calibration_limits=calibration_limits,
    )
    write_spectrum(
        directory / (e.path.stem + "_spectrum.fits"),
        result,
        dict(PRODTYPE="CUSTOM_X1D_WITH_ERRORS", EXTRPIX=options.height_native_pixels,
             BACKCORR="OMIT",
             RESBKG="APPLIED" if residual_metadata.get("applied", False) else "NOT_APPLIED"),
        [record],
        [
            "FLUX: official CALSTIS point-source calibration.",
            "FLUX_APER: x2d diffuse convention, finite detector aperture.",
            "ERROR_CALSTIS retained; ERR_STAT uses squared fractional weights.",
            "BKGERR is a correlated bound, not independent pixel noise.",
            "CALSTIS BACKCORR=OMIT; shared validated native residual already included.",
        ],
    )
    return result, record


def combine_spectra(spectra, times, output_wave):
    """Coadd each extraction method with its own mask and exposure weights.

    Whole output wavelength bins still require valid input spectral coverage;
    missing wavelengths are never interpolated. Spatially missing pixels are
    handled only by the optimal profile fit, not by this spectral coaddition.
    """
    n = len(output_wave)
    result = dict(WAVELENGTH=np.asarray(output_wave))
    for method, suffix in (("BOX", ""), ("OPT", "_OPT")):
        names = [base+suffix+aper for base in ("FLUX","ERR_STAT","BKGERR")
                 for aper in ("", "_APER")]
        sums = {key:np.zeros(n) for key in names}
        denominator=np.zeros(n);nexp=np.zeros(n,'u2')
        rejected=np.zeros(n,'u2');bkgdq=np.zeros(n,'u2');caldq=np.zeros(n,'u2')
        covariance=np.zeros(n-1)
        support={key:np.zeros(n) for key in ('PROFCOV_OPT','INFOFRAC_OPT','OPT_MASKED')}
        for spectrum,time in zip(spectra,times):
            matrix=overlap_matrix(spectrum['WAVELENGTH'],output_wave)
            footprint=matrix>0
            dq=spectrum.get('DQ_'+method,spectrum['DQ'])
            good=(dq==0)&np.isfinite(spectrum['FLUX'+suffix])&np.isfinite(spectrum['FLUX'+suffix+'_APER'])
            for key in names:
                if key.startswith(('ERR_STAT','BKGERR')) and np.any(good & (
                    ~np.isfinite(spectrum[key]) | (spectrum[key]<0))):
                    raise ValueError(f'Unknown or negative {key} in accepted {method} bins')
            matrix[:,~good]=0
            complete=matrix.sum(axis=1)>=1-1e-6
            matrix[~complete]=0
            denominator+=time*complete;nexp+=complete
            for key in names:
                values=np.where(good,spectrum[key],0)
                sums[key]+=time**2*(matrix**2@values**2) if key.startswith('ERR_STAT') else time*(matrix@values)
            variance=np.where(good,spectrum['ERR_STAT'+suffix+'_APER']**2,0)
            covariance+=np.sum(matrix[:-1]*matrix[1:]*variance,axis=1)*time**2
            for j,hits in enumerate(footprint):
                if hits.any():
                    rejected[j]|=np.bitwise_or.reduce(dq[hits])
                    caldq[j]|=np.bitwise_or.reduce(spectrum.get('CALDQ',np.zeros_like(dq))[hits])
                if complete[j]:
                    bkgdq[j]|=np.bitwise_or.reduce(spectrum.get('BKGDQ'+suffix,spectrum['BKGDQ'])[hits])
            if method=='OPT':
                for key in support:
                    default=np.zeros_like(good,float) if key=='OPT_MASKED' else np.ones_like(good,float)
                    support[key]+=time*(matrix@np.where(good,spectrum.get(key,default),0))
        accepted=denominator>0
        for key,values in sums.items():
            if key.startswith('ERR_STAT'):values=np.sqrt(values)
            result[key]=np.divide(values,denominator,out=np.full(n,np.nan),where=accepted)
        for aper in ('','_APER'):
            error_key=('ERROR' if method=='BOX' else 'ERR_OPT')+aper
            result[error_key]=np.hypot(result['ERR_STAT'+suffix+aper],result['BKGERR'+suffix+aper])
        result['DQ_'+method]=np.where(accepted,0,4).astype('u2')
        result['EXPTIME_'+method]=denominator
        result['NEXP_'+method]=nexp
        for key,values in (('REJDQ',rejected),('BKGDQ',bkgdq),('CALDQ',caldq)):
            result[key+suffix]=values
        pair=denominator[:-1]*denominator[1:]
        result['COVAR1_STAT'+suffix+'_APER']=np.r_[np.divide(covariance,pair,out=np.full(n-1,np.nan),where=pair>0),np.nan]
        if method=='OPT':
            for key,values in support.items():
                mean=np.divide(values,denominator,out=np.full(n,np.nan),where=accepted)
                result[key]=(mean>0) if key=='OPT_MASKED' else mean
    # Historical unsuffixed columns refer to boxcar only, never the mask union.
    for key in ('DQ','EXPTIME','NEXP'):result[key]=result[key+'_BOX'].copy()
    return result


def run_extraction(config):
    """Add 1D products to a validated reduction without rerunning/changing its 2D."""
    exposures = read_inputs(sorted(Path(config.project_dir).glob(config.input_glob)))
    config.resolve_name(exposures)
    directory = config.output / "intermediate/x1d/flt"
    spectra_dir = config.output / "intermediate/x1d"
    directory.mkdir(parents=True, exist_ok=True)
    spectra_dir.mkdir(exist_ok=True)
    final = config.output / "products" / config.product_name("x1d")
    outputs = [
        spectra_dir / (config.name + suffix)
        for suffix in ("_from2d_x1d.fits", "_individual_x1d_coadd.fits")
    ]
    if final.exists() or any(p.exists() for p in outputs) or any(directory.glob("*_bkgsub_counts.fits")):
        raise FileExistsError(
            "Extraction outputs exist; preserve them or select a new reduction directory"
        )
    manifest_path = config.output / "manifest.json"
    parent = json.loads(manifest_path.read_text())
    previous_config = json.loads((config.output / "config.json").read_text())
    current = json.loads(json.dumps(config.serializable()))
    for key in previous_config:
        if key != "extraction" and previous_config[key] != current[key]:
            raise ValueError(
                f"Extraction config differs from existing 2D reduction: {key}"
            )
    products = config.output / "products"
    for name, expected in parent["products"].items():
        if digest(products / name) != expected:
            raise ValueError("2D checksum differs from parent manifest")
    stage = json.loads((config.output / "intermediate/flt/registration/inputs.json").read_text())
    for e in exposures:
        if stage["exposures"][e.path.name]["input_sha256"] != e.sha:
            raise ValueError("Original FLT changed since 2D reduction")
        for key in ("XTRACTAB", "APERTAB", "PCTAB", "SPTRCTAB"):
            if not (
                config.path(config.reference_dir) / e.p[key].split("$")[-1]
            ).is_file():
                raise FileNotFoundError(e.p[key])
    native_scale = stage["native_spatial_scale"]
    widths = [v["slit_width"] for v in stage["exposures"].values()]
    if not np.allclose(widths, widths[0]):
        raise ValueError("Different slit widths unsupported")
    options = config.extraction
    if (
        abs(options.center_arcsec) / native_scale + options.height_native_pixels / 2
        > config.background.source_radius
    ):
        raise ValueError(
            "Aperture extends outside the protected background source region"
        )
    clean = products / config.product_name("s2d_drz")
    mast = products / config.product_name("x2d_coadd")
    domains = [
        sensitivity_domain(
            e,
            fits.getheader(
                config.path(config.control_dir) / (e.path.stem + "_perform.fits"), "SCI"
            ),
            config,
        )
        for e in exposures
    ]
    common_limits = max(x[0] for x in domains), min(x[1] for x in domains)
    if options.center_method == 'continuum':
        with fits.open(clean) as hd:
            hh=hd['SCI'].header;ss=hd['SCI'].data.astype(float)
            ww=hh['CRVAL1']+(np.arange(ss.shape[1])+1-hh['CRPIX1'])*hh['CD1_1']
            yy=hh['CRVAL2']+(np.arange(ss.shape[0])+1-hh['CRPIX2'])*hh['CD2_2']
            vv=(hd['DQ'].data==0)&np.isfinite(ss)&calibrated_bins(ww,common_limits)[None,:]
            centering=continuum_center(ss,vv,ww,yy,options.center_arcsec,native_scale)
        options.center_arcsec=centering['center_arcsec']
    else:
        centering=dict(method='configured_header_center',center_arcsec=options.center_arcsec)
    centering['source_file']=clean.name
    centering['source_sha256']=digest(clean)
    (spectra_dir/'centering.json').write_text(json.dumps(centering,indent=2))
    if abs(options.center_arcsec)/native_scale + options.profile_height_native_pixels/2 > config.background.source_radius:
        raise ValueError('Wide extraction/profile extends outside the protected source region')
    from2d = extract_2d(config, clean, native_scale, widths[0], common_limits)
    frommast = extract_2d(config, mast, native_scale, widths[0], common_limits)
    np.testing.assert_allclose(from2d["WAVELENGTH"], frommast["WAVELENGTH"], rtol=0, atol=1e-9)
    calibrated = calibrated_bins(from2d["WAVELENGTH"], common_limits)
    for branch in (from2d, frommast):
        branch["CALDQ"] = (~calibrated).astype("u2")
        branch["DQ"][~calibrated] |= 4
        branch["DQ_BOX"][~calibrated] |= 4
        branch["DQ_OPT"][~calibrated] |= 4
        for key in ("FLUX", "ERROR", "ERR_DIAG", "ERR_STAT_BOUND", "BKGERR",
                    "FLUX_OPT", "ERR_OPT", "ERR_STAT_OPT", "BKGERR_OPT"):
            branch[key][~calibrated] = np.nan
    extracted, records = [], []
    with official_environment(config) as x1d:
        for e in exposures:
            print("  Official x1d:", e.path.name, flush=True)
            spectrum, record = official_exposure(config, e, directory, x1d)
            extracted.append(spectrum)
            records.append(record)
    combined = combine_spectra(
        extracted, [e.time for e in exposures], from2d["WAVELENGTH"]
    )
    from .products import (
        point_conversion, morphology_conversion, branch_columns, write_consolidated
    )
    if options.flux_mode == "point":
        correction = {key:point_conversion(extracted,[e.time for e in exposures],
            from2d['WAVELENGTH'],overlap_matrix,ratio_key=ratio)
            for key,ratio in [('BOX','CAL_RATIO'),('OPT','CAL_RATIO_OPT')]}
    elif options.flux_mode == "morph":
        correction = {key:morphology_conversion(config.path(options.source_model),
            from2d['WAVELENGTH'],fraction_key='APFRAC_'+key) for key in ('BOX','OPT')}
    else:
        correction = np.ones_like(from2d["WAVELENGTH"])
    flt_columns = branch_columns(combined, options.flux_mode, correction, native=True)
    drz_columns = branch_columns(from2d, options.flux_mode, correction)
    x2d_columns = branch_columns(frommast, options.flux_mode, correction)
    metadata = dict(
        EXTRPIX=options.height_native_pixels,
        OPTPIX=options.optimal_height_native_pixels,
        PROFPX=options.profile_height_native_pixels,
        EXTRASEC=options.height_native_pixels * native_scale,
        EXTRCENT=options.center_arcsec,
        SLITWID=widths[0],
    )
    write_spectrum(
        outputs[0],
        from2d,
        dict(PRODTYPE="CUSTOM_2D_APERTURE_1D", **metadata),
        [dict(input=clean.name, sha256=digest(clean))],
        [
            "Finite aperture diffuse flux: SCI times slit width times spatial overlap.",
            "No point-source throughput correction; no multiplication by wavelength width.",
            "ERROR uses conservative spatial covariance bounds; ERR_DIAG is diagonal only.",
        ],
    )
    write_spectrum(
        outputs[1],
        combined,
        dict(PRODTYPE="CUSTOM_CALSTIS_X1D_COADD",
             BACKCORR="OMIT",
             RESBKG="APPLIED" if any(r["residual_applied"] for r in records) else "NOT_APPLIED",
             **metadata),
        records,
        [
            "FLUX uses official point-source sensitivity; shared native sky model, BACKCORR=OMIT.",
            "FLUX_APER uses matched x2d diffuse convention for finite-aperture comparison.",
            "EXPTIME weighted; overlap fractions conserve flux density; squared variance weights.",
            "BKGERR is fully correlated bound; COVAR1_STAT_APER is adjacent-bin statistical covariance.",
        ],
    )
    for p in outputs:
        with fits.open(p, checksum=True) as hd:
            if not all(
                h.verify_checksum() == 1 and h.verify_datasum() == 1 for h in hd
            ):
                raise ValueError("1D FITS checksum failed")
    unchanged = all(digest(e.path) == e.sha for e in exposures) and all(
        digest(products / n) == v for n, v in parent["products"].items()
    )
    if not unchanged:
        raise ValueError("Original data or parent 2D products changed")
    write_consolidated(
        final,
        config,
        {"FLT": directory / (exposures[0].path.stem.replace("_flt", "") + "_bkgsub_x1d.fits"),
         "DRZ": clean, "X2D": mast},
        flt_columns, drz_columns, x2d_columns,
        [directory / (e.path.stem.replace("_flt", "") + "_bkgsub_x1d.fits") for e in exposures],
        records,
    )
    report = dict(
        config=vars(options),
        inputs_preserved=unchanged,
        parent_manifest_sha256=digest(manifest_path),
        products={p.name: digest(p) for p in [*outputs, final]},
        exposures=records,
        code_sha256={p.name: digest(p) for p in Path(__file__).parent.glob("*.py")},
    )
    from .diagnostics import extraction_diagnostics, extraction_branch_plot

    extraction_branch_plot(config, {"FLT":flt_columns,"DRZ":drz_columns,"X2D":x2d_columns})

    report.update(extraction_diagnostics(
        config, from2d, combined, flt_columns, drz_columns, records
    ))
    (spectra_dir / "manifest.json").write_text(json.dumps(report, indent=2))
    return [final]
