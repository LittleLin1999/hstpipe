"""Assemble the three comparable 1D branches into one documented FITS file."""

from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from .io import digest

FLUX_UNIT = "erg / (s cm2 Angstrom)"


def point_conversion(exposures, times, output_wave, overlap_matrix, ratio_key="CAL_RATIO"):
    """Transfer CALSTIS point/diffuse response ratios to the common grid.

    The ratio comes from positive unit-SCI controls, never from noisy source
    fluxes. It includes CALSTIS's aperture and slit corrections. The 2D
    aperture geometry is still an approximation to the native x1d aperture.
    """
    numerator = np.zeros(len(output_wave))
    denominator = np.zeros(len(output_wave))
    for spectrum, time in zip(exposures, times):
        matrix = overlap_matrix(spectrum["WAVELENGTH"], output_wave)
        ratio = np.asarray(spectrum[ratio_key], float)
        # Calibration validity is independent of bad pixels in the science aperture.
        caldq = spectrum.get("CALDQ", np.zeros(len(ratio), dtype="u2"))
        good = (caldq == 0) & np.isfinite(ratio) & (ratio > 0)
        matrix[:, ~good] = 0
        complete = matrix.sum(axis=1) >= 1 - 1e-6
        matrix[~complete] = 0
        inverse = np.divide(1, ratio, out=np.zeros_like(ratio), where=good)
        numerator += time * (matrix @ inverse)
        denominator += time * complete
    return np.divide(
        numerator, denominator,
        out=np.full(len(output_wave), np.nan), where=denominator > 0,
    )


def morphology_conversion(path, output_wave, fraction_key="APFRAC"):
    """Read a user-supplied model's aperture fraction versus wavelength.

    The selected fraction must include both slit and extraction-height losses. A
    morphology correction is therefore explicit and never inferred from the
    point-source calibration or the observed background-subtracted image.
    """
    with fits.open(path) as hd:
        table = next((h.data for h in hd if isinstance(h, fits.BinTableHDU)
                      and {"WAVELENGTH", fraction_key} <= set(h.columns.names)), None)
        if table is None:
            raise ValueError(f"Source model needs WAVELENGTH and {fraction_key} columns")
        wave = np.asarray(table["WAVELENGTH"], float)
        fraction = np.asarray(table[fraction_key], float)
    if (len(wave) < 2 or np.any(np.diff(wave) <= 0)
            or not np.isfinite(fraction).all()
            or np.any((fraction <= 0) | (fraction > 1))):
        raise ValueError("Source-model aperture fraction must be finite and in (0, 1]")
    correction = np.interp(output_wave, wave, 1 / fraction, left=np.nan, right=np.nan)
    return correction


def branch_columns(data, mode, correction, native=False):
    """Select boxcar and empirical-profile measurements in one flux convention."""
    if native:
        if mode == "point":
            names = ("FLUX", "ERROR", "FLUX_OPT", "ERR_OPT", "BKGERR", "BKGERR_OPT")
            scale = 1.0
        else:
            names = ("FLUX_APER", "ERROR_APER", "FLUX_OPT_APER", "ERR_OPT_APER",
                     "BKGERR_APER", "BKGERR_OPT_APER")
            scale = correction if mode == "morph" else 1.0
    else:
        names = ("FLUX", "ERROR", "FLUX_OPT", "ERR_OPT", "BKGERR", "BKGERR_OPT")
        scale = correction if mode != "diffuse" else 1.0
    box_scale = scale.get('BOX',1.) if isinstance(scale,dict) else scale
    opt_scale = scale.get('OPT',1.) if isinstance(scale,dict) else scale
    values = {
        "WAVELENGTH": np.asarray(data["WAVELENGTH"], float),
        "FLUX_BOX": np.asarray(data[names[0]], float) * box_scale,
        "ERR_BOX": np.asarray(data[names[1]], float) * box_scale,
        "FLUX_OPT": np.asarray(data[names[2]], float) * opt_scale,
        "ERR_OPT": np.asarray(data[names[3]], float) * opt_scale,
        "BKGERR_BOX": np.asarray(data[names[4]], float) * box_scale,
        "BKGERR_OPT": np.asarray(data[names[5]], float) * opt_scale,
        "EXPTIME": np.asarray(data["EXPTIME"], float),
        "NEXP": np.asarray(data["NEXP"], "u2"),
        "DQ": np.asarray(data["DQ"], "u2").copy(),
    }
    for method in ("BOX", "OPT"):
        dq = np.asarray(data.get("DQ_"+method, data["DQ"]), "u2").copy()
        bad = ~np.isfinite(values["FLUX_"+method]) | ~np.isfinite(values["ERR_"+method])
        dq[bad] |= 4
        values["DQ_"+method] = dq
        # A failed method must never invalidate the other method.
        for key in ("FLUX_", "ERR_", "BKGERR_"):
            values[key+method] = np.where(dq == 0, values[key+method], np.nan)
        for key in ("EXPTIME", "NEXP"):
            values[key+"_"+method] = np.asarray(data.get(key+"_"+method,data[key]), dtype="u2" if key=="NEXP" else float)
    values["DQ"] = values["DQ_BOX"].copy()  # backwards-compatible boxcar alias
    for key in ("PROFCOV_OPT", "INFOFRAC_OPT", "OPT_MASKED", "COVERAGE"):
        if key in data:
            values[key] = np.asarray(data[key])
    return values


def write_consolidated(path, config, source_paths, flt, drz, x2d,
                       native_paths, residual_records):
    """Write exactly FLT, DRZ and X2D BINTABLE extensions with FITS checksums."""
    mode = config.extraction.flux_mode
    ph = fits.Header()
    ph["PRODTYPE"] = "CUSTOM_X1D_COMPARISON"
    ph["OFFICIAL"] = False
    ph["FILENAME"] = Path(path).name
    ph["TARGNAME"], ph["OPT_ELEM"] = config.name.rsplit("_", 1)
    ph["FLXMODE"] = mode.upper()
    ph["EXTRPIX"] = config.extraction.height_native_pixels
    ph['OPTPIX'] = config.extraction.optimal_height_native_pixels
    ph['PROFPX'] = config.extraction.profile_height_native_pixels
    ph['CENTMETH'] = config.extraction.center_method
    ph["EXTRCENT"] = config.extraction.center_arcsec
    ph["OPTMETH"] = "EMPIRICAL_EQUALVAR"
    ph["NINPUT"] = len(native_paths)
    ph["MINPROF"] = config.extraction.optimal_min_profile_fraction
    ph["MININFO"] = config.extraction.optimal_min_information_fraction
    ph.add_history("DQ/EXPTIME/NEXP alias BOX; use method-specific *_BOX and *_OPT columns.")
    ph.add_history("Optimal model normalized over full aperture before masking; no gap interpolation.")
    applied = [bool(record["residual_applied"]) for record in residual_records]
    if len(applied) != len(native_paths):
        raise ValueError("One residual-background decision is required per FLT")
    ph["RESBKG"] = "APPLIED" if all(applied) else "NOT_APPLIED" if not any(applied) else "MIXED"
    if mode == "morph":
        ph["SRCMODEL"] = Path(config.extraction.source_model).name
        ph.add_history("Morphology corrections from supplied APFRAC_BOX and APFRAC_OPT model.")
    ph.add_history("FLT: official CALSTIS x1d boxcar plus custom empirical-profile extraction.")
    ph.add_history("FLT/DRZ share primary plus validated native residual sky; x1d BACKCORR=OMIT.")
    ph.add_history("No additional CALSTIS background subtraction is applied.")
    ph.add_history("DRZ: extraction of native-FLT single-pass diffuse-calibrated drizzle.")
    ph.add_history("X2D: individually sky-subtracted MAST x2d, then validated coadd residual.")
    ph.add_history("Profile weights avoid inverse near-zero MAMA errors; negative SCI retained.")
    ph.add_history("2D ERR values bound unknown drizzle covariance; model errors correlated.")
    for label, source in source_paths.items():
        ph.add_history(f"{label} input {Path(source).name} SHA256 {digest(source)}")
    for source in native_paths:
        ph.add_history(f"FLT exposure {Path(source).name} SHA256 {digest(source)}")
    for record in residual_records:
        audit = record["residual_audit"]
        ph.add_history(
            f"{record['input']}: residual applied={record['residual_applied']}; "
            f"method={audit.get('method', 'legacy')}"
        )
    hdus = [fits.PrimaryHDU(header=ph)]
    for name, columns in (("FLT", flt), ("DRZ", drz), ("X2D", x2d)):
        table = Table(columns)
        table["WAVELENGTH"].unit = "Angstrom"
        for key in ("EXPTIME", "EXPTIME_BOX", "EXPTIME_OPT"):
            table[key].unit = "s"
        for key in ("FLUX_BOX", "ERR_BOX", "FLUX_OPT", "ERR_OPT",
                    "BKGERR_BOX", "BKGERR_OPT"):
            table[key].unit = FLUX_UNIT
        hdu = fits.BinTableHDU(table, name=name)
        hdu.header["FLXMODE"] = mode.upper()
        hdu.header["APHEIGHT"] = config.extraction.height_native_pixels
        hdu.header['OPTHEIGH'] = config.extraction.optimal_height_native_pixels
        hdu.header['PROFHEIG'] = config.extraction.profile_height_native_pixels
        hdu.header['EXTRCENT'] = config.extraction.center_arcsec
        hdu.header.add_history('Boxcar and optimal have different apertures and official point responses.')
        hdu.header["OPTMETH"] = "EMPIRICAL_EQUALVAR"
        hdu.header["DQALIAS"] = "DQ_BOX"
        hdu.header.add_history("PROFCOV_OPT and INFOFRAC_OPT quantify observed profile support.")
        hdu.header["BACKCORR"] = "OMIT"
        if name == "FLT":
            hdu.header["RESBKG"] = ph["RESBKG"]
        elif name == "X2D":
            source_header = fits.getheader(source_paths[name])
            hdu.header['BKG2DSUB'] = source_header.get('BKG2DSUB', False)
            hdu.header['COADDSUB'] = source_header.get('COADDSUB', False)
            hdu.header['BKGSTAGE'] = source_header.get('BKGSTAGE', 'UNKNOWN')
        hdu.header["SRCFILE"] = Path(source_paths[name]).name
        hdu.header["SRCSHA"] = digest(source_paths[name])
        hdu.header.add_history("ERR_OPT includes a conservative background-model bound.")
        if name != "FLT" and mode == "point":
            hdu.header.add_history("Point response transferred from official per-FLT unit-SCI controls.")
        hdus.append(hdu)
    fits.HDUList(hdus).writeto(path, overwrite=False, checksum=True)
    with fits.open(path, checksum=True) as hd:
        assert [h.name for h in hd] == ["PRIMARY", "FLT", "DRZ", "X2D"]
        for h in hd:
            if h.verify_checksum() != 1 or h.verify_datasum() != 1:
                raise ValueError("Consolidated 1D FITS checksum failed")
