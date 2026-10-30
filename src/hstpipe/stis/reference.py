"""Official x2d controls and coordinate exports, kept separate from SCI drizzle.

Requires a coordinate-export-enabled HSTCAL installation when controls are not
already supplied. This module never patches installed C code silently.
"""

import os
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates
from .geometry import invert_mapping
from .calibration import factor


def prepare(exposures, config):
    directory = config.path(config.control_dir)
    directory.mkdir(parents=True, exist_ok=True)
    for e in exposures:
        required = [
            directory / (e.path.stem + suffix)
            for suffix in ("_omit.fits", "_perform.fits", "_mapping.f64")
        ]
        if all(p.exists() for p in required):
            continue
        if not config.prepare_controls:
            raise FileNotFoundError(
                f"Missing official controls for {e.path.name}; provide control_dir or enable prepare_controls with coordinate-export HSTCAL"
            )
        from stistools import x2d

        old = {
            key: os.environ.get(key) for key in ("oref", "PATH", "STIS_COORD_EXPORT")
        }
        try:
            os.environ["oref"] = str(config.path(config.reference_dir)) + "/"
            if config.software_dir:
                os.environ["PATH"] = (
                    str(config.path(config.software_dir) / "hstcal-install/bin")
                    + os.pathsep
                    + os.environ["PATH"]
                )
            os.environ.pop("STIS_COORD_EXPORT", None)
            for mode, dest in zip(("omit", "perform"), required[:2]):
                if not dest.exists():
                    status = x2d.x2d(
                        str(e.path),
                        output=str(dest),
                        fluxcorr=mode,
                        trailer=str(dest.with_suffix(".log")),
                    )
                    if status:
                        raise RuntimeError(f"x2d failed: {dest}")
            if not required[2].exists():
                os.environ["STIS_COORD_EXPORT"] = str(required[2])
                temp = directory / (e.path.stem + "_coordinate_check.fits")
                if temp.exists():
                    raise FileExistsError(temp)
                status = x2d.x2d(
                    str(e.path),
                    output=str(temp),
                    fluxcorr="omit",
                    trailer=str(temp.with_suffix(".log")),
                )
                if status or not required[2].exists():
                    raise RuntimeError(
                        "HSTCAL did not export coordinates; configure the instrumented executable"
                    )
                np.testing.assert_array_equal(
                    fits.getdata(temp, "SCI"), fits.getdata(required[0], "SCI")
                )
                temp.unlink()
        finally:
            for key, value in old.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


def operator(exposure, config):
    """Validate official controls; map original native pixels without resampling SCI."""
    e = exposure
    directory = config.path(config.control_dir)
    prefix = e.path.stem
    with fits.open(directory / (prefix + "_perform.fits")) as hd:
        h = hd["SCI"].header.copy()
        official = hd["SCI"].data.astype(float)
        good = hd["DQ"].data == 0
    omit = fits.getdata(directory / (prefix + "_omit.fits"), "SCI").astype(float)
    field = np.fromfile(directory / (prefix + "_mapping.f64"), dtype="f8").reshape(
        official.shape + (2,)
    )
    interpolated = map_coordinates(
        e.counts,
        [field[:, :, 1], field[:, :, 0]],
        order=1,
        mode="constant",
        cval=0,
        prefilter=False,
    )
    g = good & (omit != 0) & np.isfinite(official)
    if not g.any():
        raise ValueError("No usable official calibration control pixels")
    geometry_error = float(
        np.max(abs(interpolated[g] - omit[g])) / max(np.max(abs(omit[g])), 1.0)
    )
    if geometry_error > 2e-6:
        raise ValueError(
            "Coordinate map does not reproduce original exposure rectification"
        )
    wave = h["CRVAL1"] + (np.arange(official.shape[1]) + 1 - h["CRPIX1"]) * h["CD1_1"]
    f, meta = factor(
        e,
        Path(config.project_dir),
        wave,
        h,
        reference_dir=config.path(config.reference_dir),
    )
    factor_error = float(
        np.max(abs(official[g] / omit[g] / np.broadcast_to(f, official.shape)[g] - 1))
    )
    if factor_error > 3e-6:
        raise ValueError("Diffuse calibration differs from official x2d")
    native, covered, residual = invert_mapping(field, e.counts.shape)
    native_scale = h["CD2_2"] * 3600
    offset = float(config.source_offsets.get(e.path.name, 0.0))
    spatial = (
        (native[:, :, 1] + 1 - h["CRPIX2"]) * native_scale - e.p["POSTARG2"] - offset
    )
    wavelength = h["CRVAL1"] + (native[:, :, 0] + 1 - h["CRPIX1"]) * h["CD1_1"]
    conversion, _ = factor(
        e,
        Path(config.project_dir),
        wavelength,
        h,
        reference_dir=config.path(config.reference_dir),
    )
    valid = covered & (e.dq == 0) & np.isfinite(e.counts)
    if np.any(valid & ~np.isfinite(e.variance)):
        raise ValueError("Unknown pipeline variance in usable SCI pixels")
    return dict(
        header=h,
        wavelength=wavelength,
        observed_wavelength=wavelength / (1 - h["V_HELIO"] / 299792.458),
        spatial=spatial,
        source_y=spatial / native_scale,
        factor=conversion,
        valid=valid,
        covered=covered,
        meta=dict(
            **meta,
            geometry_error=geometry_error,
            factor_error=factor_error,
            inverse_residual=residual,
        ),
    )
