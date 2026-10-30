"""Invert the coordinate-only export from official HSTCAL X2DCorr7."""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree


def invert_mapping(field, shape):
    """Invert exported cs7 output-to-detector map with Newton iteration.

    Input-pixel centers are mapped directly onto the final grid. Science is
    never interpolated here. Out-of-domain inverse solutions are excluded.
    """
    ny, nx = field.shape[:2]
    yy, xx = np.indices(shape, dtype=float)
    cy, cx = np.mgrid[0:ny:8, 0:nx:8]
    tree = cKDTree(field[::8, ::8].reshape(-1, 2))
    _, k = tree.query(np.stack((xx.ravel(), yy.ravel()), axis=-1))
    u = cx.ravel()[k].astype(float)
    v = cy.ravel()[k].astype(float)
    fs = [
        RegularGridInterpolator(
            (np.arange(ny), np.arange(nx)),
            field[:, :, i],
            bounds_error=False,
            fill_value=None,
        )
        for i in (0, 1)
    ]
    grads = []
    for i in (0, 1):
        dy, dx = np.gradient(field[:, :, i])
        grads.extend(
            [
                RegularGridInterpolator(
                    (np.arange(ny), np.arange(nx)),
                    a,
                    bounds_error=False,
                    fill_value=None,
                )
                for a in (dx, dy)
            ]
        )
    for _ in range(8):
        points = np.stack((v, u), axis=-1)
        rx = fs[0](points) - xx.ravel()
        ry = fs[1](points) - yy.ravel()
        a, b, c, d = [g(points) for g in grads]
        det = a * d - b * c
        u -= (d * rx - b * ry) / det
        v -= (-c * rx + a * ry) / det
    points = np.stack((v, u), axis=-1)
    res = np.hypot(fs[0](points) - xx.ravel(), fs[1](points) - yy.ravel())
    valid = (res < 0.005) & (u >= 1) & (u < nx - 2) & (v >= 1) & (v < ny - 2)
    return (
        np.stack((u.reshape(shape), v.reshape(shape)), axis=-1),
        valid.reshape(shape),
        float(np.max(res[valid])),
    )


from dataclasses import dataclass
from astropy.io import fits


@dataclass
class Grid:
    """Linear wavelength and target-relative spatial output coordinates."""

    wave_start: float
    wave_step: float
    nx: int
    spatial_start: float
    spatial_step: float
    ny: int

    @property
    def shape(self):
        return (self.ny, self.nx)

    def map(self, wavelength, spatial_arcsec):
        return np.stack(
            (
                (wavelength - self.wave_start) / self.wave_step,
                (spatial_arcsec - self.spatial_start) / self.spatial_step,
            ),
            axis=-1,
        )

    def header(self):
        h = fits.Header()
        for k, v in dict(
            CTYPE1="WAVE",
            CUNIT1="Angstrom",
            CRPIX1=1.0,
            CRVAL1=self.wave_start,
            CD1_1=self.wave_step,
            CD1_2=0.0,
            CTYPE2="LINEAR",
            CUNIT2="arcsec",
            CRPIX2=1.0,
            CRVAL2=self.spatial_start,
            CD2_1=0.0,
            CD2_2=self.spatial_step,
        ).items():
            h[k] = v
        h["CNAME2"] = "Along slit relative to commanded target"
        return h


def rectified_coordinates(header, shape, postarg2, extra_offset=0.0):
    """Wavelength and target-relative Y from an official rectified SCI header.

    SHIFTA values have already entered x2d's rectification; do not apply twice.
    extra_offset is a documented registration correction in arcseconds.
    """
    if abs(header.get("CD1_2", 0)) > 1e-10 or abs(header.get("CD2_1", 0)) > 1e-10:
        raise ValueError("Nonseparable rectified WCS is unsupported")
    scale = float(header["CD2_2"]) * 3600
    if scale <= 0 or header["CD1_1"] <= 0:
        raise ValueError("Reversed rectified axes unsupported")
    yy, xx = np.indices(shape, dtype=float)
    wave = header["CRVAL1"] + (xx + 1 - header["CRPIX1"]) * header["CD1_1"]
    spatial = (yy + 1 - header["CRPIX2"]) * scale - postarg2 - extra_offset
    return wave, spatial


def make_grid(headers, config, spatial_offsets=None):
    """Use header coverage, preserving native spectral sampling by default."""
    step = float(headers[0]["CD1_1"])
    if any(not np.isclose(h["CD1_1"], step, rtol=1e-6) for h in headers):
        raise ValueError(
            "Different spectral dispersions require an explicit regridding design"
        )
    starts = [h["CRVAL1"] + (1 - h["CRPIX1"]) * h["CD1_1"] for h in headers]
    ends = [s + (h["NAXIS1"] - 1) * step for s, h in zip(starts, headers)]
    start = min(starts)
    nx = int(np.ceil((max(ends) - start) / step - 1e-7)) + 1
    if config.spatial_limits is None:
        # Enclose all target-relative footprints, including commanded dithers.
        # A symmetric nominal-slit extent would crop large POSTARG2 offsets.
        offsets = [0.0] * len(headers) if spatial_offsets is None else spatial_offsets
        if len(offsets) != len(headers):
            raise ValueError("One spatial offset is required for each header")
        lower = min(
            (1 - h["CRPIX2"]) * h["CD2_2"] * 3600 - offset
            for h, offset in zip(headers, offsets)
        )
        upper = max(
            (h["NAXIS2"] - h["CRPIX2"]) * h["CD2_2"] * 3600 - offset
            for h, offset in zip(headers, offsets)
        )
        lower = np.floor(lower / config.spatial_scale) * config.spatial_scale
        upper = np.ceil(upper / config.spatial_scale) * config.spatial_scale
    else:
        lower, upper = config.spatial_limits
    ny = int(np.ceil((upper - lower) / config.spatial_scale - 1e-7)) + 1
    return Grid(start, step, nx, float(lower), config.spatial_scale, ny)
