"""FITS input validation, count-valued imsets, and source-file checksums."""

from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
import hashlib
import numpy as np
from astropy.io import fits


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class Exposure:
    path: Path
    extver: int
    p: dict
    h: dict
    counts: np.ndarray
    variance: np.ndarray
    dq: np.ndarray
    time: float
    sha: str
    y: float = 0.0
    shift: float = 0.0
    scale: float = 0.0
    wave: np.ndarray | None = None
    trace: np.ndarray | None = None
    background: np.ndarray | None = None

    @property
    def label(self):
        return f"{self.path.stem}[SCI,{self.extver}]"


def read_inputs(paths):
    """Read every SCI imset, using its own EXPTIME and matching ERR/DQ."""
    result = []
    for path in paths:
        sha = digest(path)
        with fits.open(path, memmap=False) as hd:
            multiple = sum(s.name == "SCI" for s in hd) > 1
            for s in (h for h in hd if h.name == "SCI"):
                v = s.header.get("EXTVER", 1)
                e = hd["ERR", v]
                d = hd["DQ", v]
                p = dict(hd[0].header)
                h = dict(s.header)
                if h.get("BUNIT", "").strip().upper() not in ("COUNTS", "COUNT"):
                    raise ValueError(f"{path}[{v}] SCI BUNIT is not COUNTS")
                if e.header.get("BUNIT", "").strip().upper() not in ("COUNTS", "COUNT"):
                    raise ValueError(f"{path}[{v}] ERR BUNIT is not COUNTS")
                if s.data.ndim != 2 or s.shape != e.shape or s.shape != d.shape:
                    raise ValueError("SCI/ERR/DQ shape mismatch")
                if multiple and "EXPTIME" not in h:
                    raise ValueError("Multi-imset SCI must have its own EXPTIME")
                t = float(h.get("EXPTIME", p.get("EXPTIME", np.nan)))
                if not np.isfinite(t) or t <= 0:
                    raise ValueError("Invalid imset exposure time")
                err = e.data.astype(float)
                var = np.where(np.isfinite(err) & (err >= 0), err**2, np.nan)
                result.append(
                    Exposure(
                        path,
                        v,
                        p,
                        h,
                        s.data.astype(float),
                        var,
                        d.data.astype(np.uint16),
                        t,
                        sha,
                    )
                )
    if not result:
        raise ValueError("No SCI imsets found")
    return result


def validate(exposures, config):
    """Reject unsupported configurations before any science products are written."""
    first = exposures[0]
    target = str(first.p.get("TARGNAME", "")).strip()
    config.resolve_name(exposures)
    unknown = set(config.source_offsets) - {e.path.name for e in exposures}
    if unknown:
        raise ValueError(
            f"Registration overrides refer to unknown inputs: {sorted(unknown)}"
        )
    for e in exposures:
        with fits.open(e.path) as hd:
            multiple = sum(x.name == "SCI" for x in hd) != 1
        if e.extver != 1 or multiple:
            raise ValueError(
                "Current coordinate-export workflow supports one imset per FLT"
            )
        if e.p["DETECTOR"] != "FUV-MAMA" or e.p["OPT_ELEM"] != "G140L":
            raise ValueError("Only STIS FUV-MAMA G140L is validated")
        if e.counts.shape != (1024, 1024):
            raise ValueError("Only full-frame unbinned FUV-MAMA is validated")
        for key in ("DARKCORR", "FLATCORR", "WAVECORR"):
            if e.p.get(key) != "COMPLETE":
                raise ValueError(f"{e.path}: {key} must be COMPLETE")
        for key in ("DETECTOR", "OPT_ELEM", "CENWAVE", "APERTURE", "PROPAPER"):
            if e.p.get(key) != first.p.get(key):
                raise ValueError(f"Incompatible {key}")
        if str(e.p.get("TARGNAME", "")).strip() != target:
            raise ValueError("Incompatible TARGNAME")
        angle = (float(e.h["PA_APER"]) - float(first.h["PA_APER"]) + 180) % 360 - 180
        if abs(angle) > 0.01:
            raise ValueError(
                "Slit orientations differ; cannot coadd on one along-slit grid"
            )
        for key in (
            "PHOTTAB",
            "TDSTAB",
            "APDESTAB",
            "SDCTAB",
            "SPTRCTAB",
            "DISPTAB",
            "INANGTAB",
        ):
            if not (
                config.path(config.reference_dir) / e.p[key].split("$")[-1]
            ).is_file():
                raise FileNotFoundError(e.p[key])
        mast = config.path(config.mast_dir) / e.path.name.replace("_flt", "_x2d")
        if not mast.is_file():
            raise FileNotFoundError(mast)


def save_background(path, e, result):
    """Write native-rate intermediates, retaining original detector WCS and DQ."""
    header = fits.getheader(e.path, ("SCI", e.extver)).copy()
    header["BUNIT"] = "count / s"
    ph = fits.getheader(e.path, 0).copy()
    ph["PRODTYPE"] = "CUSTOM_NATIVE_RATE"
    ph["HSTCALOK"] = False
    ph["BKG2DSUB"] = result.metadata["enabled"]
    ph["RESSUB"] = bool(result.metadata.get("residual", {}).get("applied", False))
    ph["RESMETH"] = "UNCLIPPED_VALIDATED"
    ph["BACKCORR"] = "OMIT"
    ph.add_history("BKG=BKGPRIMARY+BKGRESID; residual inferred from disjoint guard strips.")
    ph.add_history("No repeated CALSTIS sky subtraction. BKGERR is a correlated bound.")
    ph["SKYWIDTH"] = result.metadata.get("selected_width", 0.0)
    ph.add_history(
        "Rate=original COUNTS/EXPTIME. Residual sky only; DARKCORR not repeated."
    )
    hd = [fits.PrimaryHDU(header=ph)]
    original_error = np.sqrt(e.variance) / e.time
    for name, data in [
        ("SCI", e.counts / e.time - result.rate),
        ("ERRBASE", original_error),
        ("ERR", np.hypot(original_error, result.error)),
        ("BKG", result.rate),
        ("BKGERR", result.error),
    ]:
        hd.append(fits.ImageHDU(data.astype("f4"), header=header.copy(), name=name))
    for name, data in (result.components or {}).items():
        component_header = header.copy()
        component_header["BUNIT"] = "" if name == "RESDQ" else "count / s"
        hd.append(fits.ImageHDU(data.astype("u2" if name == "RESDQ" else "f4"),
                                header=component_header, name=name))
    header["BUNIT"] = ""
    hd.extend(
        [
            fits.ImageHDU(e.dq, header=header.copy(), name="DQ"),
            fits.ImageHDU(result.flags, header=header.copy(), name="BKGDQ"),
        ]
    )
    fits.HDUList(hd).writeto(path, overwrite=True, checksum=True)


def save_rectified_background(path, original, result):
    """Save a custom cleaned x2d copy with original WCS, DQ and ERRBASE.

    Input and model are both surface brightness. Original dark subtraction,
    rectification and flux calibration must not be repeated.
    """
    with fits.open(original) as hd:
        h = hd['SCI'].header.copy()
        sci, err = hd['SCI'].data.astype(float), hd['ERR'].data.astype(float)
        hd['SCI'].data = (sci - result.rate).astype('f4')
        hd['ERR'].data = np.hypot(err, result.error).astype('f4')
        hd[0].header['ORIGBACK'] = hd[0].header.get('BACKCORR', 'UNKNOWN')
        hd[0].header['BACKCORR'] = 'OMIT'
        hd[0].header['BKG2DSUB'] = result.metadata['enabled']
        hd[0].header['OFFICIAL'] = False
        hd[0].header['FILENAME'] = path.name
        hd[0].header['PRODTYPE'] = 'CUSTOM_X2D_BKGSUB'
        hd[0].header.add_history('Unclipped 2D sky model in existing surface-brightness units.')
        hd[0].header.add_history('Original dark/flux calibration and WCS retained; no new resampling.')
        hd[0].header.add_history('ERRBASE retains pipeline ERR; BKGERR is a correlated model bound.')
        for name, data in [('ERRBASE',err),('BKG',result.rate),('BKGERR',result.error)]:
            hd.append(fits.ImageHDU(data.astype('f4'), header=h.copy(), name=name))
        h['BUNIT'] = ''
        hd.append(fits.ImageHDU(result.flags, header=h, name='BKGDQ'))
        hd.writeto(path, overwrite=False, checksum=True)
