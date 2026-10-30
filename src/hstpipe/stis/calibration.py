"""HSTCAL cs7 diffuse-source calibration, without point-source corrections."""

import numpy as np
from astropy.io import fits

UNIT = "erg / (s cm2 Angstrom arcsec2)"


def factor(e, root, wave, h, reference_dir=None):
    """HSTCAL AbsFlux multiplier per detector count at native reference sampling."""

    def table(key):
        directory = (
            reference_dir if reference_dir is not None else root / "calibration_refs"
        )
        return fits.getdata(directory / e.p[key].split("$")[-1], 1)

    t = table("PHOTTAB")
    row = t[
        (t["OPT_ELEM"] == e.p["OPT_ELEM"])
        & (t["CENWAVE"] == e.p["CENWAVE"])
        & (t["SPORDER"] == 1)
    ][0]
    if row["REFWAV"] != 0:
        raise ValueError("Blaze-corrected modes unsupported")
    n = row["NELEM"]
    observed = wave / (1 - h["V_HELIO"] / 299792.458)
    response = np.interp(observed, row["WAVELENGTH"][:n], row["THROUGHPUT"][:n])
    td = table("TDSTAB")
    tr = td[td["OPT_ELEM"] == e.p["OPT_ELEM"]][0]
    nw = int(tr["NWL"])
    nt = int(tr["NT"])
    if "INTERCEPT" in td.names:
        raise ValueError("COS-format TDS unsupported; do not silently approximate")
    times = tr["TIME"][:nt]
    slopes = tr["SLOPE"][:nt, :nw]
    sensitivity = np.ones(nw)
    time = e.h["EXPSTART"]
    for j in range(nt):
        stop = time if j == nt - 1 or time <= times[j + 1] else times[j + 1]
        sensitivity += (stop - times[j]) * slopes[j] / 365.25 / 100
        if stop == time:
            break
    temp = e.h.get("OM1CAT", -1.0)
    if tr["REFTEMP"] > 0 and temp > 0:
        sensitivity *= 1 + tr["TEMPSENS"][:nw] * (temp - tr["REFTEMP"])
    tds = np.interp(observed, tr["WAVELENGTH"][:nw], sensitivity)
    ap = table("APDESTAB")
    ap = ap[ap["APERTURE"] == e.p["PROPAPER"]][0]
    width = float(ap["WIDTH1"])
    # Same native output photometric convention as official x2d: jacobian=1
    # in cs7. The smaller drizzle pixels carry brightness, not pixel flux.
    dispersion = h["CD1_1"]
    spatial = h["CD2_2"] * 3600
    phot = (
        6.6260755e-27
        * 29979245800
        / (45238.93416 * e.time * dispersion * spatial * width)
    )
    correction = (phot / (response * observed * 1e-8)).astype("f4")
    correction = (correction.astype(float) / tds).astype("f4")
    if np.any(~np.isfinite(correction)) or np.any(correction <= 0):
        raise ValueError("Invalid diffuse calibration")
    return correction, dict(
        exptime=e.time,
        temperature=temp,
        slit_width=width,
        native_spatial_scale=spatial,
        dispersion=dispersion,
        tds_min=float(tds.min()),
        tds_max=float(tds.max()),
    )
