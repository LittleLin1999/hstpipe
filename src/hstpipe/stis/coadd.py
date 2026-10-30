"""Fixed-weight surface-brightness drizzle for native FLT and rectified MAST data."""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from .calibration import UNIT
from .resampling import drizzle


class Coadd:
    """Accumulate SCI once; propagate photon variance and correlated model bounds.

    Input arrays are already calibrated surface brightness. EXPTIME is a weight,
    not a second rate conversion. Drizzle's secondary channel uses squared
    redistribution weights. Model-error bounds use linear positive weights.
    """

    def __init__(self, grid, pixfrac=1.0):
        from drizzle.resample import Drizzle

        self.grid = grid
        self.pixfrac = pixfrac
        self.norm = 1e-15
        self.science = Drizzle(out_shape=grid.shape, fillval="NaN", fillval2="NaN")
        self.model = Drizzle(out_shape=grid.shape, fillval="NaN")
        self.unc_num = np.zeros(grid.shape)
        self.expmap = np.zeros(grid.shape)
        self.nexp = np.zeros(grid.shape, "u2")
        self.bkgdq = np.zeros(grid.shape, "u2")
        self.rejected_dq = np.zeros(grid.shape, "u2")
        self.provenance = []

    def add(
        self,
        science,
        variance,
        mapping,
        valid,
        covered,
        exptime,
        record,
        model=None,
        model_error=None,
        model_flags=None,
        input_dq=None,
    ):
        weight = (valid * exptime).astype("f4")
        self.science.add_image(
            np.where(valid, science / self.norm, 0).astype("f4"),
            exptime=1,
            pixmap=mapping,
            data2=np.where(valid, variance / self.norm**2, 0).astype("f4"),
            weight_map=weight,
            in_units="cps",
            iscale=1,
            pixel_scale_ratio=1,
            pixfrac=self.pixfrac,
        )
        one = drizzle(
            np.ones(valid.shape),
            np.zeros(valid.shape),
            weight,
            mapping,
            self.grid.shape,
            self.pixfrac,
        )
        full = drizzle(
            np.ones(valid.shape),
            np.zeros(valid.shape),
            covered.astype("f4"),
            mapping,
            self.grid.shape,
            self.pixfrac,
        )
        self.expmap += np.divide(
            one.out_wht,
            full.out_wht,
            out=np.zeros(self.grid.shape),
            where=full.out_wht > 0,
        )
        self.nexp += one.out_wht > 0
        if model is not None:
            self.model.add_image(
                np.where(valid, model / self.norm, 0).astype("f4"),
                exptime=1,
                pixmap=mapping,
                weight_map=weight,
                in_units="cps",
                iscale=1,
                pixel_scale_ratio=1,
                pixfrac=self.pixfrac,
            )
            uncertainty = drizzle(
                np.where(valid, model_error / self.norm, 0),
                np.zeros(valid.shape),
                weight,
                mapping,
                self.grid.shape,
                self.pixfrac,
            )
            self.unc_num += (
                np.nan_to_num(uncertainty.out_img.astype(float))
                * uncertainty.out_wht
                * self.norm
            )
        for flags, destination, domain in [
            (model_flags, self.bkgdq, valid),
            (input_dq, self.rejected_dq, covered),
        ]:
            if flags is None:
                continue
            present = int(np.bitwise_or.reduce(flags[domain])) if domain.any() else 0
            for k in range(16):
                bit = 1 << k
                if not present & bit:
                    continue
                projection = drizzle(
                    np.ones(valid.shape),
                    np.zeros(valid.shape),
                    (domain & ((flags & bit) != 0)).astype("f4"),
                    mapping,
                    self.grid.shape,
                    self.pixfrac,
                )
                destination[projection.out_wht > 0] |= bit
        self.provenance.append(record)

    def arrays(self, background=False):
        w = self.science.out_wht.astype(float)
        valid = w > 0
        sci = self.science.out_img.astype(float) * self.norm
        variance = self.science.out_img2[0].astype(float) * self.norm**2
        base = np.sqrt(variance)
        bound = np.divide(self.unc_num, w, out=np.zeros_like(w), where=valid)
        sci[~valid] = np.nan
        base[~valid] = np.nan
        bound[~valid] = np.nan
        result = dict(
            SCI=sci,
            ERR=np.hypot(base, bound),
            DQ=np.where(valid, 0, 4).astype("u2"),
            EXPMAP=self.expmap,
            NEXP=self.nexp,
            WHT=w,
            REJDQ=self.rejected_dq,
        )
        if background:
            result.update(
                BKG=self.model.out_img.astype(float) * self.norm,
                BKGERR=bound,
                ERRBASE=base,
                BKGDQ=self.bkgdq,
            )
        return result

    def write(self, path, config, kind, background=False, arrays=None, metadata=None):
        ph = fits.Header()
        ph["PRODTYPE"] = kind
        ph["FILENAME"] = path.name
        ph["OFFICIAL"] = False
        ph["TARGNAME"], ph["OPT_ELEM"] = config.name.rsplit("_", 1)
        ph["BACKCORR"] = "OMIT"
        ph["RESSUB"] = bool(background and any(r.get("residual_applied", False) for r in self.provenance))
        ph["RESMETH"] = "UNCLIPPED_VALIDATED" if background else "NONE"
        bg_options = config.mast_background if kind == "CUSTOM_MAST_X2D_COADD" else config.background
        ph["BKG2DSUB"] = bool(background and bg_options.enabled)
        ph["PIXFRAC"] = self.pixfrac
        ph["TOTEXPT"] = sum(r["exptime"] for r in self.provenance)
        ph["ISCALDRZ"] = 1.0
        ph["REGWARN"] = True
        ph.add_history(
            "Header-based target registration; source-centroid corrections not inferred."
        )
        ph.add_history(
            "SCI is diffuse surface brightness; no DIFF2PT or point-source correction."
        )
        ph.add_history(
            "EXPTIME weights only; ERR propagates input variances with squared weights."
        )
        ph.add_history(
            "Resampling covariance is not contained in diagonal ERR; inspect diagnostics."
        )
        if background:
            ph["BKGERRTP"] = "CORR_UPPER_BOUND"
            ph.add_history(
                "Individual x2d sky subtracted in surface brightness before combination."
                if kind == "CUSTOM_MAST_X2D_COADD" else
                "Native count-rate background subtracted before calibration and one SCI drizzle."
            )
            ph.add_history(
                "Background uncertainty is spatially correlated; photon noise counted once."
            )
        else:
            ph.add_history(
                "Official MAST x2d inputs aligned/coadded; custom combination, not CALSTIS output."
            )
            ph.add_history(
                "ERR additionally omits pre-existing x2d interpolation covariance."
            )
        if kind == "CUSTOM_MAST_X2D_COADD":
            ph.add_history("Input x2d interpolation covariance is not fully known.")
            ph.add_history("Coadd residual is additional only; no second primary-sky subtraction.")
        for key, value in (metadata or {}).items():
            ph[key] = value
        hd = [fits.PrimaryHDU(header=ph)]
        for name, array in (self.arrays(background) if arrays is None else arrays).items():
            h = self.grid.header()
            h["BUNIT"] = (
                ""
                if name in ("DQ", "REJDQ", "NEXP", "BKGDQ", "RESDQ")
                else ("s" if name in ("EXPMAP", "WHT") else UNIT)
            )
            h.add_history(
                "REJDQ records rejected-input flags; DQ=4 marks no usable coverage."
            )
            dtype = "u2" if array.dtype.kind in "ui" else "f4"
            hd.append(fits.ImageHDU(array.astype(dtype), header=h, name=name))
        hd.append(fits.BinTableHDU(Table(rows=self.provenance), name="PROVENANCE"))
        fits.HDUList(hd).writeto(path, overwrite=False, checksum=True)
        return path
