"""Independent aperture, spectral-overlap, and count-valued-input checks."""

import tempfile
from pathlib import Path
import unittest
import numpy as np
from astropy.io import fits
from hstpipe.stis.extraction import (
    aperture_weights,
    overlap_matrix,
    extract_image,
    combine_spectra,
    count_copy,
    calibrated_bins,
    write_spectrum,
)
from hstpipe.stis.config import ExtractionConfig


class ExtractionTests(unittest.TestCase):
    def test_nested_audit_provenance_is_serialized(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "audit.fits"
            write_spectrum(
                path,
                {"WAVELENGTH": np.array([1400.]), "FLUX": np.array([1e-15])},
                {},
                [{"input": "test", "residual_audit": {"difference_sigma": 8.0}}],
            )
            with fits.open(path) as hd:
                self.assertIn('"difference_sigma": 8.0',
                              hd["PROVENANCE"].data["residual_audit"][0])

    def test_calibration_support_includes_full_bin(self):
        np.testing.assert_array_equal(
            calibrated_bins(np.arange(5.0), (0.25, 3.5)),
            [False, True, True, True, False],
        )

    def test_fractional_aperture_and_units(self):
        y = np.arange(-20, 21) * 0.0123
        shape = (len(y), 5)
        data = np.full(shape, -2e-15)
        err = np.full(shape, 1e-16)
        bkgerr = np.full(shape, 2e-17)
        a = extract_image(
            data, err, bkgerr, np.ones(shape, bool), y, 0.0123, 0.5, 0.003, 11 * 0.0246
        )
        np.testing.assert_allclose(a["FLUX"], -2e-15 * 0.5 * 11 * 0.0246)
        np.testing.assert_allclose(a["BKGERR"], 2e-17 * 0.5 * 11 * 0.0246)
        self.assertTrue((a["ERR_STAT_BOUND"] > a["ERR_DIAG"]).all())

    def test_no_area_double_scaling(self):
        results = []
        for step in (0.0246, 0.0123):
            y = np.arange(-20, 21) * step
            shape = (len(y), 3)
            result = extract_image(
                np.ones(shape),
                np.ones(shape),
                np.zeros(shape),
                np.ones(shape, bool),
                y,
                step,
                0.5,
                0,
                0.2706,
            )
            results.append(result["FLUX"])
        np.testing.assert_allclose(results[0], results[1])

    def test_missing_pixel_not_renormalized(self):
        y = np.arange(9, dtype=float)
        shape = (9, 3)
        valid = np.ones(shape, bool)
        valid[4, 1] = False
        result = extract_image(
            np.ones(shape), np.zeros(shape), np.zeros(shape), valid, y, 1, 0.5, 4, 3
        )
        self.assertTrue(np.isnan(result["FLUX"][1]))
        self.assertEqual(result["DQ"][1], 4)
        self.assertAlmostEqual(result["COVERAGE"][1], 2 / 3)

    def test_line_flux_conservation_rebin(self):
        wave = np.arange(0.5, 10, 1)
        output = np.arange(0.25, 10, 0.5)
        flux = np.exp(-0.5 * ((wave - 5) / 1.2) ** 2)
        rebinned = overlap_matrix(wave, output) @ flux
        self.assertAlmostEqual(
            float(rebinned.sum() * 0.5), float(flux.sum()), places=12
        )

    def test_unequal_times_variance_covariance_and_masks(self):
        wave = np.arange(0.5, 6, 1)
        out = np.arange(0.25, 6, 0.5)
        spectra = []
        for value in (-2.0, 4.0):
            row = dict(
                WAVELENGTH=wave,
                FLUX=np.full(6, value),
                FLUX_APER=np.full(6, value),
                FLUX_OPT=np.full(6, value),
                FLUX_OPT_APER=np.full(6, value),
                ERR_STAT=np.ones(6),
                ERR_STAT_APER=np.ones(6),
                ERR_STAT_OPT=np.ones(6),
                ERR_STAT_OPT_APER=np.ones(6),
                BKGERR=np.full(6, 0.3),
                BKGERR_APER=np.full(6, 0.3),
                BKGERR_OPT=np.full(6, 0.3),
                BKGERR_OPT_APER=np.full(6, 0.3),
                DQ=np.zeros(6, "u2"),
                BKGDQ=np.zeros(6, "u2"),
            )
            spectra.append(row)
        result = combine_spectra(spectra, [1, 3], out)
        np.testing.assert_allclose(result["FLUX"], 2.5)
        np.testing.assert_allclose(result["ERR_STAT"], np.sqrt(10) / 4)
        np.testing.assert_allclose(result["BKGERR"], 0.3)
        np.testing.assert_allclose(result["COVAR1_STAT_APER"][::2], 10 / 16)
        spectra[1]["DQ"][2] = 16
        result = combine_spectra(spectra, [1, 3], out)
        np.testing.assert_allclose(result["FLUX"][4:6], -2)
        np.testing.assert_array_equal(result["NEXP"][4:6], 1)
        np.testing.assert_array_equal(result["REJDQ"][4:6], 16)

    def test_counts_conversion_preserves_exptime_err_dq(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            hdr = fits.Header(
                {"EXPTIME": 20.0, "BUNIT": "COUNTS", "SHIFTA2": 0.25, "CRPIX2": 3.0}
            )
            data = np.arange(12, dtype="f4").reshape(3, 4) - 5
            original = fits.HDUList(
                [
                    fits.PrimaryHDU(),
                    fits.ImageHDU(data, header=hdr, name="SCI"),
                    fits.ImageHDU(np.ones_like(data), header=hdr, name="ERR"),
                    fits.ImageHDU(np.zeros_like(data, dtype="u2"), name="DQ"),
                ]
            )
            original.writeto(p / "original.fits")
            background = fits.HDUList(
                [
                    fits.PrimaryHDU(header=fits.Header({"BKG2DSUB": True})),
                    fits.ImageHDU(
                        np.full_like(data, 0.1),
                        header=fits.Header({"BUNIT": "count / s"}),
                        name="BKG",
                    ),
                ]
            )
            background.writeto(p / "background.fits")
            count_copy(p / "original.fits", p / "background.fits", p / "counts.fits")
            with fits.open(p / "counts.fits") as hd:
                np.testing.assert_allclose(hd["SCI"].data, data - 2, atol=1e-6)
                np.testing.assert_array_equal(hd["ERR"].data, original["ERR"].data)
                np.testing.assert_array_equal(hd["DQ"].data, original["DQ"].data)
                self.assertEqual(hd["SCI"].header["EXPTIME"], 20)
                self.assertEqual(hd["SCI"].header["SHIFTA2"], 0.25)
                self.assertEqual(hd[0].header["BACKCORR"], "OMIT")
                self.assertEqual(hd["SCI"].header["BUNIT"], "COUNTS")
            np.testing.assert_array_equal(
                fits.getdata(p / "original.fits", "SCI"), data
            )

    def test_invalid_aperture_rejected(self):
        with self.assertRaises(ValueError):
            ExtractionConfig(height_native_pixels=0)
        with self.assertRaises(ValueError):
            ExtractionConfig(center_arcsec=np.nan)


if __name__ == "__main__":
    unittest.main()

class HeightCalibrationTests(unittest.TestCase):
    def test_reject_nearest_height_substitution_in_point_mode(self):
        from hstpipe.stis.extraction import pctab_height
        from astropy.table import Table
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'pct.fits'
            table=Table(dict(APERTURE=['52X0.5']*4,CENWAVE=[1425]*4,
                             EXTRHEIGHT=[11,21,140,600]))
            fits.HDUList([fits.PrimaryHDU(),fits.BinTableHDU(table)]).writeto(p)
            self.assertTrue(pctab_height(p,'52X0.5',1425,140)['exact'])
            with self.assertRaisesRegex(ValueError,'would use 21'):
                pctab_height(p,'52X0.5',1425,61)
            self.assertEqual(pctab_height(p,'52X0.5',1425,61,False)['pctab_height'],21)

class IndependentCoadditionTests(unittest.TestCase):
    def test_method_masks_have_independent_exposure_normalization(self):
        wave=np.arange(4.)+1400;rows=[]
        for value in (2.,8.):
            d={key:np.full(4,value) for key in ('FLUX','FLUX_APER','FLUX_OPT','FLUX_OPT_APER')}
            d.update({key:np.ones(4) for key in ('ERR_STAT','ERR_STAT_APER','ERR_STAT_OPT','ERR_STAT_OPT_APER','BKGERR','BKGERR_APER','BKGERR_OPT','BKGERR_OPT_APER')})
            d.update(WAVELENGTH=wave,DQ=np.zeros(4,'u2'),DQ_BOX=np.zeros(4,'u2'),DQ_OPT=np.zeros(4,'u2'),BKGDQ=np.zeros(4,'u2'),PROFCOV_OPT=np.full(4,.9),INFOFRAC_OPT=np.full(4,.95))
            rows.append(d)
        rows[1]['DQ_BOX'][0]=4
        rows[1]['DQ_OPT'][1]=4
        rows[1]['FLUX'][0]=np.nan
        rows[1]['FLUX_OPT'][1]=np.nan
        out=combine_spectra(rows,[1,3],wave)
        np.testing.assert_allclose(out['FLUX'],[2,6.5,6.5,6.5])
        np.testing.assert_allclose(out['FLUX_OPT'],[6.5,2,6.5,6.5])
        np.testing.assert_array_equal(out['NEXP_BOX'],[1,2,2,2])
        np.testing.assert_array_equal(out['NEXP_OPT'],[2,1,2,2])
        np.testing.assert_array_equal(out['EXPTIME_BOX'],[1,4,4,4])
        np.testing.assert_array_equal(out['EXPTIME_OPT'],[4,1,4,4])
        np.testing.assert_allclose(out['PROFCOV_OPT'],.9)
        self.assertAlmostEqual(out['ERR_STAT_OPT'][0],np.sqrt(10)/4)
        self.assertEqual(out['ERR_STAT'][0],1)
