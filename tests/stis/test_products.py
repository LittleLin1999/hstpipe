"""Product calibration modes and final FITS structure."""

import tempfile
import unittest
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from hstpipe.stis.config import Config, ExtractionConfig
from hstpipe.stis.extraction import overlap_matrix
from hstpipe.stis.products import (
    branch_columns, morphology_conversion, point_conversion, write_consolidated
)


class ProductTests(unittest.TestCase):
    def test_calibration_modes_and_model_bounds(self):
        wave = np.array([1000.0, 1001.0, 1002.0])
        exposure = {"WAVELENGTH": wave, "CAL_RATIO": np.full(3, 0.5),
                    "DQ": np.zeros(3, "u2")}
        np.testing.assert_allclose(
            point_conversion([exposure], [5], wave, overlap_matrix), 2
        )
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "model.fits"
            fits.HDUList([
                fits.PrimaryHDU(),
                fits.BinTableHDU(Table({"WAVELENGTH": wave, "APFRAC": [0.5, 0.5, 0.5]})),
            ]).writeto(path)
            np.testing.assert_allclose(morphology_conversion(path, wave), 2)
            with self.assertRaises(ValueError):
                ExtractionConfig(flux_mode="morph")

    def test_three_extensions_and_separate_errors(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            source = d / "source.fits"
            fits.PrimaryHDU().writeto(source)
            wave = np.array([1400., 1401.])
            data = dict(WAVELENGTH=wave, FLUX=np.array([1., 2.]),
                        ERROR=np.ones(2), FLUX_OPT=np.array([1., 2.]),
                        ERR_OPT=np.ones(2), BKGERR=np.zeros(2), BKGERR_OPT=np.zeros(2),
                        DQ=np.zeros(2, "u2"), EXPTIME=np.full(2, 10.),
                        NEXP=np.ones(2, "u2"))
            data.update(FLUX_APER=data["FLUX"], ERROR_APER=data["ERROR"],
                        FLUX_OPT_APER=data["FLUX_OPT"], ERR_OPT_APER=data["ERR_OPT"],
                        BKGERR_APER=data["BKGERR"], BKGERR_OPT_APER=data["BKGERR_OPT"])
            config = Config(name="J1022+0841_G140L")
            cols = branch_columns(data, "point", np.ones(2), native=True)
            path = d / config.product_name("x1d")
            write_consolidated(path, config,
                               {"FLT": source, "DRZ": source, "X2D": source},
                               cols, cols, cols, [source],
                               [{"input": source.name, "residual_applied": True,
                                 "residual_audit": {"difference_sigma": 0.1}}])
            with fits.open(path, checksum=True) as hd:
                self.assertEqual([h.name for h in hd], ["PRIMARY", "FLT", "DRZ", "X2D"])
                self.assertEqual(hd["FLT"].header["BACKCORR"], "OMIT")
                self.assertEqual(
                    u.Unit(hd["FLT"].columns["FLUX_BOX"].unit),
                    u.erg / (u.s * u.cm**2 * u.AA),
                )


if __name__ == "__main__":
    unittest.main()

class SeparateApertureResponseTests(unittest.TestCase):
    def test_boxcar_and_optimal_use_independent_corrections(self):
        wave=np.arange(3.)+1400
        exposure=dict(WAVELENGTH=wave,CAL_RATIO=np.full(3,.5),
                      CAL_RATIO_OPT=np.full(3,.8),DQ=np.zeros(3,'u2'))
        box=point_conversion([exposure],[10],wave,overlap_matrix)
        opt=point_conversion([exposure],[10],wave,overlap_matrix,ratio_key='CAL_RATIO_OPT')
        data={key:np.ones(3) for key in ('FLUX','ERROR','FLUX_OPT','ERR_OPT','BKGERR','BKGERR_OPT','EXPTIME','NEXP','DQ')}
        data.update(WAVELENGTH=wave,DQ=np.zeros(3,'u2'))
        result=branch_columns(data,'point',{'BOX':box,'OPT':opt})
        np.testing.assert_allclose(result['FLUX_BOX'],2)
        np.testing.assert_allclose(result['FLUX_OPT'],1.25)
        np.testing.assert_allclose(result['ERR_OPT'],1.25)

class IndependentMaskTests(unittest.TestCase):
    def test_calibration_does_not_use_science_rejection(self):
        wave=np.arange(5.)+1400
        row=dict(WAVELENGTH=wave,CAL_RATIO=np.full(5,.5),DQ=np.full(5,4,'u2'),CALDQ=np.zeros(5,'u2'))
        np.testing.assert_allclose(point_conversion([row],[10],wave,overlap_matrix),2)
        row['CALDQ'][2]=4
        result=point_conversion([row],[10],wave,overlap_matrix)
        self.assertTrue(np.isnan(result[2]))
        np.testing.assert_allclose(result[[0,1,3,4]],2)

    def test_boxcar_and_optimal_masks_do_not_invalidate_each_other(self):
        d={key:np.ones(3) for key in ('FLUX','ERROR','FLUX_OPT','ERR_OPT','BKGERR','BKGERR_OPT','EXPTIME','NEXP')}
        d.update(WAVELENGTH=np.arange(3.)+1400,DQ=np.array([4,0,0],'u2'),DQ_BOX=np.array([4,0,0],'u2'),DQ_OPT=np.array([0,4,0],'u2'))
        out=branch_columns(d,'diffuse',1)
        self.assertTrue(np.isnan(out['FLUX_BOX'][0]));self.assertEqual(out['FLUX_OPT'][0],1)
        self.assertTrue(np.isnan(out['FLUX_OPT'][1]));self.assertEqual(out['FLUX_BOX'][1],1)
        np.testing.assert_array_equal(out['DQ'],out['DQ_BOX'])
