"""Target-independent configuration, mapping and scientifically meaningful drizzle checks."""

import unittest
import tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from hstpipe.stis.config import Config, BackgroundConfig
from hstpipe.stis.geometry import Grid, rectified_coordinates, make_grid
from hstpipe.stis.background import wavelength_edges
from hstpipe.stis.coadd import Coadd


class WorkflowTests(unittest.TestCase):
    def test_no_redshift_or_line_parameters(self):
        for kwargs in ({"redshift": 0.2}, {"line_windows": {}}):
            with self.assertRaises(TypeError):
                Config(**kwargs)
        self.assertEqual(Config(name="AnotherTarget").spatial_scale, 0.0123)

    def test_product_prefix_is_fits_target_and_opt_elem(self):
        exposure = SimpleNamespace(p={"TARGNAME": "J1022+0841", "OPT_ELEM": "G140L"})
        config = Config()
        self.assertEqual(config.resolve_name([exposure]), "J1022+0841_G140L")
        self.assertEqual(config.product_name("s2d_drz"), "J1022+0841_G140L_s2d_drz.fits")
        self.assertEqual(config.product_name("x2d_coadd"), "J1022+0841_G140L_x2d_coadd.fits")
        self.assertEqual(config.product_name("x1d"), "J1022+0841_G140L_x1d.fits")
        with self.assertRaises(ValueError):
            Config(name="J1022_G140L").resolve_name([exposure])

    def test_coverage_driven_edges(self):
        edges = wavelength_edges(np.array([1450.0, 1800.0]))
        self.assertGreaterEqual(edges.min(), 1400)
        self.assertGreaterEqual(edges.max(), 1800)
        self.assertFalse(np.any((edges > 1200) & (edges < 1230)))

    def test_dither_mapping(self):
        h = {
            "CD2_2": 0.0246 / 3600,
            "CD1_1": 0.584,
            "CRVAL1": 1400.0,
            "CRPIX1": 1.0,
            "CRPIX2": 5.0,
        }
        wave, spatial = rectified_coordinates(h, (10, 5), 0.235)
        self.assertAlmostEqual(spatial[4, 0], -0.235)
        grid = Grid(1400, 0.584, 5, -1.0, 0.0123, 200)
        mapping = grid.map(wave, spatial)
        self.assertAlmostEqual(mapping[4, 0, 1], (1 - 0.235) / 0.0123)

    def test_count_rate_units_and_unequal_times(self):
        grid = Grid(1000, 1, 5, 0, 1, 5)
        c = Coadd(grid)
        yy, xx = np.indices((5, 5), dtype=float)
        mapping = np.stack((xx, yy), axis=-1)
        mask = np.ones((5, 5), bool)
        for t in [10.0, 30.0]:
            # Data already surface brightness; exposure only changes its weight.
            c.add(
                np.full((5, 5), -2e-15),
                np.full((5, 5), 1e-30),
                mapping,
                mask,
                mask,
                t,
                dict(input=str(t), exptime=t),
            )
        a = c.arrays()
        np.testing.assert_allclose(a["SCI"], -2e-15, rtol=2e-7)
        np.testing.assert_allclose(
            a["ERR"], np.sqrt(10**2 + 30**2) / 40 * 1e-15, rtol=2e-7
        )
        np.testing.assert_allclose(a["EXPMAP"], 40)
        np.testing.assert_array_equal(a["NEXP"], 2)

    def test_half_pixel_sampling_uniform_brightness(self):
        grid = Grid(1000, 1, 5, 0, 0.5, 9)
        c = Coadd(grid)
        yy, xx = np.indices((5, 5), dtype=float)
        mapping = np.stack((xx, yy * 2), axis=-1)
        mask = np.ones((5, 5), bool)
        c.add(
            np.full((5, 5), 3e-15),
            np.full((5, 5), 1e-30),
            mapping,
            mask,
            mask,
            10,
            dict(input="test", exptime=10),
        )
        a = c.arrays()
        np.testing.assert_allclose(a["SCI"][a["NEXP"] > 0], 3e-15, rtol=2e-7)

    def test_invalid_background_geometry_rejected(self):
        with self.assertRaises(ValueError):
            BackgroundConfig(sky_inner=20)

    def test_default_grid_encloses_actual_pointings(self):
        h = dict(
            CD1_1=0.584,
            CRVAL1=1400.0,
            CRPIX1=1.0,
            NAXIS1=20,
            CD2_2=0.0246 / 3600,
            CRPIX2=6.0,
            NAXIS2=10,
        )
        grid = make_grid([h, h], Config(), [3.0, 3.235])
        for offset in (3.0, 3.235):
            wave, spatial = rectified_coordinates(h, (10, 20), offset)
            mapped = grid.map(wave, spatial)
            self.assertGreaterEqual(mapped[..., 1].min(), -1e-9)
            self.assertLessEqual(mapped[..., 1].max(), grid.ny - 1 + 1e-9)

    def test_fractional_variance_and_masked_coverage(self):
        grid = Grid(1000, 1, 5, 0, 1, 7)
        c = Coadd(grid)
        yy, xx = np.indices((5, 5), dtype=float)
        mapping = np.stack((xx, yy + 0.25), axis=-1)
        mask = np.ones((5, 5), bool)
        mask[:, 2] = False
        dq = np.zeros((5, 5), "u2")
        dq[:, 2] = 16
        c.add(
            np.full((5, 5), -1e-15),
            np.full((5, 5), 1e-30),
            mapping,
            mask,
            np.ones_like(mask),
            10,
            dict(input="test", exptime=10),
            input_dq=dq,
        )
        a = c.arrays()
        np.testing.assert_allclose(
            a["ERR"][1:5, 1], np.sqrt(0.75**2 + 0.25**2) * 1e-15, rtol=1e-6
        )
        self.assertTrue(np.isnan(a["SCI"][:, 2]).all())
        self.assertTrue((a["DQ"][:, 2] == 4).all())
        self.assertTrue((a["NEXP"][:, 2] == 0).all())
        self.assertTrue((a["REJDQ"][:5, 2] == 16).all())

    def test_subtraction_commutes_with_fixed_drizzle(self):
        grid = Grid(1000, 1, 5, 0, 0.5, 12)
        base, clean = Coadd(grid), Coadd(grid)
        yy, xx = np.indices((5, 5), dtype=float)
        mask = np.ones((5, 5), bool)
        for shift, time in [(0.0, 10), (0.35, 25)]:
            mapping = np.stack((xx, (yy + shift) * 2), axis=-1)
            science = (2 + np.sin(yy)) * 1e-15
            model = (0.1 + yy * 0.2) * 1e-15
            variance = np.ones_like(science) * 1e-30
            args = (
                variance,
                mapping,
                mask,
                mask,
                time,
                dict(input=str(time), exptime=time),
            )
            base.add(science, *args)
            clean.add(
                science - model, *args, model=model, model_error=np.zeros_like(model)
            )
        a, b = base.arrays(), clean.arrays(background=True)
        use = b["NEXP"] > 0
        np.testing.assert_allclose(b["SCI"][use], (a["SCI"] - b["BKG"])[use], rtol=1e-6)
