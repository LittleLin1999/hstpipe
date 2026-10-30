"""Independent synthetic checks of continuous residual-sky interpolation."""

import unittest
import numpy as np
from hstpipe.stis.background import fit, assess_x1d_residual


class ContinuousTests(unittest.TestCase):
    def test_signed_constant_and_boundaries(self):
        y, x = np.indices((400, 350))
        sy = y - 200.0
        wave = 1100 + x * 1.8
        valid = np.ones(y.shape, bool)
        data = np.full(y.shape, -0.001)
        b, err, outside, records = fit(data, np.zeros(y.shape), wave, sy, valid)
        np.testing.assert_allclose(b, data, atol=1e-15)
        self.assertTrue(np.isfinite(err).all())
        # The former +/-98 boundary receives exactly the same model.
        np.testing.assert_allclose(b[102], b[101], atol=1e-15)
        self.assertTrue(outside.any())

    def test_source_protection_and_gradient(self):
        y, x = np.indices((400, 350))
        sy = y - 200.0
        wave = 1100 + x * 1.8
        valid = np.ones(y.shape, bool)
        sky = 0.002 + sy * 0.000002
        signal = np.where(abs(sy) < 80, 0.01, 0)
        b, _, _, _ = fit(sky + signal, np.zeros(y.shape), wave, sy, valid)
        roi = abs(sy) < 80
        np.testing.assert_allclose(b[roi], sky[roi], atol=1e-14)
        np.testing.assert_allclose((sky + signal - b)[roi], signal[roi], atol=1e-14)

    def test_sparse_mama_residual_rejects_negative_clipped_pedestal(self):
        # One positive photon among many negative dark-subtracted sky pixels
        # can leave the arithmetic mean at zero while a clipped estimator
        # approaches the negative pedestal.
        rate = np.full((450, 300), -1e-5)
        trace = np.full(300, 225.)
        rate[113] = 4.9e-4
        args = (rate, np.zeros(rate.shape, "u2"), np.ones(rate.shape, bool),
                trace, np.linspace(1400, 1600, 300))
        biased = assess_x1d_residual(
            *args, np.full(300, -11e-5), 11, 112.5, 25
        )
        self.assertTrue(biased["reject_calstis_residual"])
        unbiased = assess_x1d_residual(
            *args, np.full(300, biased["unclipped_mean_rate"] * 11), 11, 112.5, 25
        )
        self.assertFalse(unbiased["reject_calstis_residual"])


if __name__ == "__main__":
    unittest.main()


class SparseAirglowEdgeTests(unittest.TestCase):
    def test_thin_edge_strips_preserve_airglow(self):
        # Spatial edges leave only 9 rows in the outermost nominal 50-row strip.
        # Those bins fail n>=12 while their broad continuum bins pass. The line
        # must come from interior rows at the same wavelength, not continuum.
        y, x = np.indices((418, 1110))
        sy = y - 209.0
        wave = 1100 + (x + 0.5) * 0.584
        valid = np.ones(y.shape, bool)
        data = 0.001 + 0.1 * np.exp(-0.5 * ((wave - 1216) / 2) ** 2)
        b, error, outside, _ = fit(data, np.ones(y.shape) * 0.001, wave, sy, valid)
        peak = abs(wave[0] - 1216) < 1
        interior = b[209, peak].mean()
        self.assertGreater(interior, 0.08)
        np.testing.assert_allclose(
            b[:9, peak], np.broadcast_to(b[209, peak], b[:9, peak].shape), rtol=0.02
        )
        np.testing.assert_allclose(
            b[-9:, peak], np.broadcast_to(b[209, peak], b[-9:, peak].shape), rtol=0.02
        )
        self.assertTrue(outside[:9, peak].all())


class AdditionalResidualTests(unittest.TestCase):
    def setUp(self):
        from hstpipe.stis.config import BackgroundConfig
        self.options = BackgroundConfig()
        y, x = np.indices((240, 900))
        self.sy = y - 120.
        self.wave = 1150 + x * .7
        self.valid = np.ones(y.shape, bool)
        self.variance = np.full(y.shape, 1e-8)
        self.bound = np.zeros(y.shape)
        self.noise = np.random.default_rng(14).normal(0, 1e-4, y.shape)

    def correction(self, rate, valid=None):
        from hstpipe.stis.background import residual_correction
        return residual_correction(rate, self.variance, self.bound,
            self.wave, self.sy, self.valid if valid is None else valid, self.options, .7)

    def test_signed_constant_recovery_and_source_protection(self):
        for level in (-.001, .001):
            data = level + self.noise
            model, error, flags, report = self.correction(data)
            self.assertTrue(report["components"][0]["applied"])
            np.testing.assert_allclose(model, level, atol=1e-5)
            protected = data + np.where(abs(self.sy) < 80, .1, 0)
            np.testing.assert_array_equal(self.correction(protected)[0], model)
            self.assertTrue((error > 0).all())

    def test_zero_residual_does_not_force_subtraction(self):
        model, _, _, report = self.correction(self.noise)
        np.testing.assert_array_equal(model, 0)
        self.assertFalse(report["applied"])

    def test_validation_failure_and_gradient_are_not_subtracted(self):
        training_only = self.noise + np.where(
            (abs(self.sy) >= 84) & (abs(self.sy) < 90), .01, 0)
        model, _, _, report = self.correction(training_only)
        np.testing.assert_array_equal(model, 0)
        self.assertFalse(report["components"][0]["applied"])
        gradient = self.noise + .001 + self.sy * .00005
        model, _, _, report = self.correction(gradient)
        np.testing.assert_array_equal(model, 0)

    def test_no_valid_sky_is_unconstrained(self):
        mask = abs(self.sy) < 80
        model, _, flags, report = self.correction(self.noise, mask)
        self.assertFalse(report["applied"])
        self.assertEqual(report["components"][0]["reason"], "unconstrained")
        np.testing.assert_array_equal(model, 0)

    def test_sparse_photons_and_zero_pipeline_errors(self):
        rng = np.random.default_rng(3)
        self.variance[:] = 0
        data = (rng.poisson(.02, self.sy.shape) - .02) / 1000
        model, _, _, report = self.correction(data)
        np.testing.assert_array_equal(model, 0)
        self.assertGreater(report["components"][0]["measurements"]["fit"]["error"], 0)

    def test_airglow_profile_not_continuum_or_second_dark(self):
        line = .02 * np.exp(-.5 * ((self.wave - 1216) / 1.5)**2)
        model, _, _, report = self.correction(self.noise + line)
        self.assertFalse(report["components"][0]["applied"])
        self.assertTrue(report["components"][1]["applied"])
        self.assertTrue(np.all(model[self.wave > 1400] == 0))
        np.testing.assert_allclose(model[120].sum(), line[120].sum(), rtol=.04)
        self.assertLess(abs(model[120].max() / line[120].max() - 1), .1)

    def test_disabled_primary_disables_all_background(self):
        from types import SimpleNamespace
        from hstpipe.stis.background import BackgroundResult, add_residual
        zeros = np.zeros(self.sy.shape)
        primary = BackgroundResult(zeros, zeros, zeros.astype('u2'), {'enabled': False})
        exposure = SimpleNamespace(counts=(self.noise + .001) * 30,
                                   variance=self.variance * 30**2, time=30)
        result = add_residual(exposure, primary, self.wave, self.sy, self.valid, self.options, .7)
        np.testing.assert_array_equal(result.rate, 0)
        self.assertFalse(result.metadata['residual']['enabled'])


class RectifiedBackgroundTests(unittest.TestCase):
    def test_physical_units_signed_sky_and_source_preservation(self):
        from hstpipe.stis.background import estimate_rectified
        from hstpipe.stis.config import BackgroundConfig
        y,x=np.indices((400,350));sy=y-200.;wave=1400+x*.7
        valid=np.ones(y.shape,bool)
        sky=(-1+sy*.001)*1e-15
        source=np.where(abs(sy)<80,3e-15,0.)
        options=BackgroundConfig(spatial_widths=(50,),bootstrap_draws=10)
        first=estimate_rectified(sky+source,np.full(y.shape,1e-34),wave,sy,valid,options,.7)
        np.testing.assert_allclose(first.rate[abs(sy)<80],sky[abs(sy)<80],rtol=1e-12)
        np.testing.assert_allclose((sky+source-first.rate)[abs(sy)<80],source[abs(sy)<80],rtol=1e-12)
        scaled=estimate_rectified((sky+source)*7,np.full(y.shape,49e-34),wave,sy,valid,options,.7)
        np.testing.assert_allclose(scaled.rate,first.rate*7,rtol=1e-12)
        np.testing.assert_allclose(scaled.error,first.error*7,rtol=1e-9,atol=1e-30)

    def test_rectified_residual_zero_and_physical_sampling(self):
        from hstpipe.stis.background import rectified_residual
        from hstpipe.stis.config import BackgroundConfig
        y,x=np.indices((480,400));sy=(y-240)/2;wave=1400+x*.7
        noise=np.random.default_rng(21).normal(0,1e-18,y.shape)
        valid=np.ones(y.shape,bool);err=np.full(y.shape,1e-18)
        r=rectified_residual(noise,err,err*0,wave,sy,valid,BackgroundConfig(),.7,.5)
        np.testing.assert_array_equal(r.rate,0)
        self.assertEqual(r.metadata['cluster_rows'],8)
        r=rectified_residual(noise-1e-15,err,err*0,wave,sy,valid,BackgroundConfig(),.7,.5)
        self.assertTrue(r.metadata['applied'])
        np.testing.assert_allclose(r.rate,-1e-15,atol=1e-18)
