"""Candidate spatial models preserve source exclusion and signed sky means."""

import unittest
import numpy as np
from hstpipe.stis.background import fit
from hstpipe.stis.background import scoring_cells as cells


class ExposureModelTests(unittest.TestCase):
    def test_candidate_widths_do_not_use_source(self):
        y, x = np.indices((400, 350))
        sy = y - 200.0
        wave = 1100 + x * 1.8
        sky = -0.001 + sy * 1e-6
        variance = np.full(y.shape, 1e-6)
        valid = np.ones(y.shape, bool)
        source = np.where(abs(sy) < 80, 10.0, 0)
        training = valid & ((y // 8) % 2 == 0)
        for width in [25.0, 50.0, 100.0]:
            first = fit(sky, variance, wave, sy, training, width=width)[0]
            second = fit(sky + source, variance, wave, sy, training, width=width)[0]
            np.testing.assert_array_equal(first, second)

    def test_zero_errors_do_not_remove_scoring_pixels(self):
        y, x = np.indices((300, 350))
        sy = y - 150.0
        wave = 1100 + x * 1.8
        data = -0.01 + np.sin(x) * 0.0001
        mask = abs(sy) > 100
        result = cells(data, np.zeros(y.shape), wave, sy, mask)
        self.assertGreater(len(result), 0)
        self.assertTrue(all(scale > 0 for _, scale in result))

    def test_heldout_predictions_are_evaluated(self):
        y, x = np.indices((400, 350))
        sy = y - 200.0
        wave = 1100 + x * 1.8
        valid = np.ones(y.shape, bool)
        train = valid & ((y // 8) % 2 == 0)
        data = np.full(y.shape, 0.01)
        model = fit(
            data, np.ones(y.shape) * 0.001, wave, sy, train, prediction_valid=valid
        )[0]
        np.testing.assert_allclose(model[~train], 0.01, atol=1e-14)
