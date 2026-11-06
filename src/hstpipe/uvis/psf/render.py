"""Native empirical PSF rendering, preserving the finite stamp's center."""
import numpy as np
from scipy.ndimage import map_coordinates


def render_empirical(model, dx, dy):
    model = np.asarray(model, dtype=float)
    if (model.ndim != 2 or model.shape[0] != model.shape[1]
            or model.shape[0] % 2 != 1 or not np.isfinite(model).all()
            or model.sum() <= 0):
        raise ValueError('Expected a finite, odd square empirical PSF with positive flux')
    center = (model.shape[0] - 1) / 2
    return map_coordinates(model, [dy + center, dx + center], order=3,
                           mode='constant', cval=0., prefilter=True)
