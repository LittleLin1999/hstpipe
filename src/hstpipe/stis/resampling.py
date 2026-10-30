"""Drizzle surface brightness with iscale=1 and squared-weight variance."""

import numpy as np


def drizzle(data, variance, weight, mapping, shape, pixfrac=1):
    from drizzle.resample import Drizzle

    d = Drizzle(out_shape=shape, fillval="NaN", fillval2="NaN")
    d.add_image(
        np.asarray(data, "f4"),
        exptime=1,
        pixmap=mapping,
        data2=np.asarray(variance, "f4"),
        weight_map=np.asarray(weight, "f4"),
        in_units="cps",
        iscale=1,
        pixel_scale_ratio=1,
        pixfrac=pixfrac,
    )
    return d
