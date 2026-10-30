import numpy as np
from astropy.io import fits
from ..common.fits import validate_product

def inspect_flc(path):
    """DQ/artifact assessment per chip, without altering detector calibration."""
    result=[]
    with fits.open(path) as hd:
        for h in hd:
            if h.name!='SCI':continue
            v=h.header.get('EXTVER',1);dq=hd['DQ',v].data
            good=(dq==0)&np.isfinite(h.data)
            result.append(dict(extver=v,shape=h.data.shape,flagged_fraction=float(np.mean(dq!=0)),
                finite_fraction=float(np.mean(np.isfinite(h.data))),
                mean=float(np.mean(h.data[good])) if good.any() else None))
    return result

def validate_drizzle(path):
    validate_product(path,('SCI','WHT','CTX'))
    with fits.open(path) as hd:
        if hd['SCI'].data.shape!=hd['WHT'].data.shape:raise ValueError('SCI/WHT geometry mismatch')
        if np.any(hd['WHT'].data<0):raise ValueError('Negative weight')
        if not np.any(hd['WHT'].data>0):raise ValueError('Empty drizzle coverage')
