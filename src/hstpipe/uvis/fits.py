"""Materialize HST constant ERR/DQ arrays without modifying downloaded HDUs."""
import numpy as np

def image_data(hdu):
    if hdu.data is not None:return hdu.data
    h=hdu.header
    if 'PIXVALUE' not in h:raise ValueError('Empty image extension without PIXVALUE')
    return np.full((h['NPIX2'],h['NPIX1']),h['PIXVALUE'],dtype='u2' if hdu.name=='DQ' else 'f4')

