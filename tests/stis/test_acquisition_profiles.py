import numpy as np
import pytest
from hstpipe.stis.acquisition_profiles import slit_bin_weights,polygon_area,clip_rectangle


def test_polygon_strip_collapse_conserves_and_masks_partial_pixels():
    squares=np.array([[[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]],[[-.5,.5],[.5,.5],[.5,1.5],[-.5,1.5]]])
    w=slit_bin_weights(squares,np.array([-.5,0,.5,1.,1.5]),0,.5)
    assert np.allclose(w,[[.25,0],[.25,0],[0,.25],[0,.25]])
    assert np.allclose(w.sum(axis=0),[.5,.5])
    assert (w@np.array([8.,4.])).sum()==pytest.approx(6.)


def test_polygon_clipping_retains_rotated_pixel_area():
    p=np.array([[0.,-1.],[1.,0.],[0.,1.],[-1.,0.]])
    assert polygon_area(p)==pytest.approx(2.)
    assert polygon_area(clip_rectangle(p,-10,10,-10,10))==pytest.approx(2.)
    assert polygon_area(clip_rectangle(p,0,10,-10,10))==pytest.approx(1.)
    assert polygon_area(clip_rectangle(p,2,3,-10,10))==0
