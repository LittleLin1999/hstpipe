import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from hstpipe.stis.acquisition import slit_sky_vertices,target_anchored_wcs


def test_slit_width_length_and_orientation():
    center=SkyCoord(155*u.deg,8*u.deg)
    vertices=slit_sky_vertices(center,0,52,.5)
    assert vertices[0].separation(vertices[1]).arcsec==pytest.approx(.5,rel=1e-7)
    assert vertices[1].separation(vertices[2]).arcsec==pytest.approx(52,rel=1e-7)
    east,north=center.spherical_offsets_to(vertices)
    assert np.allclose(north.arcsec,[-26,-26,26,26,-26])


def test_target_anchor_does_not_mutate_header():
    w=WCS(naxis=2);w.wcs.crpix=[51,51];w.wcs.crval=[155,8];w.wcs.cdelt=[-.05/3600,.05/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    h=w.to_header();copy=h.copy();target=SkyCoord(155.1*u.deg,8.1*u.deg)
    anchored=target_anchored_wcs(h,[47.1,50.9],target)
    assert np.allclose(anchored.world_to_pixel(target),[47.1,50.9])
    assert h==copy


def test_cutout_coordinates_put_target_at_zero_with_pixel_edges():
    from hstpipe.stis.acquisition import target_offset_geometry
    w=WCS(naxis=2);w.wcs.crpix=[51,51];w.wcs.crval=[155,8];w.wcs.cdelt=[-.05/3600,.05/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    scale,extent=target_offset_geometry(w.to_header(),[47.1,50.9],(100,100))
    assert np.allclose(scale,[.05,.05])
    assert np.allclose(extent,[-2.38,2.62,-2.57,2.43])
    assert extent[0]+(47.1+.5)*scale[0]==pytest.approx(0.)
    assert extent[2]+(50.9+.5)*scale[1]==pytest.approx(0.)
