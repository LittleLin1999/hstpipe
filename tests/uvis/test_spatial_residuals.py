import numpy as np
from hstpipe.uvis.astrometry.relative_fit import quality
from hstpipe.uvis.config import settings

def case(amplitude):
    rng=np.random.default_rng(72);xy=rng.uniform(-1000,1000,(30,2))
    residual=np.column_stack([amplitude*xy[:,0]/1000,amplitude*xy[:,1]/1000])
    return xy,xy+residual+rng.normal(0,.001,xy.shape)

def test_high_correlation_small_amplitude_preserves_valid_registration():
    o=settings({})['astrometry']['relative'];o=dict(o,pixel_scale_arcsec=.04)
    x,y=case(.03);r=quality(x,y,np.full(len(x),.01),o)
    assert r['status']=='validated'
    t=r['spatial_systematic'];assert t['max_abs_correlation']>.9
    assert t['amplitude_pix']<.06 and not t['exceeds_accuracy']
    assert np.isclose(t['amplitude_mas'],40*t['amplitude_pix'])

def test_significant_local_systematic_not_hidden_by_small_average_rms():
    o=settings({})['astrometry']['relative'];o=dict(o,max_scatter_pix=.5,scatter_floor_pix=.5,max_transform_error_pix=1.,max_current_rms_pix=1.,pixel_scale_arcsec=.04)
    rng=np.random.default_rng(4);x=rng.uniform(-1000,1000,(60,2))
    y=x+np.column_stack([x[:,0]*.0007,np.zeros(len(x))])
    r=quality(x,y,np.full(len(x),.05),o)
    assert r['postfit_scatter_pix']<.5
    assert r['spatial_systematic']['significant'] and r['spatial_systematic']['exceeds_accuracy']
    assert r['status'] not in ('validated','refinement_supported')

def test_sparse_or_collinear_trend_is_underconstrained_not_automatic_failure():
    o=settings({})['astrometry']['relative']
    for x in (np.array([[0.,0.],[100.,100.],[200.,200.]]),np.array([[0.,0.],[100.,0.]])):
        r=quality(x,x+.01,np.full(len(x),.08),o,identities_verified=True)
        assert r['spatial_systematic']['status']=='underconstrained'
        assert r['status']=='validated'

def test_large_uncorrelated_residuals_still_fail():
    rng=np.random.default_rng(3);x=rng.uniform(0,1000,(30,2));o=settings({})['astrometry']['relative']
    r=quality(x,x+rng.normal(0,3,x.shape),np.full(len(x),.08),o)
    assert r['status'] not in ('validated','refinement_supported')
