import numpy as np
import pytest
from hstpipe.stis.psf.extraction import overlap_profile,subtract_background,normalize,measures,fractional_aperture
from hstpipe.stis.psf.validation import extension_direction,extension_fit


def test_overlap_flux_covariance_and_dq():
    image=np.zeros((31,2));image[15]=[10,20];v=np.ones_like(image);good=np.ones_like(image,dtype=bool)
    y,p,cov,coverage=overlap_profile(image,v,good,np.array([15.25,15.25]),[0,1],10)
    assert np.isclose(p.sum(),30)
    assert np.isclose((y*p).sum()/p.sum(),-.25)
    # Full variance of sum equals original sum of independent contributing pixel variances.
    assert np.isclose(cov.sum(),2*(20+.75**2+.25**2))
    assert np.isclose(cov[9,10],2*.25*.75)
    good[15,1]=False
    _,q,_,_=overlap_profile(image,v,good,np.array([15.25,15.25]),[0,1],10)
    assert np.isclose(q.sum(),10)


def test_normalization_covariance():
    p=np.array([2.,6.,2.]);cov=np.diag([2.,6.,2.]);n,c=normalize(p,cov)
    assert n.sum()==pytest.approx(1)
    assert np.max(abs(c.sum(axis=0)))<1e-15
    assert np.allclose(c,c.T)
    with pytest.raises(ValueError):normalize(-p,cov)


def test_spatial_background_models():
    y=np.arange(-190,191,dtype=float);source=300*np.exp(-y*y/8);bg=10+.02*y
    p=source+bg;cov=np.eye(len(y))*4
    q,cv,model,mask,_,_=subtract_background(y,p,cov,np.ones(len(y)),[100,160],1)
    assert np.max(abs(q-source))<1e-10
    assert np.max(abs(model-bg))<1e-10
    assert np.all(abs(y[mask])>=100)
    assert cv[190,190]>cov[190,190]
    assert np.min(np.linalg.eigvalsh(cv))>-1e-10


def test_centroid_width_and_apertures():
    y=np.arange(-70,71,dtype=float);p=np.exp(-.5*((y-.3)/1.8)**2)*1000;cov=np.diag(np.maximum(p,1))
    m,n,cv=measures(y,p,cov,[3,5,7,9,11,15,21,41,81,140])
    assert abs(m['centroid']-.3)<.01
    assert abs(m['fwhm']-2.35482*1.8)<.15
    assert m['centroid_error']>0
    assert 0<m['ee_3']<m['ee_11']<=1
    assert fractional_aperture(y,3).sum()==3


def test_extended_component_recovery():
    y=np.arange(-70,71,dtype=float);p=np.exp(-.5*(y/1.8)**2);p/=p.sum()
    kernel=np.exp(-abs(y)/(.2/.0246));kernel/=kernel.sum();broad=np.convolve(p,kernel,mode='same');broad/=broad.sum();extended=.9*p+.1*broad;cov=np.eye(len(y))*1e-7
    f,e,chi=extension_fit(y,extended,p,cov,.2)
    assert f==pytest.approx(.1,abs=1e-5)
    assert e>0 and chi<1e-8
    f,e,chi=extension_fit(y,p,p,cov,.2)
    assert abs(f)<1e-10


def test_psf_config_uses_reduction_schema(tmp_path):
    import yaml
    from hstpipe.stis.psf.common import config
    base=dict(input=str(tmp_path/"data"),output=str(tmp_path/"redux"),instrument="stis")
    for name,options in [("selection.yaml",dict(psf_selection=dict(priority_programs=[18169]))),("reduction.yaml",dict(reference_dir=str(tmp_path/"refs"),spatial_scale=.0123,pixfrac=1)),("psf.yaml",dict(psf=dict(bands=[[1330,1380]])))]:
        (tmp_path/name).write_text(yaml.safe_dump(dict(base,options=options)))
    c=config(tmp_path/"psf.yaml")
    assert c["priority_programs"]==[18169]
    assert c["workdir"]==str(tmp_path)
    assert c["hstpipe_spatial_scale"]==.0123


def test_nested_stis_stage_paths(tmp_path):
    from hstpipe.common.paths import ReductionPaths
    paths=ReductionPaths(tmp_path,'GD71','stis','G140L')
    assert paths.area('intermediate','flt/inputs')==paths.directory/'intermediate/flt/inputs'
    with pytest.raises(ValueError):paths.area('intermediate','../escape')
    with pytest.raises(ValueError):paths.area('intermediate','/tmp/escape')


def test_half_pixel_aperture_boundaries():
    y=np.arange(-1.75,2,.5)
    a=fractional_aperture(y,1.)
    assert np.all(a[abs(y)==.25]==1)
    assert np.all(a[abs(y)>.5]==0)


def test_repeat_scatter_is_template_uncertainty(tmp_path):
    from hstpipe.stis.psf.psf import combine
    items=[]
    for i,p in enumerate([np.array([.1,.8,.1]),np.array([.2,.6,.2])]):
        path=tmp_path/f'p{i}.npz'
        np.savez(path,normalized=p,normalized_covariance=np.eye(3)*1e-6)
        items.append(dict(target='one-star',snr=100,profile_file=str(path)))
    mean,cm,pred,scatter,ns=combine(items)
    assert ns==1
    assert np.allclose(mean,[.15,.7,.15])
    assert pred[1,1]>.01
    assert scatter[1,1]>.01
    assert cm[1,1]<pred[1,1]


def test_psf_interpolation_requires_matching_sampling_and_bracketing(tmp_path,monkeypatch):
    import json
    from hstpipe.stis.psf.psf import interpolate_library
    monkeypatch.chdir(tmp_path);directory=tmp_path/'products/combined_psfs';directory.mkdir(parents=True)
    idx=[]
    for number,band in enumerate(['1330-1380','1420-1480']):
        name=f'PSF{number:03d}';np.savez(directory/(name+'.npz'),y=np.arange(3),profile=np.array([.2,.6,.2]),covariance_prediction=np.eye(3)*.01)
        idx.append(dict(extension=name,aperture='52X2',band=band,epoch='2026',position_group=4,sampling='hstpipe_native',normalization_radius=70))
    (directory/'index.json').write_text(json.dumps(idx))
    with pytest.raises(ValueError):interpolate_library('52X2',1200,2026,4)
    with pytest.raises(ValueError):interpolate_library('52X2',1400,2026,4,sampling='drizzled')
    y,p,e=interpolate_library('52X2',1400,2026,4)
    assert p.sum()==pytest.approx(1)
    assert np.all(e>0)


def test_fractional_wavelength_tilt_and_overlapping_masks():
    from hstpipe.stis.psf.product_profiles import wavelength_weights
    w=np.array([[0.,1.,2.,3.],[.25,1.25,2.25,3.25]])
    weights=wavelength_weights(w,.5,2.5,[])
    assert np.allclose(weights[0],[0,1,1,0])
    assert np.allclose(weights[1],[.25,1,.75,0])
    masked=wavelength_weights(w,.5,2.5,[[1.,2.],[1.5,2.5]])
    assert np.allclose(masked[0],[0,.5,0,0])


def test_profile_metrics_do_not_depend_on_flux_units():
    y=np.arange(-70,71,dtype=float);p=np.exp(-.5*(y/2)**2)*1e4;cv=np.diag(np.maximum(p,1.))
    a,_,_=measures(y,p,cv,[3,11,21])
    b,_,_=measures(y,p*1e-15,cv*1e-30,[3,11,21])
    for key in ['snr','wing_snr','fwhm','ee_11','wing_fraction']:
        assert a[key]==pytest.approx(b[key],rel=1e-5)


def test_requested_d1_aperture_is_not_regular_hardware_slit():
    from hstpipe.stis.psf.common import aperture_configuration
    assert aperture_configuration(dict(APERTURE='52X0.2',PROPAPER='52X0.2D1'))=='52X0.2D1'
    assert aperture_configuration(dict(APERTURE='52X0.5',PROPAPER='52X0.5'))=='52X0.5'


def test_wing_boundary_is_configurable():
    y=np.arange(-70,71,dtype=float);p=np.exp(-abs(y)/8)*1000;cv=np.diag(np.maximum(p,1))
    inner,_,_=measures(y,p,cv,[11],wing_inner=5)
    outer,_,_=measures(y,p,cv,[11],wing_inner=15)
    assert inner['wing_fraction']>outer['wing_fraction']>0


def test_normalized_point_source_background_is_not_extension():
    y=np.arange(-70,71,dtype=float);p=np.exp(-.5*(y/2)**2);p/=p.sum()
    background=1e-4;data=(p+background)/(1+len(p)*background)
    f,e,chi=extension_fit(y,data,p,np.eye(len(y))*1e-7,.2)
    assert abs(f)<1e-5
    assert e>0 and chi<1e-7


def test_library_plot_common_aperture_and_density():
    from hstpipe.stis.psf.library_plots import common_window_density
    y=np.arange(-2.,3.);p=np.array([.1,.2,.4,.2,.1]);cov=np.eye(5)*.001
    x,density,dcov=common_window_density(y,p,cov,2.)
    a=fractional_aperture(y,4.)
    assert a@(density*.0246)==pytest.approx(1.)
    assert np.allclose(x,y*.0246)
    assert np.allclose(dcov,dcov.T)
    assert np.linalg.eigvalsh(dcov).min()>-1e-9
    assert abs(a@dcov@a)<1e-9
    with pytest.raises(ValueError,match='not covered'):
        common_window_density(y,p,cov,3.)


def test_plot_window_respects_drizzle_pixel_edges():
    from hstpipe.stis.psf.library_plots import common_window_density
    y=np.arange(-1.75,2.,.5)+.002
    p=np.ones(len(y))/len(y);cv=np.eye(len(y))*.001
    with pytest.raises(ValueError,match='not covered'):
        common_window_density(y,p,cv,2.)
    x,density,cov=common_window_density(y,p,cv,1.5)
    assert fractional_aperture(y,3.)@(density*.5*.0246)==pytest.approx(1.)


def test_plot_selection_preserves_slit_grade_and_supported_radius():
    from astropy.io import fits
    from hstpipe.stis.psf.library_plots import select_templates
    hdus=[]
    for number,ap,grade,radius,stars in [(1,'52X0.5','WING',70,1),(2,'52X0.5','WING',60,2),(3,'52X0.5D1','WING',70,9),(4,'52X0.5','CORE',70,9),(5,'52X0.5','WING',40,9)]:
        h=fits.BinTableHDU(name=f'PSF{number:03d}')
        h.header.update(dict(APERTURE=ap,GRADE=grade,NORMRAD=radius,NSTARS=stars,NEXP=stars,EPOCH='2023',SAMPLING='hstpipe_native',WAVEMIN=1520,WAVEMAX=1580))
        hdus.append(h)
    selected=select_templates(hdus,(1520,1580),'hstpipe_native',['52X0.5'],50)
    assert [h.name for h in selected]==['PSF002']


def test_position_aperture_ratio_covariance_and_scale():
    from hstpipe.stis.psf.position_plots import aperture_ratio
    p=np.array([1.,2.,3.]);cov=np.eye(3)*.1
    a=np.array([1.,0.,0.]);b=np.ones(3)
    value,error=aperture_ratio(p,cov,a,b)
    value2,error2=aperture_ratio(p*5,cov*25,a,b)
    assert value==pytest.approx(1/6)
    assert value2==pytest.approx(value) and error2==pytest.approx(error)
    assert aperture_ratio(p,cov,b,b)==pytest.approx((1.,0.))
    with pytest.raises(ValueError):aperture_ratio(-p,cov,a,b)
