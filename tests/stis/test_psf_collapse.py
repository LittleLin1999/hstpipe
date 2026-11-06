import numpy as np
import pytest
from hstpipe.stis.psf.collapse import pixel_overlap_matrix


def test_projection_conserves_flux_and_integrates_pixel_overlap():
    source=np.arange(-5.,6.);target=np.arange(-4.75,5.,.5)
    operator=pixel_overlap_matrix(source,target)
    p=np.zeros(11);p[5]=1
    assert (operator@p).sum()==pytest.approx(1.)
    assert np.allclose((operator@p)[8:12],[0,.5,.5,0])
    covariance=operator@np.eye(11)@operator.T
    assert np.linalg.eigvalsh(covariance).min()>-1e-10
    with pytest.raises(ValueError,match='support'):pixel_overlap_matrix(source,np.arange(-10.,11.))


def test_collapse_spectral_interpolation_and_aperture_guard(tmp_path):
    from astropy.io import fits
    from astropy.table import Table
    from hstpipe.stis.psf.collapse import collapse_library_profile
    hdus=[fits.PrimaryHDU()]
    profiles=[np.array([.05,.2,.5,.2,.05]),np.array([.1,.2,.4,.2,.1])]
    for i,(center,p) in enumerate(zip([1450,1550],profiles),1):
        h=fits.BinTableHDU(Table(dict(offset_arcsec=np.arange(-2.,3.),profile=p)),name=f'PSF{i:03d}')
        h.header.update(dict(APERTURE='52X0.5',EPOCH='2023',POSGROUP=4,SAMPLING='drizzled',NORMRAD=60,WAVEMIN=center-10,WAVEMAX=center+10,YMIN=423.,YMAX=423.,GRADE='WING',SCATSTAT='UNMEASURED'))
        hdus.extend([h,fits.ImageHDU(np.eye(5)*1e-5,name=f'PCOV{i:03d}')])
    lib=tmp_path/'library.fits';fits.HDUList(hdus).writeto(lib,checksum=True)
    req=dict(templates=['PSF001','PSF002'],science_aperture='52X0.5D1',sampling='drizzled',wavelength_angstrom=[1480.,1490.],position_arcsec=[-1.,0.,1.],trace_center_arcsec=[0.,0.],center_at_reference_arcsec=0.)
    with pytest.raises(ValueError,match='proxy flag'):collapse_library_profile(lib,req,tmp_path/'out.ecsv')
    req['allow_aperture_proxy']=True
    out=collapse_library_profile(lib,req,tmp_path/'out.ecsv')
    expected=(.65*profiles[0]+.35*profiles[1])[1:4];expected/=expected.max()
    assert np.allclose(out['profile_peak_normalized'],expected)
    assert np.all(np.isfinite(out['template_error_bound']))
    req['wavelength_angstrom']=[1400.,1490.]
    with pytest.raises(ValueError,match='extrapolation'):collapse_library_profile(lib,req,tmp_path/'out.ecsv')
