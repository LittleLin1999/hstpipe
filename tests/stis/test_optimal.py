"""Profile-weighted extraction checks for signed, low-count STIS data."""

import unittest
import numpy as np
from hstpipe.stis.optimal import empirical_profile, profile_extract


class OptimalTests(unittest.TestCase):
    def test_known_profile_and_negative_background(self):
        profile = np.array([0.1, 0.2, 0.4, 0.2, 0.1])
        aperture = np.ones((5, 3))
        scene = 10 * profile[:, None] + np.array([-1.0, 0.0, 1.0])[None, :]
        result = profile_extract(
            scene, np.ones_like(scene), np.zeros_like(scene),
            np.ones_like(scene, bool), aperture, profile,
        )
        self.assertAlmostEqual(result["FLUX_OPT"][1], 10.0)
        self.assertLess(result["FLUX_OPT"][0], result["FLUX_OPT"][1])
        self.assertGreater(result["FLUX_OPT"][2], result["FLUX_OPT"][1])
        self.assertTrue(np.all(np.isfinite(result["ERR_OPT"])))

    def test_bad_pixel_and_zero_error(self):
        image = np.ones((5, 20))
        valid = np.ones_like(image, bool)
        valid[2, 8] = False
        aperture = np.ones_like(image)
        profile = empirical_profile(image, valid, np.arange(20) + 1400, aperture)
        result = profile_extract(
            image, np.zeros_like(image), np.zeros_like(image),
            valid, aperture, profile,
        )
        self.assertEqual(result["DQ_OPT"][8], 0)
        self.assertAlmostEqual(result["FLUX_OPT"][8], 5)
        self.assertAlmostEqual(result["PROFCOV_OPT"][8], .8)
        self.assertEqual(result["ERR_OPT"][7], 0)


if __name__ == "__main__":
    unittest.main()


class ForwardModelTests(unittest.TestCase):
    def test_exact_fractional_aperture_model_and_chi(self):
        from hstpipe.stis.optimal import profile_extract, profile_fit_diagnostics
        p=np.array([.05,.2,.5,.2,.05]);area=np.array([.3,1,1,1,.8])[:,None]
        amplitude=np.array([-3.,1.,20.])
        # Model is in full-pixel units; aperture overlap enters only once.
        science=p[:,None]*amplitude[None,:]/np.sum(p*area[:,0])
        valid=np.ones(science.shape,bool);err=np.ones_like(science)*.2
        out=profile_extract(science,err,err*0,valid,area,p)
        np.testing.assert_allclose(out['FLUX_OPT'],amplitude,atol=1e-14)
        diag=profile_fit_diagnostics(science,err,valid,area,p,out['FLUX_OPT'])
        np.testing.assert_allclose(diag['MODEL'],science,atol=1e-14)
        np.testing.assert_allclose(diag['CHI2_DIAG'],0,atol=1e-25)
        np.testing.assert_array_equal(diag['NDOF_NOMINAL'],4)
        err[0]=0
        new=profile_extract(science,err,err*0,valid,area,p)
        np.testing.assert_array_equal(new['FLUX_OPT'],out['FLUX_OPT'])
        d=profile_fit_diagnostics(science,err,valid,area,p,new['FLUX_OPT'])
        self.assertEqual(d['ERROR_UNUSABLE_PIXELS'],3)
        self.assertTrue(np.isnan(d['CHI_DIAG'][0]).all())

    def test_profile_excludes_unsupported_blue_edge(self):
        from hstpipe.stis.optimal import empirical_profile
        y,x=np.indices((15,100));wave=np.linspace(1100,1700,100)
        source=np.exp(-.5*((y-7)/1.5)**2)
        area=np.ones_like(source);valid=np.broadcast_to((wave>1150)[None,:],source.shape)
        before=empirical_profile(source,valid,wave,area)
        corrupted=source.copy();corrupted[1,wave<=1150]=1e12
        after=empirical_profile(corrupted,valid,wave,area)
        np.testing.assert_array_equal(before,after)

    def test_display_noise_floor_ignores_source_model_and_keeps_errors(self):
        from hstpipe.stis.optimal import diagnostic_noise_scale
        y,x=np.indices((400,96));spatial=np.arange(400)-200
        rng=np.random.default_rng(43);science=rng.normal(0,2.,y.shape)
        err=np.ones(y.shape)*1e-6;original=err.copy();valid=np.ones(y.shape,bool)
        scaled,floor=diagnostic_noise_scale(science,err,valid,spatial)
        self.assertTrue(np.all((floor>1.8)&(floor<2.2)))
        science[abs(spatial)<80]+=100
        _,same=diagnostic_noise_scale(science,err,valid,spatial)
        np.testing.assert_array_equal(floor,same)
        np.testing.assert_array_equal(err,original)

class WideProfileTests(unittest.TestCase):
    def test_asymmetric_wings_survive_beyond_boxcar(self):
        from hstpipe.stis.optimal import wide_profile, profile_fit_diagnostics
        from hstpipe.stis.extraction import aperture_weights
        y=np.arange(-50,51,dtype=float)
        wave=np.linspace(1400,1700,512)
        # An asymmetric core plus unequal wings: no assumed Gaussian fit.
        truth=np.exp(-.5*((y-1.35)/1.4)**2)
        truth+=.08*np.exp(-.5*((y+12)/4)**2)+.13*np.exp(-.5*((y-17)/5)**2)
        truth/=truth.sum()
        amplitude=30*(1+.2*np.sin(wave/23))
        science=truth[:,None]*amplitude
        valid=np.ones_like(science,bool)
        fit=aperture_weights(y,1,1.35,81)[:,None]
        optimal=aperture_weights(y,1,1.35,61)[:,None]
        box=aperture_weights(y,1,1.35,11)[:,None]
        profile,metadata=wide_profile(science,valid,wave,fit)
        self.assertGreater(profile[y < -5].sum(),.05)
        self.assertGreater(profile[y > 7].sum(),.05)
        out=profile_extract(science,science*0,science*0,valid,optimal,profile)
        expected=np.sum(science*optimal,axis=0)
        np.testing.assert_allclose(out['FLUX_OPT'],expected,rtol=1e-11)
        self.assertGreater(np.mean(out['FLUX_OPT']/np.sum(science*box,axis=0)),1.3)
        diag=profile_fit_diagnostics(science,np.ones_like(science),valid,optimal,profile,out['FLUX_OPT'])
        trained=fit[:,0]>0
        np.testing.assert_allclose(diag['MODEL_FULL'][trained],science[trained],atol=1e-12)
        self.assertTrue(metadata['profile_fitting_independent_of_boxcar'])

    def test_continuum_center_excludes_offset_bright_line(self):
        from hstpipe.stis.optimal import continuum_center
        y=np.arange(-30,31,.5);wave=np.linspace(1140,1725,1024)
        science=np.broadcast_to(np.exp(-.5*((y-1.35)/1.5)**2)[:,None],(len(y),len(wave))).copy()
        science+=100*np.exp(-.5*((y[:,None]+3)/1.5)**2)*np.exp(-.5*((wave[None,:]-1487)/1.2)**2)
        out=continuum_center(science,np.ones_like(science,bool),wave,y*.0246)
        self.assertAlmostEqual(out['shift_native_pixels'],1.35,places=2)
        self.assertLess(out['block_error_native_pixels'],.05)
        self.assertLess(abs(out['excluded_feature_wavelength']-1487),1)

    def test_wide_height_validation(self):
        from hstpipe.stis.config import ExtractionConfig
        for kwargs in ({'optimal_height_native_pixels':5},
                       {'profile_height_native_pixels':20},
                       {'profile_height_native_pixels':float('inf')}):
            with self.assertRaises(ValueError):ExtractionConfig(**kwargs)

class MaskedOptimalTests(unittest.TestCase):
    def test_missing_wing_recovers_full_flux_and_increases_variance(self):
        from hstpipe.stis.optimal import profile_fit_diagnostics
        p=np.array([.02,.08,.2,.4,.2,.08,.02]);flux=np.array([-10.,0.,20.])
        data=p[:,None]*flux;err=np.ones_like(data)*.3;valid=np.ones_like(data,bool)
        before=profile_extract(data,err,err*0,valid,np.ones_like(data),p)
        valid[0]=False;data[0]=np.nan;err[0]=np.nan
        after=profile_extract(data,err,err*0,valid,np.ones_like(data),p)
        np.testing.assert_allclose(after['FLUX_OPT'],flux,atol=1e-14)
        np.testing.assert_allclose(after['PROFCOV_OPT'],.98)
        self.assertTrue(np.all(after['ERR_STAT_OPT']>before['ERR_STAT_OPT']))
        np.testing.assert_allclose(after['ERR_STAT_OPT'],.3/np.sqrt(np.sum(p[1:]**2)))
        diag=profile_fit_diagnostics(data,err,valid,np.ones_like(data),p,after['FLUX_OPT'])
        np.testing.assert_allclose(diag['MODEL_FULL'],p[:,None]*flux,atol=1e-14)
        self.assertTrue(diag['ACCEPTED_COLUMNS'].all())
        self.assertFalse(diag['FULL_COLUMNS'].any())

    def test_missing_core_is_not_fitted_from_tiny_wings(self):
        p=np.array([.02,.08,.2,.4,.2,.08,.02]);data=np.tile(p[:,None],(1,3))
        valid=np.ones_like(data,bool);valid[2:5]=False
        out=profile_extract(data,data*0,data*0,valid,np.ones_like(data),p)
        self.assertTrue(np.isnan(out['FLUX_OPT']).all())
        self.assertTrue((out['DQ_OPT']==4).all())
        self.assertTrue((out['INFOFRAC_OPT']<.5).all())

    def test_fractional_masked_flux_matches_closed_form_and_noise_variance(self):
        rng=np.random.default_rng(7);p=np.array([.04,.16,.4,.3,.1]);area=np.array([.2,1,1,1,.6])[:,None]
        q=p*area[:,0];q/=q.sum();valid=np.ones((5,16000),bool);valid[0]=False
        truth=12.;sigma=.5;science=truth*p[:,None]/np.sum(p*area[:,0])+rng.normal(0,sigma,valid.shape)
        result=profile_extract(science,np.full_like(science,sigma),science*0,valid,area,p)
        # Full-pixel noise propagates through fractional aperture coefficients.
        coeff=q[1:]*area[1:,0]/np.sum(q[1:]**2)
        expected=sigma*np.sqrt(np.sum(coeff**2))
        np.testing.assert_allclose(result['ERR_OPT'],expected)
        self.assertLess(abs(np.mean(result['FLUX_OPT'])-truth),.03)
        self.assertAlmostEqual(np.std(result['FLUX_OPT'])/expected,1,delta=.025)
