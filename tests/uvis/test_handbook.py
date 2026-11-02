"""Focused UVIS DQ/sky policy regressions, with no external catalog network."""
import shutil
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity
from hstpipe.uvis.astrometry.centroids import pixel_masks,validate_centroids

@pytest.mark.parametrize('bits,valid',[(4096,True),(8192,False),(4096|256,False),(4096|2048,False)])
def test_inherited_drizzle_flags_reassessed_but_calwf3_and_saturation_stay_bad(flc,bits,valid):
    p=flc(size=128)
    with fits.open(p,mode='update') as hd:hd['DQ',7].data[64,64]=bits
    before=identity(p);cat=Table(dict(x=[64.],y=[64.],extver=[7],snr=[100.],width_pix=[1.5],minor_width_pix=[1.5]))
    out,audit=validate_centroids(p,cat,{},return_all=True)
    assert bool(len(out))==valid and identity(p)==before
    if valid:assert audit['recovered_from_inherited_cr'][0]

def test_no_blanket_dilation_and_fresh_cr_flags_are_excluded(flc):
    p=flc()
    with fits.open(p,mode='update') as hd:
        hd['DQ',7].data[40,40]=16;hd['DQ',7].data[60,60]=4096
        _,_,mask=pixel_masks(hd,7,{})
        assert mask[40,40] and not mask[40,41] and not mask[60,60]
        hd[0].header['HPCRRUN']=True
        _,_,fresh=pixel_masks(hd,7,{})
        assert fresh[60,60] and fresh[60,61]

def test_resetbits_is_new_pass_only_and_never_clears_calwf3(flc,tmp_path,monkeypatch):
    from hstpipe.uvis.drizzle import drizzle
    from drizzlepac import astrodrizzle
    p=flc();calls=[]
    monkeypatch.setattr(astrodrizzle,'AstroDrizzle',lambda *a,**kw:calls.append(kw))
    drizzle([p,p],tmp_path/'new',dict(driz_cr=True));assert calls[-1]['resetbits']==4096
    with fits.open(p,mode='update') as hd:hd[0].header['HPCRRUN']=True
    drizzle([p,p],tmp_path/'combine',dict(driz_cr=False));assert calls[-1]['resetbits']==0
    assert calls[-1]['skysub'] is False and calls[-1]['skyuser']=='HPSKY'

def test_current_cr_mask_and_provenance_preserve_other_detector_bits(flc,tmp_path):
    from hstpipe.uvis.qa import cr_report
    p=flc(size=128);work=tmp_path/'working.fits';mask=tmp_path/'mask.fits'
    with fits.open(p,mode='update') as hd:
        hd['DQ',7].data[64,64]=4096;hd['DQ',7].data[30,30]=8192;hd['DQ',7].data[31,31]=256
    shutil.copy2(p,work)
    with fits.open(work,mode='update') as hd:
        hd['DQ',7].data[64,64]=0;hd['DQ',7].data[40,40]=4096
    before=identity(p);r=cr_report(p,work,mask,tmp_path,.02)
    assert r[0]['cleared_inherited_4096_fraction']>0 and r[0]['new_cr_fraction']>0
    assert r[0]['calwf3_8192_fraction']>0 and identity(p)==before
    with fits.open(mask) as hd:assert hd['SCI',7].data[40,40]==1 and hd['SCI',7].data[64,64]==0

@pytest.mark.parametrize('filt,mode',[('F275W','none'),('F656N','auto')])
def test_disabled_sky_keeps_measurement_without_double_subtraction(flc,tmp_path,filt,mode):
    from hstpipe.uvis.background import prepare_background,restore_background
    from hstpipe.uvis.config import settings
    p=flc(size=256,filt=filt);paths=ReductionPaths(tmp_path/'redux','TEST','uvis',filt)
    bg=paths.area('intermediate','background',create=True);work=paths.exposure_file('intermediate','working','test','flc')
    before=identity(p);options=settings({})['background'];options['mode']=mode
    result=prepare_background(p,work,bg,options,paths=paths)
    with fits.open(p) as hd,fits.open(work) as wh:
        assert np.array_equal(hd['SCI',7].data,wh['SCI',7].data)
        assert wh['SCI',7].header['HPSKY']==0 and 'HPSKYBG' in wh['SCI',7].header
        assert wh['SCI',7].header['MDRIZSKY']==0
    restored=work.with_name('restored.fits');restore_background(p,restored,paths)
    assert np.array_equal(fits.getdata(restored,('SCI',7)),fits.getdata(work,('SCI',7)))
    assert identity(p)==before

def test_native_resetbits_on_both_working_chips_preserves_other_flags(flc,tmp_path):
    from drizzlepac.resetbits import reset_dq_bits
    p=flc(chips=(1,2));before=identity(p);work=tmp_path/'working_flc.fits';shutil.copy2(p,work)
    with fits.open(work,mode='update') as hd:
        for v in (1,2):hd['DQ',v].data[30,30]=4096|8192|256
    from hstpipe.uvis.dq import native_compatible_dq
    with fits.open(work,mode='update') as hd:native_compatible_dq(hd)
    reset_dq_bits(str(work),4096,extname='DQ')
    with fits.open(work) as hd:
        for v in (1,2):assert hd['DQ',v].data[30,30]==8192|256
    assert identity(p)==before

def test_no_cr_mosaic_detection_measures_filtered_noise(flc):
    from hstpipe.uvis.astrometry.sources import sources
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:
        rng=np.random.default_rng(113);a=rng.normal(0,1,(256,256))
        a[rng.random(a.shape)<.035]+=80
        yy,xx=np.indices(a.shape);a+=600*np.exp(-((xx-128)**2+(yy-128)**2)/8)
        hd['SCI',7].data=a.astype('f4')
        del hd['ERR',7]
        hd.append(fits.ImageHDU(np.ones(a.shape,'f4'),name='WHT'))
    before=identity(p);cat=sources(p,{'snr':5,'min_area':9},True,6)
    measured=cat.meta['detection_chips'][0]
    assert measured['filtered_noise']>0
    assert measured['detection_threshold']>=5*measured['filtered_noise']
    assert 'correlated noise' in measured['detection_noise_method']
    assert len(cat) and np.min(np.hypot(cat['x']-128,cat['y']-128))<1
    assert identity(p)==before
