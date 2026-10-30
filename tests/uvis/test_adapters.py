import sys,types
from hstpipe.uvis.drizzle import drizzle
from hstpipe.uvis.alignment import align

def test_official_api_dispatch(monkeypatch,tmp_path):
    calls={}
    fake=types.SimpleNamespace(astrodrizzle=types.SimpleNamespace(AstroDrizzle=lambda *a,**k:calls.update(drizzle=(a,k))),tweakreg=types.SimpleNamespace(TweakReg=lambda *a,**k:calls.update(align=(a,k))))
    monkeypatch.setitem(sys.modules,'drizzlepac',fake)
    paths=[tmp_path/'a_flc.fits',tmp_path/'b_flc.fits'];align(paths);drizzle(paths,tmp_path/'result')
    assert calls['align'][1]['updatehdr'] is True
    assert calls['drizzle'][1]['driz_cr'] is True
    assert calls['drizzle'][1]['final_wht_type']=='IVM'

def test_signed_sky_assessment_does_not_modify_flc(tmp_path):
    from astropy.io import fits
    import numpy as np
    from hstpipe.uvis.background import assess_regions
    from hstpipe.common.provenance import identity
    p=tmp_path/'flc.fits';m=tmp_path/'mask.fits'
    shape=(10,12);values=np.arange(120).reshape(shape)/100-2
    hd=[fits.PrimaryHDU()]
    for name,array in [('SCI',values),('ERR',np.ones(shape)),('DQ',np.zeros(shape,'u2'))]:
        h=fits.ImageHDU(array,name=name);h.header['EXTVER']=1;hd.append(h)
    fits.HDUList(hd).writeto(p)
    fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(np.ones(shape,'u2'),name='SCI')]).writeto(m)
    before=identity(p);report=assess_regions(p,m)
    assert report[0]['mean']<0
    assert np.isclose(report[0]['mean'],values.mean())
    assert identity(p)==before
