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
