"""Focused compact-reference selection checks; no survey or reduction jobs."""
import numpy as np
from astropy.table import Table
from hstpipe.uvis.astrometry.matching import detect_sources
from hstpipe.uvis.astrometry.catalogs import standardize_legacy


def catalog(x, valid=None, ext=None, snr=None):
    n=len(x)
    return Table(dict(x=x,y=np.zeros(n),extver=ext or [1]*n,snr=snr or [20.]*n,
        centroid_valid=valid or [True]*n,quality_flags=['']*n,rejection_reason=['']*n,
        centroid_status=['accepted']*n,selected_for_astrometry=[True]*n))


def test_both_kernel_quality_and_stability_required(monkeypatch):
    from hstpipe.uvis.astrometry import sources
    narrow=catalog([10.,20.,30.,40.],valid=[True,False,True,True],snr=[20.,20.,20.,8.])
    wide=catalog([10.1,20.1,30.8,40.1])
    monkeypatch.setattr(sources,'sources',lambda *a,width,**k:(None,(narrow if width==3.5 else wide).copy()))
    clean,audit=detect_sources('unused',dict(method='compact_stable'),return_all=True)
    assert list(clean['x'])==[10.1]
    assert list(audit['selected_for_astrometry'])==[True,False,False,False]
    assert 'quality_failed' in audit['rejection_reason'][1]
    assert 'missing_or_ambiguous' in audit['rejection_reason'][2]
    assert audit['rejection_reason'][3]=='compact_low_hst_snr'
    assert not clean.meta['compact_selection']['external_positions_used']


def test_chip_identity_and_ambiguity(monkeypatch):
    from hstpipe.uvis.astrometry import sources
    narrow=catalog([10.,20.,20.2],ext=[2,1,1]);wide=catalog([10.,20.1])
    monkeypatch.setattr(sources,'sources',lambda *a,width,**k:(None,(narrow if width==3.5 else wide).copy()))
    assert len(detect_sources('unused',dict(method='compact_stable')))==0


def test_precision_selection_rejects_missing_or_large_errors():
    n=3
    raw=Table(dict(ls_id=np.arange(n),ra=10+np.arange(n)*.01,dec=np.zeros(n),type=['REX']*n,
        brick_primary=[True]*n,maskbits=[0]*n,flux_r=[9.5]*n,flux_z=[9.5]*n,
        flux_ivar_r=[1.]*n,flux_ivar_z=[1.]*n,ra_ivar=[(3600000/60)**2,(3600000/300)**2,0.],
        dec_ivar=[(3600000/60)**2,(3600000/300)**2,0.]))
    assert len(standardize_legacy(raw))==0
    out=standardize_legacy(raw,dict(min_snr=0,max_position_error_mas=200))
    assert list(out['source_id'])==['0']
    assert out.meta['position_error_rejected']==2
