"""Focused maintenance checks: selection, configuration and cache boundaries."""
import numpy as np
import pytest
from astropy.table import Table
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry.matching import detect_sources


def test_unknown_detector_does_not_silently_fall_back_to_dao():
    with pytest.raises(ValueError,match='Unknown'):
        detect_sources('not_opened.fits',{'method':'typo'})
    with pytest.raises(ValueError,match='Unknown'):
        settings({'astrometry':{'detection':{'method':'typo'}}})


def test_compact_settings_validate_without_changing_relative_policy():
    old=settings({});new=settings({'astrometry':{'detection':{'method':'compact_stable'}}})
    assert old['astrometry']['relative']==new['astrometry']['relative']
    assert new['drizzle']['final_pixfrac']==.8
    for value in ([6.,3.5],[3.5,np.inf]):
        with pytest.raises(ValueError,match='kernels'):
            settings({'astrometry':{'detection':{'method':'compact_stable','compact_kernel_widths':value}}})


def test_compact_selection_respects_catalog_cap(monkeypatch):
    from hstpipe.uvis.astrometry import sources
    t=Table(dict(x=[10.,20.],y=[10.,20.],extver=[1,1],snr=[20.,20.],
        centroid_valid=[True,True],selected_for_astrometry=[True,False],
        quality_flags=['',''],rejection_reason=['',''],centroid_status=['accepted','accepted']))
    monkeypatch.setattr(sources,'sources',lambda *a,**k:(None,t.copy()))
    clean,audit=detect_sources('unused',{'method':'compact_stable'},return_all=True)
    assert len(clean)==1 and list(audit['selected_for_astrometry'])==[True,False]
    assert clean.meta['compact_selection']['pixel_quality_accepted']==2
    assert clean.meta['compact_selection']['accepted']==1


def test_report_package_not_a_scientific_cache_dependency(tmp_path,monkeypatch):
    from hstpipe.uvis import storage
    from hstpipe.common.cache import StageCache
    root=tmp_path/'uvis';(root/'reporting').mkdir(parents=True);(root/'astrometry').mkdir()
    science=root/'astrometry'/'matching.py';science.write_text('science v1')
    report=root/'reporting'/'report.py';report.write_text('report v1')
    (root/'report_decisions.py').write_text('old presentation')
    (root/'._metadata.py').write_bytes(b'not Python')
    monkeypatch.setattr(storage,'__file__',str(root/'storage.py'))
    common=tmp_path/'common';common.mkdir();calibration=common/'calibration.py';calibration.write_text('calibration v1')
    cache=StageCache(tmp_path/'cache');first=cache.key(storage.code_inputs(),{})
    assert storage.code_inputs()==[calibration,science]
    report.write_text('report v2');assert first==cache.key(storage.code_inputs(),{})
    science.write_text('science v2');assert first!=cache.key(storage.code_inputs(),{})


def test_absolute_selection_options_do_not_invalidate_relative_or_background(flc,tmp_path):
    from hstpipe.uvis.storage import stage_keys
    from hstpipe.uvis.inventory import discover
    from hstpipe.common.paths import ReductionPaths
    flc();rows=discover(tmp_path);cache=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W').cache()
    before=settings({});after=settings({'astrometry':{'detection':{'method':'compact_stable'},'catalogs':{'legacy_dr10':{'quality':{'min_snr':0,'max_position_error_mas':200}}}}})
    assert stage_keys(cache,rows,before,[])[:2]==stage_keys(cache,rows,after,[])[:2]


def test_compact_report_uses_configuration_without_manual_resume_field(tmp_path):
    import json
    from hstpipe.uvis.reporting import generate_report
    p=tmp_path/'T'/'uvis';p.mkdir(parents=True)
    record={'status':'astrometry_unvalidated','filter':'F689M','configuration':{'options':settings({'astrometry':{'detection':{'method':'compact_stable','compact_min_snr':12.}}})},'scientific_assessment':{}}
    (p/'manifest.json').write_text(json.dumps({'target':'T','filters':{'F689M':record}}))
    generate_report(p)
    report=p/'index.html'
    text=report.read_text();assert 'Absolute reference selection:' in text and '12.0' in text
    assert report.stat().st_size<1000000
