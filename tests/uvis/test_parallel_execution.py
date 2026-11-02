"""Focused process-isolation and shared-target transaction checks."""
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from dataclasses import replace
from pathlib import Path
import time,json
from hstpipe.common.config import ReductionConfig
from hstpipe.common.target_runner import run_uvis

def _overlapping_engine(config,paths,rows):
    marker=paths.file('intermediate','alignment','parallel_marker','txt')
    marker.parent.mkdir(parents=True,exist_ok=True);marker.write_text('running')
    deadline=time.monotonic()+15
    while len(list(marker.parent.glob('*_parallel_marker.txt')))<2:
        if time.monotonic()>deadline:raise RuntimeError('Filter computations were serialized')
        time.sleep(.05)
    product=paths.file('products',None,'parallel_test','txt');product.parent.mkdir(exist_ok=True);product.write_text(config.element)
    return [product]

def _run(config):return run_uvis(config,_overlapping_engine,_retry_anchors=False)

def test_filters_compute_concurrently_without_lost_manifest_rows(flc,tmp_path):
    from hstpipe.common.provenance import identity
    a=flc(name='first',filt='F275W');b=flc(name='second',filt='F336W');raw=tmp_path/'raw';raw.mkdir();a=a.rename(raw/a.name);b=b.rename(raw/b.name);before=[identity(a),identity(b)]
    config=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',target='TEST',resume=True)
    with ProcessPoolExecutor(max_workers=2,mp_context=get_context('spawn')) as pool:
        jobs=[pool.submit(_run,replace(config,element=f)) for f in ('F275W','F336W')]
        for job in jobs:job.result(timeout=40)
    root=tmp_path/'redux/TEST/uvis';state=json.loads((root/'manifest.json').read_text())
    assert set(state['filters'])=={'F275W','F336W'}
    assert all(r['status']=='complete' for r in state['filters'].values())
    assert [identity(a),identity(b)]==before
    assert (root/'products/TEST_uvis_report.html').stat().st_size<1_000_000

def test_nested_catalog_cache_lock_is_reentrant_and_reuses_raw(tmp_path):
    from astropy.table import Table
    from hstpipe.uvis.astrometry.catalogs import _cached
    calls=[]
    def raw():calls.append(True);return Table({'ra':[10.],'dec':[20.]}),{}
    def normalized():return _cached('raw',{},raw,tmp_path,False,False),{}
    first=_cached('normalized',{},normalized,tmp_path,False,False)
    second=_cached('normalized',{},normalized,tmp_path,False,True)
    assert len(first)==len(second)==1 and len(calls)==1
