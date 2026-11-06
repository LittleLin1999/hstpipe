import json
import numpy as np
from pathlib import Path
from astropy.io import fits
from hstpipe.common.config import ReductionConfig
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.inventory import discover
from hstpipe.common.cache import StageCache
import pytest

def exposure(path,instrument='WFC3',detector='UVIS',subarray=False):
    header=fits.Header(dict(INSTRUME=instrument,DETECTOR=detector,TARGNAME='TEST',FILTER='F606W',OPT_ELEM='G140L',OBSTYPE='IMAGING',RA_TARG=10.,DEC_TARG=20.,ROOTNAME=path.stem,BIASCORR='COMPLETE',BLEVCORR='COMPLETE',DARKCORR='COMPLETE',FLATCORR='COMPLETE',PCTECORR='COMPLETE',SUBARRAY=subarray,APERTURE='UVIS',EXPTIME=100.))
    hd=[fits.PrimaryHDU(header=header)]
    for name in ('SCI','ERR','DQ'):
        h=fits.ImageHDU(np.ones((8,9),dtype='u2' if name=='DQ' else 'f4'),name=name);h.header['EXTVER']=1
        if name=='SCI':
            from astropy.wcs import WCS
            w=WCS(naxis=2);w.wcs.crpix=[4,4];w.wcs.cdelt=[-0.04/3600,0.04/3600];w.wcs.crval=[10,20];w.wcs.ctype=['RA---TAN','DEC--TAN'];h.header.update(w.to_header());h.header['BUNIT']='ELECTRONS'
        hd.append(h)
    fits.HDUList(hd).writeto(path)

def test_paths(tmp_path):
    p=ReductionPaths(tmp_path,'J1022+0841','stis','G140L')
    assert p.product('x1d')==tmp_path/'J1022+0841/stis/G140L/products/J1022+0841_G140L_x1d.fits'
    assert not p.directory.exists()
    with pytest.raises(ValueError):ReductionPaths(tmp_path,'..','uvis','F606W').directory

def test_cache_invalidates(tmp_path):
    source=tmp_path/'input';source.write_text('first');output=tmp_path/'out';output.write_text('result')
    c=StageCache(tmp_path/'cache');key=c.key([source],{'mode':1});c.record('step',key,[output]);assert c.valid('step',key)
    assert c.key([source],{'mode':2})!=key
    source.write_text('changed');assert c.key([source],{'mode':1})!=key
    output.write_text('stale');assert not c.valid('step',key)

def test_uvis_inventory_and_dry_run(tmp_path):
    from hstpipe.uvis import reduce_uvis
    raw=tmp_path/'raw';raw.mkdir();exposure(raw/'a_flc.fits',subarray=True)
    rows=discover(raw,'uvis');assert rows[0].metadata['SUBARRAY'];assert rows[0].metadata['imsets']==[1]
    c=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis')
    plan=reduce_uvis(c,dry_run=True);assert len(plan)==1;assert not (tmp_path/'redux').exists()

def test_duplicate_detection(tmp_path):
    import shutil
    exposure(tmp_path/'a_flc.fits');shutil.copy2(tmp_path/'a_flc.fits',tmp_path/'b_flc.fits')
    with pytest.raises(ValueError,match='Duplicate'):discover(tmp_path)

def test_config_paths_and_raw_protection(tmp_path):
    path=tmp_path/'config.yaml';path.write_text('input: raw/stis\noutput: redux\ninstrument: stis\n')
    c=ReductionConfig.read(path);assert c.input==str(tmp_path/'raw/stis')
    with pytest.raises(ValueError):ReductionConfig(str(tmp_path),str(tmp_path/'redux'),'stis')

def test_constant_hst_imset(tmp_path):
    from hstpipe.common.fits import imsets
    path=tmp_path/'acquisition_raw.fits'
    hd=[fits.PrimaryHDU(),fits.ImageHDU(np.ones((7,8)),name='SCI')]
    hd[1].header['EXTVER']=1
    for name in ('ERR','DQ'):
        h=fits.ImageHDU(name=name)
        for key,value in dict(EXTVER=1,PIXVALUE=0,NPIX1=8,NPIX2=7).items():h.header[key]=value
        hd.append(h)
    fits.HDUList(hd).writeto(path)
    assert imsets(path)==[1]

def test_resume_in_existing_directory_and_cache(tmp_path):
    from hstpipe.common.runner import run
    raw=tmp_path/'raw';raw.mkdir();exposure(raw/'a_flc.fits')
    c=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',resume=True)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F606W')
    paths.directory.mkdir(parents=True)
    calls=[]
    def engine(config,paths,rows):
        calls.append(True)
        p=paths.area('products',create=True)/'test.txt';p.write_text('complete')
        return [p]
    run(c,engine)
    assert len(calls)==1
    run(c,engine)
    assert len(calls)==1

def test_force_rebuilds_even_with_resume(tmp_path):
    from hstpipe.common.runner import run
    raw=tmp_path/'raw';raw.mkdir();exposure(raw/'a_flc.fits')
    c=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',resume=True)
    calls=[]
    def engine(config,paths,rows):
        calls.append(True);p=paths.area('products',create=True)/'result.txt';p.write_text(str(len(calls)));return [p]
    run(c,engine);c.force=True;run(c,engine)
    assert len(calls)==2
    assert list((tmp_path/'redux/uvis').glob('.uvis-TEST-history-TEST_F606W-*'))
