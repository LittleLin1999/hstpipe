import json,shutil
from pathlib import Path
import numpy as np
import pytest
import yaml
from astropy.io import fits
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.config import ReductionConfig
from hstpipe.common.cache import StageCache
from hstpipe.common.provenance import identity,write_manifest
from hstpipe.common.runner import run
from hstpipe.uvis.config import settings
from hstpipe.uvis.inventory import discover
from hstpipe.uvis.grouping import group_exposures
from hstpipe.uvis.background import prepare_background,restore_background
from hstpipe.uvis.migration import migrate_uvis


def test_target_paths_and_names():
    a=ReductionPaths(Path('/data/redux'),'T','uvis','F275W')
    b=ReductionPaths(a.root,'T','uvis','F336W')
    c=ReductionPaths(a.root,'T','uvis','F336W','field-abc')
    assert a.directory==b.directory==c.directory==Path('/data/redux/uvis/T')
    assert a.product('drc').name=='T_F275W_drc.fits'
    assert c.product('drc').name=='T_F336W_field-abc_drc.fits'
    assert a.file('diagnostics','drizzle','coverage','png')!=b.file('diagnostics','drizzle','coverage','png')
    assert a.cache_directory==b.cache_directory
    assert a.cache().receipt('reduction')!=b.cache().receipt('reduction')
    assert c.exposure_file('intermediate','background','root','sky').name=='root_field-abc_sky.fits'
    stis=ReductionPaths(a.root,'T','stis','G140L')
    assert stis.directory==Path('/data/redux/T/stis/G140L') and stis.product('x1d').name=='T_G140L_x1d.fits'
    assert 'filter_overrides' not in ReductionConfig('/raw','/redux','stis').effective()


def test_filter_override_and_independent_rerun(flc,tmp_path):
    raw=tmp_path/'raw';raw.mkdir()
    for name,filt in [('a','F275W'),('b','F336W')]:
        p=flc(name,filt=filt);p.rename(raw/p.name)
    calls=[]
    def engine(config,paths,rows):
        calls.append((config.element,config.options['drizzle']['final_pixfrac']))
        p=paths.product('drc');p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(str(calls[-1]).encode());return [p]
    cfg=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',resume=True,filter_overrides={'F336W':{'drizzle':{'final_pixfrac':.7}}})
    [out]=run(cfg,engine);assert len(calls)==2 and calls[1][1]==.7
    second=identity(out/'products/TEST_F336W_drc.fits');record=json.loads((out/'manifest.json').read_text())['filters']['F336W']
    cfg.element='F275W';cfg.options={'drizzle':{'final_scale':.03}}
    run(cfg,engine)
    assert len(calls)==3 and identity(out/'products/TEST_F336W_drc.fits')==second
    assert json.loads((out/'manifest.json').read_text())['filters']['F336W']==record
    assert not (out/'F275W').exists()
    assert yaml.safe_load((out/'config.yaml').read_text())['element'] is None
    cfg.force=True;run(cfg,engine)
    assert identity(out/'products/TEST_F336W_drc.fits')==second


def test_stable_setup_ids(flc,tmp_path):
    flc('a');flc('b',ra=11)
    rows=discover(tmp_path);first=group_exposures(rows);second=group_exposures(list(reversed(rows)))
    assert set(first)==set(second) and all(k[2].startswith('field-') for k in first)
    assert len({ReductionPaths(tmp_path,'TEST','uvis',k[1],k[2]).product('drc') for k in first})==2


def test_background_checkpoint_replay_is_numerically_identical(flc,tmp_path):
    p=flc(size=256);paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    bg=paths.area('intermediate','background',create=True);paths.area('diagnostics','background',create=True)
    first=tmp_path/'first.fits';again=tmp_path/'again.fits'
    opts=settings({'background':{'target_radius':.2,'box_size':32}})['background']
    prepare_background(p,first,bg,opts,paths=paths);restore_background(p,again,paths)
    with fits.open(first) as a,fits.open(again) as b:
        for name in ('SCI','ERR','DQ'):assert np.array_equal(a[name,7].data,b[name,7].data)
        for key in ('BUNIT','CRVAL1','CRVAL2','PHOTFLAM','PHOTPLAM','HPSKY','HPSKYERR'):
            assert a['SCI',7].header[key]==b['SCI',7].header[key]


def legacy_group(flc,raw,out,name,filt,pixfrac=1):
    p=flc(name,filt=filt);p.rename(raw/p.name);p=raw/p.name
    d=out/'TEST/uvis'/filt
    for stage in ('background','alignment/prepared','drizzle'):(d/'intermediate'/stage).mkdir(parents=True)
    options=settings({'astrometry':{'method':'relative_only'},'background':{'mode':'none'},'alignment':{'mode':'original'},'drizzle':{'final_pixfrac':pixfrac}})
    cfg=ReductionConfig(str(raw),str(out),'uvis',options=options,resume=True)
    prepared=d/'intermediate/alignment/prepared'/p.name
    background=prepare_background(p,prepared,d/'intermediate/background',options['background'])
    alignment=dict(status='original_wcs',reference=str(prepared),fits=[],warnings=[])
    upstream=d/'intermediate/alignment/upstream.json';upstream.write_text(json.dumps(dict(exposures={},background={p.name:background},alignment=alignment)))
    mask=d/'intermediate/drizzle'/f'{name}_crmask.fits'
    with fits.open(p) as hd:
        mh=fits.ImageHDU(np.zeros(hd['SCI',7].shape,'u1'),name='SCI');mh.header['EXTVER']=7
        fits.HDUList([fits.PrimaryHDU(),mh]).writeto(mask)
        hdr=hd['SCI',7].header.copy();hdr['BUNIT']='ELECTRONS/S'
        product=d/f'products/TEST_{filt}_drc.fits';product.parent.mkdir()
        fits.HDUList([fits.PrimaryHDU(header=hd[0].header),fits.ImageHDU(hd['SCI',7].data/100,header=hdr,name='SCI'),
                      fits.ImageHDU(np.ones(mh.shape),name='WHT'),fits.ImageHDU(np.ones(mh.shape,'i4'),name='CTX')]).writeto(product,checksum=True)
    cr=d/'intermediate/drizzle/cr.json';cr.write_text(json.dumps({p.name:[]}))
    qa=d/'diagnostics/drizzle/qa_summary.json';qa.parent.mkdir(parents=True)
    qa.write_text(json.dumps(dict(alignment=alignment,status='review_required')))
    (d/'config.yaml').write_text(yaml.safe_dump(cfg.effective()))
    manifest=d/'manifest.json'
    write_manifest(manifest,status='complete',inputs=[identity(p)],configuration=cfg.effective(),products=[identity(product)],scientific_assessment=dict(inventory={p.name:{}}))
    cache=StageCache(d/'.cache')
    cache.record('uvis_upstream','old-upstream',[prepared,upstream]+list((d/'intermediate/background').glob('*.fits')))
    cache.record('uvis_cr','old-cr',[mask,cr])
    cache.record('reduction','old-reduction',[q for q in d.rglob('*') if q.is_file() and '.cache' not in q.parts])
    return d


def test_migration_bytes_cache_config_and_recovery(flc,tmp_path):
    raw=tmp_path/'raw';raw.mkdir();out=tmp_path/'redux'
    a=legacy_group(flc,raw,out,'a','F275W');b=legacy_group(flc,raw,out,'b','F336W',.7)
    before={p.name:identity(p)['sha256'] for d in (a,b) for p in d.rglob('*.fits')}
    plans=migrate_uvis(out,dry_run=True);assert plans[0]['status']=='planned' and a.exists()
    [report]=migrate_uvis(out);d=out/'TEST/uvis';backup=Path(report['backup'])
    assert backup.exists() and not (d/'F275W').exists() and not (d/'F336W').exists()
    for p in d.rglob('*.fits'):
        if p.name in before:assert identity(p)['sha256']==before[p.name]
    cfg=yaml.safe_load((d/'config.yaml').read_text())
    assert cfg['filter_overrides']['F336W']['drizzle']['final_pixfrac']==.7
    state=json.loads((d/'manifest.json').read_text());assert set(state['filters'])=={'F275W','F336W'}
    for f in ('F275W','F336W'):
        receipt=StageCache(d/'.cache'/f'TEST_{f}')
        for stage in ('reduction','uvis_background','uvis_alignment','uvis_cr'):
            key=json.loads((receipt.directory/f'{stage}.json').read_text())['key'];assert receipt.valid(stage,key)
        assert (d/f'intermediate/alignment/TEST_{f}_registration.ecsv').exists()
    # Resume each migrated filter must not call any scientific processing.
    def forbidden(*args):raise AssertionError('Migration unnecessarily reprocessed a valid product')
    run(ReductionConfig(**cfg),forbidden)
    assert identity(backup/'F275W/products/TEST_F275W_drc.fits')['sha256']==before['TEST_F275W_drc.fits']
    # Replacing a prepared file in the new layout cannot alter its recovery link.
    prepared=d/'intermediate/alignment/a_flc.fits';saved=identity(backup/'F275W/intermediate/alignment/prepared/a_flc.fits')
    paths=ReductionPaths(out,'TEST','uvis','F275W');paths.detach(prepared);shutil.copy2(raw/'a_flc.fits',prepared)
    assert identity(backup/'F275W/intermediate/alignment/prepared/a_flc.fits')==saved


def test_migration_conflict_is_non_destructive(flc,tmp_path):
    raw=tmp_path/'raw';raw.mkdir();out=tmp_path/'redux';old=legacy_group(flc,raw,out,'a','F275W')
    product=out/'TEST/uvis/products/TEST_F275W_drc.fits';product.parent.mkdir();product.write_bytes(b'conflict')
    with pytest.raises(FileExistsError):migrate_uvis(out)
    assert old.exists() and product.read_bytes()==b'conflict'


def test_stage_invalidation_is_local(flc,tmp_path,monkeypatch):
    from hstpipe.uvis import pipeline
    raw=tmp_path/'raw';raw.mkdir();p=flc();p.rename(raw/p.name)
    calls={'background':0,'alignment':0,'cr':0,'final':0}
    from hstpipe.uvis.astrometry import relative_engine
    original_bg=pipeline.prepare_background;original_align=relative_engine.refine_relative;original_cr=pipeline.cr_report
    def bg(*a,**k):calls['background']+=1;return original_bg(*a,**k)
    def al(*a,**k):calls['alignment']+=1;return original_align(*a,**k)
    def cr(*a,**k):calls['cr']+=1;return original_cr(*a,**k)
    def drizzle(inputs,prefix,options):
        calls['final']+=1
        with fits.open(inputs[0]) as hd:
            h=hd['SCI',7].header.copy();h['BUNIT']='ELECTRONS/S';shape=hd['SCI',7].shape
            fits.HDUList([fits.PrimaryHDU(header=hd[0].header.copy()),fits.ImageHDU(hd['SCI',7].data/100,header=h,name='SCI'),
                fits.ImageHDU(np.ones(shape),name='WHT'),fits.ImageHDU(np.ones(shape,'i4'),name='CTX')]).writeto(str(prefix)+'_drc.fits')
    monkeypatch.setattr(relative_engine,'refine_relative',al)
    for name,value in [('prepare_background',bg),('cr_report',cr),('drizzle',drizzle)]:monkeypatch.setattr(pipeline,name,value)
    cfg=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',resume=True,options={'alignment':{'mode':'original'},'astrometry':{'method':'relative_only'}})
    pipeline.reduce_uvis(cfg);assert calls==dict(background=1,alignment=1,cr=1,final=1)
    cfg.options['drizzle']['final_scale']=.03;pipeline.reduce_uvis(cfg)
    assert calls==dict(background=1,alignment=1,cr=1,final=2)
    cfg.options['drizzle']['driz_cr_snr']='5 4';pipeline.reduce_uvis(cfg)
    assert calls==dict(background=1,alignment=1,cr=2,final=3)
    cfg.options['alignment']['max_rms']=.4;pipeline.reduce_uvis(cfg)
    assert calls['background']==1 and calls['alignment']==2 and calls['final']==4
    cfg.options['background']['threshold']=3;pipeline.reduce_uvis(cfg)
    assert calls['background']==2 and calls['alignment']==3


def test_top_level_astrometry_yaml_is_normalized(tmp_path):
    path=tmp_path/'config.yaml';path.write_text('input: raw/uvis\noutput: redux\ninstrument: uvis\nastrometry:\n  method: two_pass\n')
    cfg=ReductionConfig.read(path)
    assert cfg.options['astrometry']['method']=='two_pass'


@pytest.mark.parametrize('disjoint',[False,True])
def test_failed_selected_filter_can_retry_later_validated_anchor(flc,tmp_path,disjoint):
    raw=tmp_path/'raw';raw.mkdir()
    for name,filt in [('a','F467M'),('b','F689M'),('c','F275W')]:
        p=flc(name,filt=filt)
        if filt=='F467M':
            with fits.open(p,mode='update') as hd:hd[0].header['EXPTIME']=200.
        p.rename(raw/p.name)
    if disjoint:
        p=flc('d',filt='F467M',ra=11)
        with fits.open(p,mode='update') as hd:hd[0].header['EXPTIME']=200.
        p.rename(raw/p.name)
    calls=[]
    def engine(config,paths,rows):
        calls.append(config.element)
        state=json.loads((paths.directory/'manifest.json').read_text())
        anchor=any(v.get('products') for v in state['filters'].values())
        ok=config.element=='F689M' or anchor
        qa=paths.file('diagnostics','drizzle','qa_summary')
        save=dict(status='complete' if ok else 'astrometry_unvalidated',
            astrometry=dict(status='validated' if ok else 'unvalidated',attempts=[{'reference':'legacy_dr10'}]))
        qa.parent.mkdir(parents=True,exist_ok=True);qa.write_text(json.dumps(save))
        if not ok:return []
        product=paths.product('drc');product.parent.mkdir(parents=True,exist_ok=True);product.write_bytes(config.element.encode());return [product]
    cfg=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',resume=True,
        options=settings({'astrometry':{'method':'two_pass'}}))
    [out]=run(cfg,engine)
    assert calls==(['F467M','F467M','F275W','F689M','F467M','F467M','F275W'] if disjoint else ['F467M','F275W','F689M','F467M','F275W'])
    state=json.loads((out/'manifest.json').read_text())
    assert all(r['scientific_assessment']['astrometry']['status']=='validated' for r in state['filters'].values())
    calls.clear();cfg.element='F275W';run(cfg,engine)
    assert calls==[]
