"""Focused multi-chip, shallow paths, recovery, cache and standalone-report checks."""
import json,shutil
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.cache import StageCache
from hstpipe.common.config import ReductionConfig
from hstpipe.common.provenance import identity
from hstpipe.common.target_runner import save,run_uvis
from hstpipe.uvis.config import settings

def test_shallow_namespaces_and_unchanged_stis(tmp_path):
    a=ReductionPaths(tmp_path,'T','uvis','F275W');b=ReductionPaths(tmp_path,'T','uvis','F336W')
    assert a.catalog_directory()==b.catalog_directory()==a.directory/'intermediate/catalogs'
    assert a.area('diagnostics','alignment')==a.area('diagnostics','drizzle')==a.directory/'diagnostics'
    assert a.cache().receipt('uvis_sources')!=b.cache().receipt('uvis_sources')
    assert a.exposure_file('intermediate','alignment','root','flc').parent==a.directory/'intermediate/working'
    with pytest.raises(ValueError):a.area('intermediate','alignment/filters')
    assert not a.directory.exists()
    stis=ReductionPaths(tmp_path,'T','stis','G140L')
    assert stis.directory==tmp_path/'T/stis/G140L' and stis.cache().receipt('reduction')==stis.directory/'.cache/reduction.json'

def test_recoverable_migration_preserves_fits_and_background_cache(flc,tmp_path):
    from hstpipe.uvis.shallow import migrate_target
    from hstpipe.uvis.inventory import discover
    from hstpipe.uvis.storage import stage_keys
    raw=tmp_path/'raw';raw.mkdir();source=flc();source.rename(raw/source.name);source=raw/source.name
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');old=paths.directory/'intermediate/alignment/test_flc.fits';old.parent.mkdir(parents=True)
    shutil.copy2(source,old);before=identity(old)
    bg=paths.directory/'intermediate/background/TEST_F275W_background.json';save(bg,{'exposures':{},'background':{},'working':str(old)})
    opts=settings({});rows=discover(raw);cfg=ReductionConfig(str(raw),str(paths.root),'uvis',options=opts,resume=True)
    key=stage_keys(paths.cache(),rows,opts,[old])[0]
    StageCache(paths.directory/'.cache'/paths.group_id).record('uvis_background',key,[bg])
    nested=paths.directory/'diagnostics/alignment/TEST_F275W_check.json';save(nested,{'input':str(old)})
    collision=paths.directory/'diagnostics/background/TEST_F275W_check.json';save(collision,{'background':'different figure'})
    save(paths.directory/'manifest.json',dict(schema_version=2,target='TEST',instrument='uvis',filters={'F275W':dict(filter='F275W',configuration=cfg.effective(),inputs=[identity(source)],products=[],artifacts=[identity(old)],scientific_assessment={})}))
    result=migrate_target(paths.directory);new=paths.exposure_file('intermediate','working','test','flc')
    assert identity(new)['sha256']==before['sha256'] and identity(Path(result['backup'])/'intermediate/alignment/test_flc.fits')['sha256']==before['sha256']
    assert paths.cache().valid('uvis_background',key)
    assert not any(p.is_dir() for p in (paths.directory/'intermediate').glob('*/*'))
    assert json.loads((paths.directory/'diagnostics/alignment_TEST_F275W_check.json').read_text())['input']==str(new)
    assert json.loads((paths.directory/'diagnostics/background_TEST_F275W_check.json').read_text())['background']=='different figure'
    assert migrate_target(paths.directory)['status']=='already_shallow'

def test_covering_catalog_shared_across_filters_and_epochs(tmp_path,monkeypatch):
    from hstpipe.uvis.astrometry import catalogs as c
    raw=Table(dict(source_id=[123456789012345678],ra=[10.],dec=[20.],ra_error=[1.],dec_error=[1.],pmra=[2.],pmdec=[3.],pmra_error=[1.],pmdec_error=[1.],parallax=[-1.],parallax_error=[1.],ref_epoch=[2016.],phot_g_mean_mag=[18.],ruwe=[1.],astrometric_params_solved=[31],visibility_periods_used=[12],duplicated_source=[False]))
    calls=[]
    def request(*a,**kw):calls.append(a);return raw
    monkeypatch.setattr(c,'request',request)
    a=ReductionPaths(tmp_path,'T','uvis','F275W');b=ReductionPaths(tmp_path,'T','uvis','F336W')
    c.query_gaia_dr3(10,20,.01,2020.,cache_dir=a.catalog_directory())
    t=c.cached_footprint('gaia_dr3',10,20,.005,cache_dir=b.catalog_directory(),observation_epoch=2026.)
    assert len(calls)==1 and t.meta['complete_query_coverage'] and t['position_epoch'][0]==2026.
    assert not any(p.is_dir() for p in a.catalog_directory().iterdir())
    incompatible=c.cached_footprint('gaia_dr3',10,20,.02,cache_dir=b.catalog_directory(),observation_epoch=2026.)
    assert not incompatible.meta['complete_query_coverage']

@pytest.mark.parametrize('status',['failed','relative_astrometry_unvalidated','complete'])
def test_report_after_success_or_failure_is_offline_and_bounded(tmp_path,status):
    from hstpipe.uvis.reporting import generate_report
    paths=ReductionPaths(tmp_path,'T','uvis','F275W')
    save(paths.directory/'manifest.json',dict(target='T',filters={'F275W':dict(status=status,error='masked source core' if status=='failed' else '',products=[],inputs=[])}))
    p=generate_report(paths.directory);data=p.read_bytes()
    assert 0<len(data)<1_000_000 and p==paths.directory/'index.html'
    assert b'<link' not in data and b'<script src' not in data
    assert status.encode() in data
    assert json.loads((paths.directory/'manifest.json').read_text())['report']['bytes']==len(data)

def test_report_compresses_large_diagnostics(tmp_path):
    from PIL import Image
    from hstpipe.uvis.reporting import generate_report
    paths=ReductionPaths(tmp_path,'T','uvis','F275W');paths.area('diagnostics',create=True)
    for i in range(6):Image.fromarray(np.random.default_rng(i).integers(0,256,(1400,1800,3),dtype='u1')).save(paths.area('diagnostics')/f'T_F275W_relative_sources_{i}.png')
    save(paths.directory/'manifest.json',dict(target='T',filters={'F275W':dict(status='failed')}))
    p=generate_report(paths.directory);assert p.stat().st_size<1_000_000 and b'data:image/jpeg;base64,' in p.read_bytes()

def test_workflow_failure_never_labels_old_drc_as_current(tmp_path):
    from hstpipe.uvis.reporting import generate_report
    paths=ReductionPaths(tmp_path,'T','uvis','F275W')
    save(paths.directory/'manifest.json',dict(target='T',workflow_error='Legacy mapping conflict',filters={'F275W':dict(status='complete',products=[{'path':'/old/T_F275W_drc.fits'}])}))
    output=generate_report(paths.directory);text=output.read_text()
    assert 'blocked_workflow' in text and 'Legacy mapping conflict' in text
    assert 'T_F275W_drc.fits' not in text and output.stat().st_size<1_000_000

def test_report_explicitly_warns_for_unreliable_background(tmp_path):
    from hstpipe.uvis.reporting import generate_report
    paths=ReductionPaths(tmp_path,'T','uvis','F275W')
    sky=dict(extver=1,reliable=False,measured=.1,applied=0.,warnings=['Insufficient uncontaminated sky'])
    save(paths.directory/'manifest.json',dict(target='T',filters={'F275W':dict(status='relative_astrometry_unvalidated',scientific_assessment=dict(background={'a_flc.fits':[sky]},background_discontinuity_range=.003))}))
    text=generate_report(paths.directory).read_text()
    assert 'Background assessment</td><td>WARNING' in text
    assert 'Insufficient uncontaminated sky' in text and '0.003' in text
    assert 'are not confirmed stars' in text

def test_invalid_target_does_not_block_another_and_raw_is_immutable(flc,tmp_path):
    raw=tmp_path/'raw';raw.mkdir()
    for name,target in [('a','A'),('b','B')]:
        p=flc(name)
        with fits.open(p,mode='update') as hd:
            hd[0].header['TARGNAME']=target
            if target=='B':hd['ERR',7].data=np.ones((2,2),'f4')
        p.rename(raw/p.name)
    before={p:identity(p) for p in raw.glob('*.fits')};calls=[]
    def engine(config,paths,rows):calls.append(config.target);return []
    cfg=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',options=settings({}),resume=True)
    run_uvis(cfg,engine,_retry_anchors=False)
    assert calls==['A'] and all(identity(p)==d for p,d in before.items())
    for target in ('A','B'):assert (tmp_path/f'redux/uvis/{target}/index.html').exists()

@pytest.mark.parametrize('chips',[(1,2),(7,)])
def test_full_field_sources_and_catalogs_on_every_sci(flc,tmp_path,chips):
    from hstpipe.uvis.astrometry.full_field import detect_full_field
    p=flc(size=192,chips=chips)
    with fits.open(p,mode='update') as hd:
        for k,v in enumerate(chips):
            hd['SCI',v].header['CCDCHIP']=2-k;hd['SCI',v].header['CRVAL2']+=k*.025
            yy,xx=np.indices(hd['SCI',v].shape)
            for x,y in ((35,35),(150,35),(35,150),(150,150)):hd['SCI',v].data+=150*np.exp(-((xx-x)**2+(yy-y)**2)/4)
    before=identity(p);paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    cat,report=detect_full_field([p],paths,settings({})['astrometry']['relative'])
    assert set(cat[0]['extver'])==set(chips)
    for c in report['exposures'][0]['chips']:
        assert c['valid_unique']>=4 and Path(c['catalog']).is_file()
        table=Table.read(c['catalog'],format='ascii.ecsv');assert set(table['ccdchip'])=={c['ccdchip']}
    assert identity(p)==before

def test_exposure_fit_uses_and_updates_both_chips(flc,tmp_path):
    from hstpipe.uvis.astrometry.relative_engine import native_graph_fit,evaluate
    from hstpipe.uvis.astrometry.centroids import validate_centroids
    files=[flc(n,size=192,chips=(1,2)) for n in ('a','b')];cats=[];options=settings({})['astrometry']['relative']
    for i,p in enumerate(files):
        parts=[]
        with fits.open(p,mode='update') as hd:
            for k,v in enumerate((1,2)):
                h=hd['SCI',v].header;h['CCDCHIP']=2-k;h['CRVAL2']+=k*.025
                h['CD1_1']=h.pop('CDELT1');h['CD2_2']=h.pop('CDELT2');h['CD1_2']=0.;h['CD2_1']=0.;h['CRPIX1']+=3*i
                yy,xx=np.indices(hd['SCI',v].shape)
                for j,(x,y) in enumerate(((35,35),(150,35),(35,150),(150,150))):
                    hd['SCI',v].data+=200*np.exp(-((xx-x)**2+(yy-y)**2)/3)
                    parts.append(dict(x=float(x),y=float(y),extver=v,common_source_id=f'{v}-{j}',snr=100.,width_pix=1.3,minor_width_pix=1.3))
        cats.append(validate_centroids(p,Table(rows=parts),options['detection']))
    native_graph_fit(files,cats,options)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('intermediate','alignment',create=True)
    result=evaluate(files,cats,options,paths,'chip_fit')
    assert result['may_combine']
    assert {r['extver'] for r in result['pairs'][0]['chip_matches'] if r['robust_inliers']}=={1,2}
    with fits.open(files[1]) as hd:assert all(hd['SCI',v].header['WCSNAME']=='HP-RELATIVE' for v in (1,2))

def test_no_restrictive_drizzle_group(monkeypatch,tmp_path):
    from hstpipe.uvis.drizzle import drizzle
    from drizzlepac import astrodrizzle
    calls=[];monkeypatch.setattr(astrodrizzle,'AstroDrizzle',lambda *a,**kw:calls.append((a,kw)))
    drizzle([tmp_path/'full_flc.fits'],tmp_path/'out')
    assert 'group' not in calls[0][1] and calls[0][0][0]==[str(tmp_path/'full_flc.fits')]
    with pytest.raises(ValueError):drizzle([tmp_path/'full_flc.fits'],tmp_path/'out',dict(group=1))

def test_standard_failure_invokes_target_fallback_without_propagation(flc,tmp_path,monkeypatch):
    from hstpipe.uvis.astrometry import relative_engine as e,target_fallback
    files=[flc(n) for n in ('a','b')];paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    paths.area('intermediate','alignment',create=True)
    options=settings({})['astrometry']['relative'];options['methods']=[]
    monkeypatch.setattr(e,'inspect_wcs',lambda inputs:[dict(chips=[dict(solutions={' ':'BASE'},primary='BASE')]) for p in inputs])
    calls=[]
    def fallback(*args,**kwargs):calls.append(True);return dict(status='unvalidated',reason='No reliable compact source'),[]
    monkeypatch.setattr(target_fallback,'target_alignment',fallback)
    def mosaic(inputs,prefix,reference=None):
        p=prefix.with_name(prefix.name+'_drc.fits');fits.PrimaryHDU().writeto(p);return p
    monkeypatch.setattr(e,'internal_mosaic',mosaic)
    before=[identity(p) for p in files]
    report=e.refine_relative(files,paths,options)
    assert calls==[True] and not report['may_combine'] and report['result_state']=='unvalidated'
    assert [identity(p) for p in files]==before

def test_native_tweakreg_uses_clean_reference_and_both_chips(flc,tmp_path,monkeypatch):
    from drizzlepac import catalogs as native_catalogs
    from hstpipe.uvis.astrometry.full_field import detect_full_field
    from hstpipe.uvis.astrometry.relative_engine import tweak,evaluate
    files=[flc(n,size=192,chips=(1,2)) for n in ('a','b')]
    for i,p in enumerate(files):
        with fits.open(p,mode='update') as hd:
            for k,v in enumerate((1,2)):
                h=hd['SCI',v].header;h['CCDCHIP']=2-k;h['CRVAL2']+=k*.025
                h['CD1_1']=h.pop('CDELT1');h['CD2_2']=h.pop('CDELT2');h['CD1_2']=0.;h['CD2_1']=0.;h['CRPIX1']+=3*i
                yy,xx=np.indices(hd['SCI',v].shape)
                for x,y in ((35,35),(150,35),(35,150),(150,150)):hd['SCI',v].data+=200*np.exp(-((xx-x)**2+(yy-y)**2)/3)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');options=settings({})['astrometry']['relative']
    cats,_=detect_full_field(files,paths,options)
    def forbidden(*a,**kw):raise AssertionError('Native automatic reference detection bypasses validated catalogs')
    monkeypatch.setattr(native_catalogs.ImageCatalog,'generateXY',forbidden)
    with paths.scratch('native-tweak-') as directory:
        tweak(files,cats,Path(directory),options,'shift',True)
    result=evaluate(files,cats,options,paths,'native_tweak')
    assert result['may_combine']
    assert {c['extver'] for c in result['pairs'][0]['chip_matches'] if c['robust_inliers']}=={1,2}
