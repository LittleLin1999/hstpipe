import json,shutil,types
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry.matching import science_header
from hstpipe.uvis.inventory import discover
from hstpipe.uvis.astrometry.matching import validate
from hstpipe.uvis.astrometry.catalogs import CatalogQueryError


@pytest.fixture(autouse=True)
def no_unmocked_catalog_queries(monkeypatch):
    from hstpipe.uvis.astrometry import alignment
    def forbidden(*a,**k):raise AssertionError('Unmocked catalog provider called in a unit test')
    monkeypatch.setattr(alignment,'query_gaia_dr3',forbidden)
    monkeypatch.setattr(alignment,'query_legacy_dr10',forbidden)


def grid(image):
    with fits.open(image) as hd:
        w=WCS(hd['SCI'].header,hd);xy=np.array([[20,20],[100,20],[230,20],[20,100],[230,100],[20,230],[100,230],[230,230]])
        world=w.all_pix2world(xy,0)
    return Table(dict(x=xy[:,0].astype(float),y=xy[:,1].astype(float),ra=world[:,0],dec=world[:,1],flux=np.full(8,100.)))


def preparation(flc,tmp_path):
    p=flc(size=256);rows=discover(tmp_path)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    for a,s in [('intermediate','alignment'),('intermediate','drizzle'),('diagnostics','alignment'),('diagnostics','drizzle'),('logs',None)]:paths.area(a,s,create=True)
    ref=grid(p);ref['source_id']=[str(i) for i in range(8)]
    catalog=paths.catalog_directory(create=True)/'fixture.ecsv';ref.write(catalog,format='ascii.ecsv');ref.meta['cache_file']=str(catalog)
    return p,rows,paths,ref,settings({})['astrometry']


def fake_tweakreg(files,**kw):
    for p in files:
        with fits.open(p,mode='update') as hd:
            hd['SCI'].header['WCSNAME']=kw['wcsname']
            hd['SCI'].header['CRDER1']=.1/3600
            hd['SCI'].header['CRDER2']=.15/3600
    Path(kw['outshifts']).write_text(f'{files[0]} 0 0 0 1 0 0\n')
    Path('fixture_catalog_fit.match').write_text('1 2 3 4\n'*8)


def test_existing_gaia_is_validated_without_new_fit(flc,tmp_path,monkeypatch):
    pytest.importorskip('stwcs');import drizzlepac.tweakreg
    from hstpipe.uvis.astrometry import alignment as a
    p,rows,paths,ref,opts=preparation(flc,tmp_path)
    with fits.open(p,mode='update') as hd:hd['SCI'].header['WCSNAME']='IDC-GAIA-DR3'
    monkeypatch.setattr(a,'detect_sources',lambda *args,**kwargs:grid(args[0]));monkeypatch.setattr(a,'query_gaia_dr3',lambda *args,**kw:ref)
    monkeypatch.setattr(drizzlepac.tweakreg,'TweakReg',lambda *a,**k:pytest.fail('Validated existing Gaia WCS must not be refitted'))
    before_signature=a.signature(p)
    result,_=a.absolute_alignment(p,rows,paths,opts)
    assert a.signature(p)==before_signature
    result,_=a.absolute_alignment(p,rows,paths,opts)
    assert result['status']=='validated' and not result['correction_applied'] and result['validation']['matched']==8, result
    assert science_header(p).get('HPORIGK')


@pytest.mark.parametrize('query_failure',[True,False])
def test_legacy_fallback_after_query_or_science_failure(flc,tmp_path,monkeypatch,query_failure):
    pytest.importorskip('stwcs');import drizzlepac.tweakreg
    from hstpipe.uvis.astrometry import alignment as a
    p,rows,paths,ref,opts=preparation(flc,tmp_path)
    def gaia(*args,**kw):
        if query_failure:raise CatalogQueryError('service down')
        return ref[:1]
    monkeypatch.setattr(a,'query_gaia_dr3',gaia);monkeypatch.setattr(a,'query_legacy_dr10',lambda *a,**k:ref)
    monkeypatch.setattr(a,'legacy_at_epoch',lambda r,e:r);monkeypatch.setattr(a,'detect_sources',lambda p,*a,**kwargs:grid(p))
    monkeypatch.setattr(drizzlepac.tweakreg,'TweakReg',fake_tweakreg)
    result,_=a.absolute_alignment(p,rows,paths,opts)
    assert result['status']=='validated' and result['reference']=='legacy_dr10', result
    assert result['attempts'][0]['status']==('external_catalog_unavailable' if query_failure else 'diagnostic_single_reference')
    assert result['correction_applied'] and result['transformations']
    assert science_header(p)['WCSNAME']=='HP-LEGACY'


def test_sparse_field_does_not_claim_absolute_calibration(flc,tmp_path,monkeypatch):
    pytest.importorskip('stwcs')
    from hstpipe.uvis.astrometry import alignment as a
    p,rows,paths,ref,opts=preparation(flc,tmp_path)
    monkeypatch.setattr(a,'query_gaia_dr3',lambda *a,**k:ref[:1]);monkeypatch.setattr(a,'query_legacy_dr10',lambda *a,**k:ref[:1])
    monkeypatch.setattr(a,'legacy_at_epoch',lambda r,e:r)
    result,refs=a.absolute_alignment(p,rows,paths,opts)
    assert result['status']=='unvalidated' and refs is None
    assert all(r['status'] in ('insufficient_verified_references','diagnostic_single_reference') for r in result['attempts'])


def test_matching_outliers_and_distribution():
    src=Table(dict(x=np.arange(8)*20.,y=np.arange(8)*20.,ra=10+np.arange(8)*.001,dec=np.full(8,20.)))
    ref=src.copy();ref['ra'][-1]+=.0002
    opts=dict(min_matches=6,search_radius_arcsec=2,max_rms_arcsec=.1,max_offset_arcsec=.1,systematic_floor_arcsec=.01,min_span_fraction=.1)
    result,matched=validate(src,ref,opts,(256,256))
    assert result['matched']==7 and result['rejected']==1 and result['status']=='validated'
    src['x']=np.ones(8);src['y']=np.ones(8)
    result,_=validate(src,ref,opts,(256,256));assert result['status']!='validated'


def test_cross_filter_uses_validated_anchor(flc,tmp_path,monkeypatch):
    pytest.importorskip('stwcs');import drizzlepac.tweakreg
    from hstpipe.uvis.astrometry import alignment as a
    p,rows,paths,ref,opts=preparation(flc,tmp_path)
    anchor=tmp_path/'anchor.fits';shutil.copy2(p,anchor)
    state=dict(filters={'F814W':dict(filter='F814W',products=[identity(anchor)],scientific_assessment={'relative_astrometry':{'status':'validated'},'astrometry':dict(status='validated',observation_epoch=2017.7,reference='gaia_dr3',final_validation={'rms_arcsec':.02})})})
    state['filters']['F814W']['scientific_assessment']['cross_filter_alignment']=dict(status='anchor',anchor_group='F814W')
    state['filters'][paths.group_key]=dict(scientific_assessment=dict(cross_filter_alignment=dict(status='validated',anchor_group='F814W')))
    (paths.directory/'manifest.json').write_text(json.dumps(state))
    (paths.file('intermediate','alignment','cross_filter')).write_text(json.dumps(dict(status='validated',anchor_group='F814W')))
    from hstpipe.uvis.astrometry import sources as source_module
    monkeypatch.setattr(source_module,'sources',lambda p,*args:grid(p))
    monkeypatch.setattr(a,'detect_sources',lambda p,*a,**kwargs:grid(p));monkeypatch.setattr(drizzlepac.tweakreg,'TweakReg',fake_tweakreg)
    monkeypatch.setattr(a,'query_gaia_dr3',lambda *a,**k:pytest.fail('Cross-filter anchor should avoid a separate external query'))
    result,_=a.absolute_alignment(p,rows,paths,opts)
    assert result['reference_mode']=='cross_filter' and result['cross_filter_anchor']==str(anchor), result


@pytest.mark.parametrize('chips',[(1,7),(7,)])
def test_tweakback_dispatch_extensions_and_safe_resume(flc,tmp_path,monkeypatch,chips):
    pytest.importorskip('stwcs');import drizzlepac.tweakback
    from stwcs.wcsutil import headerlet
    from hstpipe.uvis.astrometry.workflow import tweakback_inputs
    source=flc(chips=chips);before=identity(source);work=tmp_path/'work_flc.fits';shutil.copy2(source,work)
    provisional=tmp_path/'provisional.fits';shutil.copy2(source,provisional)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('intermediate','alignment',create=True)
    calls=[]
    def apply(image,**kw):
        assert image==str(provisional)+f'[SCI,{chips[0]}]'
        calls.append(kw)
        for p in kw['input_files'].split(','):
            with fits.open(p,mode='update') as hd:
                for h in hd:
                    if h.name=='SCI':h.header['WCSNAME']=kw['output_wcs_name']
    monkeypatch.setattr(drizzlepac.tweakback,'apply_tweak',apply,raising=False)
    class Hlet:
        def tofile(self,p,**k):fits.PrimaryHDU().writeto(p,overwrite=True)
    monkeypatch.setattr(headerlet,'create_headerlet',lambda *a,**k:Hlet())
    restored=[];monkeypatch.setattr(headerlet,'apply_headerlet_as_primary',lambda *a,**k:restored.append(a))
    solution=dict(correction_applied=True,original_wcs_name='ORIG',original_wcs_key='A',corrected_wcs_name='HP-GAIA-DR3')
    report=tweakback_inputs(provisional,[source],[work],solution,paths,dict(enabled=True))
    assert report['updated_extensions'][0]['extensions']==list(chips) and calls[0]['default_extname']=='SCI'
    assert identity(source)==before
    shutil.copy2(source,work);tweakback_inputs(provisional,[source],[work],solution,paths,dict(enabled=True))
    assert len(calls)==1 and len(restored)==1


def test_final_drizzle_receives_corrected_flcs(flc,tmp_path,monkeypatch):
    from hstpipe.uvis.astrometry import workflow as w
    p,rows,paths,ref,opts=preparation(flc,tmp_path);before=identity(p)
    provisional=paths.file('intermediate','drizzle','provisional','fits');shutil.copy2(p,provisional)
    solution=dict(status='validated',reference='gaia_dr3',reference_mode='gaia_dr3',observation_epoch=2017.7,exposure_epoch_span_years=0,
        correction_applied=True,corrected_wcs_name='HP-GAIA-DR3',thresholds=dict(opts['absolute'],**opts['catalogs']['gaia_dr3']['validation']))
    monkeypatch.setattr(w,'absolute_alignment',lambda *a:(dict(solution),ref));monkeypatch.setattr(w,'gaia_at_epoch',lambda r,e:r)
    def tb(provisional,relative,working,solution,paths,options):
        for file in working:
            with fits.open(file,mode='update') as hd:hd['SCI'].header['WCSNAME']='HP-GAIA-DR3'
        h=paths.file('intermediate','alignment','fixture_headerlet','fits');fits.PrimaryHDU().writeto(h,overwrite=True)
        return dict(headerlets=[identity(h)],status='applied')
    monkeypatch.setattr(w,'tweakback_inputs',tb)
    called=[]
    def drizzle(inputs,prefix,options):
        with fits.open(inputs[0]) as hd:
            assert hd['SCI'].header['WCSNAME']=='HP-GAIA-DR3';called.append(inputs)
            header=hd['SCI'].header.copy();header['BUNIT']='ELECTRONS/S';shape=hd['SCI'].shape
            fits.HDUList([fits.PrimaryHDU(header=hd[0].header.copy()),fits.ImageHDU(hd['SCI'].data/100,header=header,name='SCI'),fits.ImageHDU(np.ones(shape),header=header,name='WHT'),fits.ImageHDU(np.ones(shape,'i4'),header=header,name='CTX')]).writeto(str(prefix)+'_drc.fits')
    monkeypatch.setattr(w,'drizzle',drizzle);monkeypatch.setattr(w,'detect_sources',lambda image,*a,**kw:grid(image))
    products,report=w.finalize(provisional,[p],rows,paths,opts,settings({})['drizzle'],dict(status='single_exposure'))
    assert called and products==[paths.product('drc')] and report['astrometry']['final_validation']['status']=='validated'
    with fits.open(products[0]) as hd:assert hd[0].header['ASTSTAT']=='VALIDATED' and hd['SCI'].header['BUNIT']=='ELECTRONS/S'
    assert identity(p)==before
    with fits.open(products[0]) as final,fits.open(provisional) as initial:assert not np.array_equal(final['SCI'].data,initial['SCI'].data)


def test_new_gaia_fit_cache_and_reference_invalidation(flc,tmp_path,monkeypatch):
    pytest.importorskip('stwcs');import drizzlepac.tweakreg
    from hstpipe.uvis.astrometry import alignment as a
    p,rows,paths,ref,opts=preparation(flc,tmp_path)
    monkeypatch.setattr(a,'detect_sources',lambda *args,**kwargs:grid(args[0]))
    monkeypatch.setattr(a,'query_gaia_dr3',lambda *args,**kw:ref)
    calls=[]
    def fit(files,**kw):calls.append(files);fake_tweakreg(files,**kw)
    monkeypatch.setattr(drizzlepac.tweakreg,'TweakReg',fit)
    before=a.signature(p)
    solution,_=a.absolute_alignment(p,rows,paths,opts)
    assert solution['reference']=='gaia_dr3' and solution['correction_applied']
    assert solution['validation']['matched']==8 and solution['corrected_wcs_name']=='HP-GAIA-DR3'
    assert a.signature(p)==before
    a.absolute_alignment(p,rows,paths,opts);assert len(calls)==1
    ref.meta['revision']=2;ref.write(ref.meta['cache_file'],format='ascii.ecsv',overwrite=True)
    a.absolute_alignment(p,rows,paths,opts)
    assert len(calls)==2 and a.signature(p)==before


def test_observation_epoch_can_use_fits_date_without_expstart():
    from astropy.time import Time
    from hstpipe.uvis.astrometry.alignment import observation_epoch
    row=types.SimpleNamespace(metadata={'DATE-OBS':'2017-09-04','TIME-OBS':'12:00:00'})
    epoch,span=observation_epoch([row])
    assert abs(epoch-Time('2017-09-04T12:00:00',scale='utc').tcb.jyear)<1e-10 and span==0


def test_wcs_fit_error_metadata_do_not_change_transform_cache(flc,tmp_path):
    from hstpipe.uvis.astrometry.alignment import signature
    p=flc();original=signature(p)
    with fits.open(p,mode='update') as hd:
        hd['SCI'].header['CRDER1']=1e-5;hd['SCI'].header['CSYER1']=2e-5
    assert signature(p)==original
    with fits.open(p,mode='update') as hd:hd['SCI'].header['CRVAL1']+=.001
    assert signature(p)!=original


def test_absolute_stage_reuses_complete_identification_catalog(flc,tmp_path,monkeypatch):
    from hstpipe.uvis.astrometry import alignment as a
    p,rows,paths,ref,opts=preparation(flc,tmp_path)
    ref.meta['complete_query_coverage']=True
    with fits.open(p,mode='update') as hd:hd['SCI'].header['WCSNAME']='IDC-GAIA-DR3'
    monkeypatch.setattr(a,'cached_footprint',lambda *args,**kw:ref)
    monkeypatch.setattr(a,'query_gaia_dr3',lambda *args,**kw:pytest.fail('A covering cached identification catalog must be reused'))
    monkeypatch.setattr(a,'detect_sources',lambda image,*args,**kw:grid(image))
    result,_=a.absolute_alignment(p,rows,paths,opts)
    assert result['status']=='validated' and not result['correction_applied']
