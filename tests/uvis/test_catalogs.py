import numpy as np
import pytest
import astropy.units as u
from astropy.table import Table,MaskedColumn
from astropy.io import fits
from astropy.wcs import WCS
from hstpipe.uvis.astrometry import catalogs as c


def gaia(n=8):
    t=Table(dict(source_id=np.arange(n,dtype='i8')+123456789123456789,ra=10+np.arange(n)*.001,dec=np.full(n,20.),
        ra_error=np.ones(n),dec_error=np.ones(n),pmra=np.full(n,50.),pmdec=np.full(n,-10.),pmra_error=np.ones(n),pmdec_error=np.ones(n),
        parallax=np.full(n,-1.),parallax_error=np.ones(n),ref_epoch=np.full(n,2016.),phot_g_mean_mag=np.full(n,18.),
        ruwe=np.ones(n),astrometric_params_solved=np.full(n,31),visibility_periods_used=np.full(n,12),duplicated_source=np.zeros(n,bool)))
    return t


def legacy(n=8):
    return Table(dict(ls_id=np.arange(n),ra=10+np.arange(n)*.001,dec=np.full(n,20.),type=np.full(n,'PSF'),brick_primary=np.ones(n,bool),
        maskbits=np.zeros(n,'i4'),flux_r=np.full(n,100.),flux_z=np.full(n,100.),flux_ivar_r=np.ones(n),flux_ivar_z=np.ones(n),
        photsys=np.array(['N','S']*(n//2)+(['N'] if n%2 else [])),release=np.full(n,9999),nobs_r=np.ones(n,'i4'),nobs_z=np.ones(n,'i4')))


def test_gaia_adql_and_numeric_safety():
    q=c.gaia_adql(360.01,-20,.01)
    assert 'gaiadr3.gaia_source' in q and 'CONTAINS' in q and 'pmra_error' in q and '0.010000000000' in q
    for args in [(float('nan'),20,.01),(10,91,.01),(10,20,-1),('10;DROP',20,.1)]:
        with pytest.raises(ValueError):c.gaia_adql(*args)


def test_gaia_motion_and_uncertainty():
    t=gaia();out=c.standardize_gaia(t,2026.)
    assert len(out)==8 and out['ra'].unit==u.deg and out['ra_error_mas'].unit==u.mas
    assert np.allclose((out['ra']-t['ra'])*3600000*np.cos(np.deg2rad(20)),500,atol=.1)
    assert np.allclose((out['dec']-t['dec'])*3600000,-100,atol=.1)
    assert np.all(out['pm_applied']) and np.all(out['position_epoch']==2026.)
    assert np.all(out['ra_error_mas']>10) and out['source_id'][0]==str(t['source_id'][0])


def test_gaia_masked_pm_and_invalid_parallax():
    t=gaia();t['pmra']=MaskedColumn(t['pmra'],mask=[True]+[False]*7)
    t['parallax']=[-5,0,np.nan,10,10,1,1,1];t['parallax_error']=[1,1,1,1,100,1,1,1]
    out=c.standardize_gaia(t,2020.)
    assert len(out)==7 and np.isfinite(out['ra']).all()
    assert 'space_motion_distance' in out['propagation_method'][2]
    allowed=c.standardize_gaia(t,2020.,dict(allow_missing_pm=True))
    assert not allowed['pm_applied'][0] and allowed['position_epoch'][0]==2016.


def test_gaia_quality_duplicate_and_empty():
    t=gaia();t['ruwe'][0]=3;t['source_id'][2]=t['source_id'][1]
    assert len(c.standardize_gaia(t))==6
    empty=c.standardize_gaia(t[:0],2020.);assert len(empty)==0 and empty['ra'].unit==u.deg
    with pytest.raises(c.CatalogSchemaError):c.standardize_gaia(Table(dict(ra=[np.nan],dec=[0.])))


def test_legacy_schema_query_and_releases(monkeypatch):
    calls=[]
    def request(kind,q,**kw):
        calls.append((kind,q,kw))
        if kind=='legacy_schema':return Table(dict(column_name=legacy().colnames,datatype=['double']*len(legacy().colnames)))
        return legacy()
    monkeypatch.setattr(c,'request',request)
    out=c.query_legacy_dr10(10,20,.01)
    assert 'Q3C_RADIAL_QUERY' in calls[1][1] and 'ls_dr10.tractor' in calls[1][1]
    assert calls[1][2]['tap_url']=='https://datalab.noirlab.edu/tap'
    assert out['release'][0]=='DR9-North' and out['release'][1]=='DR10-South'
    assert out['ra'].unit==u.deg and len(out)==8
    assert 'TAP_SCHEMA.columns' in calls[0][1]


def test_legacy_selection_and_types():
    t=legacy();t['type'][0]='EXP';t['maskbits'][1]=2;t['nobs_r'][2]=0;t['nobs_z'][2]=0
    t['flux_r'][3]=0;t['flux_z'][3]=0
    assert len(c.standardize_legacy(t))==4
    t['brick_primary']=np.array(['f']+['t']*7)
    assert len(c.standardize_legacy(t))==4
    assert len(c.standardize_legacy(t[:0]))==0
    del t['photsys'];out=c.standardize_legacy(t)
    assert set(out['release'])=={'DR10-family-unresolved'}


def test_missing_or_malformed_responses():
    t=legacy();del t['ra']
    with pytest.raises(c.CatalogSchemaError):c.standardize_legacy(t)
    t=legacy();t['ra']=['x']*len(t)
    with pytest.raises(c.CatalogSchemaError):c.standardize_legacy(t)
    t=legacy();t['ra']=np.full(len(t),999.)
    with pytest.raises(c.CatalogSchemaError):c.standardize_legacy(t)


def test_cache_offline_refresh_and_epoch_invalidation(monkeypatch,tmp_path):
    calls=[]
    def request(*a,**kw):calls.append(a);return gaia()
    monkeypatch.setattr(c,'request',request)
    first=c.query_gaia_dr3(10,20,.01,2020.,cache_dir=tmp_path)
    c.query_gaia_dr3(10,20,.01,2020.,cache_dir=tmp_path,offline=True)
    assert len(calls)==1 and first.meta['cache_file']
    c.query_gaia_dr3(10,20,.01,2021.,cache_dir=tmp_path,offline=True);assert len(calls)==1
    c.query_gaia_dr3(10,20,.01,2021.,cache_dir=tmp_path,refresh=True);assert len(calls)==2
    with pytest.raises(c.CatalogOfflineError):c.query_gaia_dr3(10,20,.02,2021.,cache_dir=tmp_path,offline=True)
    assert len(calls)==2


def test_transient_retry_and_exhaustion(monkeypatch):
    attempts=[]
    def flaky(*a):
        attempts.append(a)
        if len(attempts)<3:raise OSError('temporary service failure')
        return gaia()
    monkeypatch.setattr(c,'_request',flaky)
    assert len(c.request('gaia_dr3','query',timeout=None,retries=2,backoff=0))==8 and len(attempts)==3
    monkeypatch.setattr(c,'_request',lambda *a:(_ for _ in ()).throw(OSError('down')))
    with pytest.raises(c.CatalogQueryError,match='2 attempts'):c.request('gaia_dr3','q',timeout=None,retries=1,backoff=0)


@pytest.mark.parametrize('ra,dec',[(359.999,20),(0.001,-20),(12,89.8),(250,-89.8)])
def test_actual_footprint_wrap_and_poles(tmp_path,ra,dec):
    w=WCS(naxis=2);w.wcs.crpix=[150,50];w.wcs.cdelt=[-.04/3600,.04/3600];w.wcs.crval=[ra,dec];w.wcs.ctype=['RA---TAN','DEC--TAN']
    p=tmp_path/'image.fits';fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(np.ones((100,300)),header=w.to_header(),name='SCI'),
        fits.ImageHDU(np.ones((100,300)),name='WHT')]).writeto(p)
    cone=c.query_footprint(p,margin_arcsec=5)
    assert 0<cone['radius']<.01 and np.isfinite(list(cone.values())).all()


def test_epoch_re_evaluation_uses_native_gaia_positions():
    first=c.standardize_gaia(gaia(),2020.)
    again=c.gaia_at_epoch(first,2026.)
    direct=c.standardize_gaia(gaia(),2026.)
    assert np.allclose(again['ra'],direct['ra'],atol=1e-11,rtol=0)


def test_release_codes_without_photsys():
    t=legacy();del t['photsys'];t['release']=[9011,10000,10002,9999,9011,10000,10002,9999]
    out=c.standardize_legacy(t)
    assert list(out['release'][:4])==['DR9-North','DR10-South','DR10-South','DR10-family-unresolved']


def test_documented_client_dispatch(monkeypatch):
    import sys,types
    calls=[]
    class Job:
        def get_results(self):return gaia()
    class GaiaClient:
        def __init__(self,**kw):calls.append(('gaia_client',kw))
        def launch_job_async(self,q,**kw):calls.append(('gaia',q,kw));return Job()
    monkeypatch.setitem(sys.modules,'astroquery.gaia',types.SimpleNamespace(GaiaClass=GaiaClient))
    assert len(c._request('gaia_dr3','SELECT test',c.TAP_URL))==8
    class Service:
        def __init__(self,url):calls.append(('url',url))
        def run_sync(self,q):calls.append(('pyvo',q));return types.SimpleNamespace(to_table=lambda:legacy())
    monkeypatch.setitem(sys.modules,'pyvo',types.SimpleNamespace(dal=types.SimpleNamespace(TAPService=Service)))
    assert len(c._request('legacy_dr10','SELECT tractor',c.TAP_URL))==8
    assert calls[0]==('gaia_client',dict(gaia_tap_server='https://gea.esac.esa.int/',tap_server_context='tap-server',show_server_messages=False))
    assert calls[1][2]==dict(output_format='votable',dump_to_file=False) and calls[2]==('url',c.TAP_URL)


def test_malformed_refresh_preserves_previous_cache(monkeypatch,tmp_path):
    from hstpipe.common.provenance import identity
    monkeypatch.setattr(c,'request',lambda *a,**k:gaia())
    first=c.query_gaia_dr3(10,20,.01,2020.,cache_dir=tmp_path)
    path=next(tmp_path.glob('gaia_dr3_raw_*.ecsv'));before=identity(path)
    bad=gaia();bad['ra']=np.full(len(bad),200.)
    monkeypatch.setattr(c,'request',lambda *a,**k:bad)
    with pytest.raises(c.CatalogSchemaError):c.query_gaia_dr3(10,20,.01,2020.,cache_dir=tmp_path,refresh=True)
    assert identity(path)==before
    assert len(c.query_gaia_dr3(10,20,.01,2020.,cache_dir=tmp_path,offline=True))==len(first)


def test_masked_coordinates_and_quality_flags_never_become_valid():
    t=gaia();t['ra']=MaskedColumn(t['ra'],mask=[True]+[False]*7,unit=u.deg)
    assert len(c.standardize_gaia(t))==7
    t=legacy();t['maskbits']=MaskedColumn(t['maskbits'],mask=[True]+[False]*7)
    assert len(c.standardize_legacy(t))==7


def test_angular_units_are_preserved_in_native_gaia_coordinates():
    t=gaia();t['ra']=np.deg2rad(t['ra'])*u.rad;t['dec']=np.deg2rad(t['dec'])*u.rad
    out=c.standardize_gaia(t,2020.)
    assert np.allclose(out['ra_reference'],gaia()['ra'])
    assert np.allclose(c.gaia_at_epoch(out,2026.)['ra'],c.standardize_gaia(gaia(),2026.)['ra'])


def test_no_usable_legacy_coverage_is_reported():
    t=legacy();t['nobs_r'][:]=0;t['nobs_z'][:]=0
    out=c.standardize_legacy(t)
    assert not len(out) and out.meta['coverage_source_count']==0


def test_timeout_terminates_client_worker_without_network(monkeypatch):
    from types import SimpleNamespace
    calls=[]
    read=SimpleNamespace(poll=lambda timeout:False,close=lambda:calls.append('close_read'))
    write=SimpleNamespace(close=lambda:calls.append('close_write'))
    process=SimpleNamespace(start=lambda:calls.append('start'),is_alive=lambda:True,
        terminate=lambda:calls.append('terminate'),join=lambda:calls.append('join'))
    context=SimpleNamespace(Pipe=lambda **kw:(read,write),Process=lambda **kw:process)
    monkeypatch.setattr(c.multiprocessing,'get_context',lambda method:context)
    with pytest.raises(c.CatalogQueryError,match='exceeded'):c.request('gaia_dr3','q',timeout=.01,retries=0)
    assert 'terminate' in calls and 'join' in calls
