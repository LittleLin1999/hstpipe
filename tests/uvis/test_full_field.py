"""Focused full-field catalogs, HST-only associations and complete WCS transforms."""
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS,Sip
from hstpipe.common.paths import ReductionPaths
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry.full_field import detect_full_field
from hstpipe.uvis.astrometry.common_sources import build_common


def star_field(flc,name):
    p=flc(name,size=256)
    xy=np.array([[35,35],[220,35],[35,220],[220,220]],float)
    with fits.open(p,mode='update') as hd:
        y,x=np.indices(hd['SCI'].shape)
        for xx,yy in xy:hd['SCI'].data+=100*np.exp(-((x-xx)**2+(y-yy)**2)/4)
    return p,xy


def test_full_field_keeps_sources_far_from_target_and_deduplicates(flc,tmp_path,monkeypatch):
    p,xy=star_field(flc,'one');paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    options=settings({})['astrometry']['relative'];cats,report=detect_full_field([p],paths,options)
    cat=cats[0]
    for x,y in xy:assert np.min(np.hypot(cat['x']-x,cat['y']-y))<.5
    assert len(set(cat['hst_source_id']))==len(cat)
    assert len(cat)<sum(m['accepted'] for m in report['exposures'][0]['methods'])
    from hstpipe.uvis.astrometry import full_field
    monkeypatch.setattr(full_field,'sources',lambda *a,**k:(_ for _ in ()).throw(AssertionError('Cached full-field detection must resume')))
    cached,_=detect_full_field([p],paths,options);assert len(cached[0])==len(cat)


def test_hst_only_common_sources_without_any_external_catalog(flc,tmp_path):
    files=[star_field(flc,n)[0] for n in ('a','b','c')];paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');options=settings({})['astrometry']['relative']
    cats,_=detect_full_field(files,paths,options)
    external=Table(names=['source_id','exposure','status'],dtype=['U64']*3)
    common,table,report=build_common(files,cats,external,options)
    assert min(map(len,common))>=4 and report['external_positions_used_for_fit'] is False
    assert any(str(r['source_id']).startswith('HST-') and r['status']=='common_candidate' for r in table)


def test_distortion_projection_uses_stwcs_sip_model(flc):
    from stwcs.wcsutil import HSTWCS
    from hstpipe.uvis.astrometry.sources import project
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:hd['SCI'].header['A_2_0']=2e-5;hd['SCI'].header['B_0_2']=-1e-5
    xy=np.array([[20.,20.],[220.,220.]])
    with fits.open(p) as hd:
        full=HSTWCS(hd,ext=('SCI',7));tangent=WCS(hd['SCI'].header,hd)
        expected=tangent.wcs_world2pix(full.all_pix2world(xy,0),0)
    cat=Table(dict(x=xy[:,0],y=xy[:,1],extver=[7,7]))
    actual=project(cat,p,tangent)
    assert np.allclose(actual,expected,atol=1e-6) and np.max(abs(actual-xy))>.1


def test_external_hst_identity_rejects_conflicting_histogram_peak(flc,tmp_path,monkeypatch):
    from hstpipe.uvis.astrometry import common_sources as c
    files=[star_field(flc,n)[0] for n in ('a','b')];paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');options=settings({})['astrometry']['relative']
    cats,_=detect_full_field(files,paths,options)
    external=Table(dict(source_id=['L1','L1'],exposure=[p.name for p in files],status=['valid_detection']*2,hst_source_id=[str(x['hst_source_id'][0]) for x in cats]))
    monkeypatch.setattr(c,'match',lambda *a,**k:(np.arange(4),np.arange(4)))
    monkeypatch.setattr(c,'quality',lambda *a,**k:dict(status='refinement_supported',shift=[15,0],fit_inlier_mask=[True]*4))
    _,_,r=c.build_common(files,cats,external,options)
    assert r['pairs'][0]['status']=='contradicts_common_HST_identifications'


def test_single_source_local_fallback_never_claims_full_field_validation(flc,tmp_path):
    from hstpipe.uvis.astrometry.relative_engine import single_source_fallback
    files=[flc(n,size=256) for n in ('a','b')]
    with fits.open(files[1],mode='update') as hd:hd['SCI'].header['CRPIX1']+=5
    cats=[Table(dict(x=[128.],y=[128.],extver=[7],centroid_error_pix=[.1],common_source_id=['target'])) for p in files]
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('intermediate','alignment',create=True)
    r=single_source_fallback(files,cats,settings({})['astrometry']['relative'],paths)
    assert r['status']=='local_only' and not r['full_field_validated'] and r['independent_scatter_pix'] is None
    assert np.allclose(r['transformations'][0]['shift'],[5,0],atol=.001)


def test_identical_failed_requests_are_not_reissued(monkeypatch):
    from hstpipe.uvis.astrometry import catalogs as c
    c._FAILED_REQUESTS.clear();calls=[]
    def fail(*a):calls.append(1);raise RuntimeError('offline')
    monkeypatch.setattr(c,'_request',fail)
    import pytest
    with pytest.raises(c.CatalogQueryError):c.request('legacy_dr10','unique_failure',timeout=None,retries=0)
    with pytest.raises(c.CatalogQueryError,match='already failed'):c.request('legacy_dr10','unique_failure',timeout=None,retries=0)
    assert calls==[1];c._FAILED_REQUESTS.clear()


def test_alternate_gaia_endpoint_uses_documented_client(monkeypatch):
    import sys,types
    from hstpipe.uvis.astrometry import catalogs as c
    seen={}
    class Client:
        def __init__(self,**kw):seen.update(kw)
        def launch_job_async(self,q,**kw):
            seen['query']=q;seen['output_format']=kw.get('output_format');return types.SimpleNamespace(get_results=lambda:Table(dict(ra=[10.],dec=[20.])))
    monkeypatch.setitem(sys.modules,'astroquery.gaia',types.SimpleNamespace(GaiaClass=Client))
    table=c._request('gaia_dr3','SELECT ra, dec FROM gaiadr3.gaia_source','https://gaia.ari.uni-heidelberg.de/tap')
    assert len(table)==1 and seen['tap_server_context']=='' and not seen['show_server_messages'] and seen['output_format']=='votable'


def test_real_gaia_masked_integer_response_normalizes_at_hst_epoch():
    from astropy.table import MaskedColumn
    from hstpipe.uvis.astrometry.catalogs import standardize_gaia
    t=Table(dict(source_id=[3887236262361041792],ra=[156.39199],dec=[14.00973],ra_error=[.1],dec_error=[.1],pmra=[-5.8],pmdec=[-.25],pmra_error=[.1],pmdec_error=[.1],parallax=[-1.],parallax_error=[1.],ref_epoch=[2016.],phot_g_mean_mag=[16.],ruwe=[1.],duplicated_source=[False]))
    t['astrometric_params_solved']=MaskedColumn([31],dtype='i2',mask=[False])
    t['visibility_periods_used']=MaskedColumn([12],dtype='i2',mask=[False])
    normalized=standardize_gaia(t,2026.08)
    assert len(normalized)==1 and normalized['source_id'][0]=='3887236262361041792'
    assert normalized['pm_applied'][0] and np.isfinite(normalized['ra'][0])


def test_native_graph_fits_chain_without_direct_reference_matches(flc,tmp_path):
    from hstpipe.uvis.astrometry import relative_engine as e
    from hstpipe.uvis.astrometry.centroids import validate_centroids
    from hstpipe.uvis.astrometry.wcs import detector_wcs
    files=[flc(n,size=256) for n in ('a','b','c')]
    a=np.array([[20,20],[100,20],[20,100],[100,100]],float);b=a+130
    catalogs=[]
    for index,p in enumerate(files):
        xy=a if index==0 else (np.vstack([a,b]) if index==1 else b)
        ids=['A'+str(k) for k in range(4)] if index==0 else (['A'+str(k) for k in range(4)]+['B'+str(k) for k in range(4)] if index==1 else ['B'+str(k) for k in range(4)])
        with fits.open(p,mode='update') as hd:
            h=hd['SCI'].header;h['CD1_1']=h.pop('CDELT1');h['CD2_2']=h.pop('CDELT2');h['CD1_2']=0.;h['CD2_1']=0.;h['CRPIX1']+=[0,3,7][index]
            yy,xx=np.indices(hd['SCI'].shape)
            for x,y in xy:hd['SCI'].data+=200*np.exp(-((xx-x)**2+(yy-y)**2)/3)
        raw=Table(dict(x=xy[:,0],y=xy[:,1],extver=np.full(len(xy),7),common_source_id=ids,snr=np.full(len(xy),100.),width_pix=np.full(len(xy),1.3),minor_width_pix=np.full(len(xy),1.3)))
        catalogs.append(validate_centroids(p,raw,settings({})['astrometry']['relative']['detection']))
    options=settings({})['astrometry']['relative']
    records=e.native_graph_fit(files,catalogs,options)
    assert len(records)==2 and records[1]['parent_reference']==files[1].name
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('intermediate','alignment',create=True);paths.area('diagnostics','alignment',create=True)
    measured=e.evaluate(files,catalogs,options,paths,'chain')
    assert measured['may_combine'] and measured['network']['status']=='connected'
    assert measured['pairs'][1]['status']=='insufficient_matches'


def test_image_qa_allows_a_real_source_to_become_centered_after_alignment(tmp_path):
    from hstpipe.uvis.astrometry.relative_engine import image_quality
    w=WCS(naxis=2);w.wcs.crpix=[128,128];w.wcs.crval=[10,20];w.wcs.cdelt=[-.04/3600,.04/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    xy=np.array([[50,50],[200,50],[50,200],[200,200]]);yy,xx=np.indices((256,256));files=[]
    for label in ('before','after'):
        data=np.zeros((256,256))
        for x,y in xy:
            if label=='before':data+=50*np.exp(-((xx-x)**2+(yy-y)**2)/4)+50*np.exp(-((xx-x-4)**2+(yy-y)**2)/4)
            else:data+=100*np.exp(-((xx-x)**2+(yy-y)**2)/4)
        p=tmp_path/(label+'.fits');fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(data,header=w.to_header(),name='SCI'),fits.ImageHDU(np.ones(data.shape),name='WHT')]).writeto(p);files.append(p)
    sky=w.all_pix2world(xy,0);cat=Table(dict(ra=sky[:,0],dec=sky[:,1]))
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('diagnostics','alignment',create=True)
    r=image_quality(*files,cat,settings({})['astrometry']['relative'],paths,'centering')
    assert r['status']=='validated' and r['after']['fwhm_pix']<r['before']['fwhm_pix']


def test_cached_detector_models_preserve_full_distortion_projection(flc):
    from hstpipe.uvis.astrometry.identification import projection_models,projections
    p=flc(size=256,ra=359.999)
    with fits.open(p,mode='update') as hd:hd['SCI'].header['A_2_0']=2e-5
    options=settings({})['astrometry']['source_identification']
    models=projection_models(p)
    for ra in (359.9989,359.999,359.9991):
        uncached=projections(p,ra,20.,options);cached=projections(p,ra,20.,options,models=models)
        assert uncached==cached
