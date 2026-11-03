"""Focused external-identification / HST-only relative-fitting tests."""
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from hstpipe.common.paths import ReductionPaths
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry import identification as i
from hstpipe.uvis.astrometry.relative_fit import catalog_match

@pytest.fixture
def options():return settings({})['astrometry']

def candidate(ra,dec):
    return dict(source_id='G1',gaia_id='1',legacy_id='',ra=ra,dec=dec,classification='gaia_stellar_candidate',release='Gaia DR3')


def test_identification_policy_is_independent_of_absolute_cuts(options):
    from copy import deepcopy
    prior=deepcopy(options['catalogs']['legacy_dr10']['quality'])
    options['source_identification']['catalog_quality']={'legacy_dr10':prior}
    options['catalogs']['legacy_dr10']['quality'].update(min_snr=0.,max_position_error_mas=200.)
    selected=i.identification_quality(options,'legacy_dr10')
    assert selected['min_snr']==10 and 'max_position_error_mas' not in selected
    assert options['catalogs']['legacy_dr10']['quality']['min_snr']==0.
    assert prior==options['source_identification']['catalog_quality']['legacy_dr10']


def test_identification_default_preserves_shared_policy(options):
    options['catalogs']['legacy_dr10']['quality']['min_snr']=12.
    assert i.identification_quality(options,'legacy_dr10')['min_snr']==12.
    with pytest.raises(ValueError,match='catalog_quality'):
        settings({'astrometry':{'source_identification':{'catalog_quality':{'unknown':{}}}}})

@pytest.mark.parametrize('chips',[(7,),(1,7)])
def test_footprint_and_projection_use_actual_chips(flc,chips,options):
    p=flc(chips=chips,ra=359.999,size=256)
    cone=i.footprint([p],10);assert 0<=cone['ra']<360 and cone['radius']<.02 and len(cone['geometry'])==len(chips)
    projected=i.projections(p,359.999,20,options['source_identification'])
    assert {r['extver'] for r in projected}==set(chips) and all(r['covered'] for r in projected)


def test_bounded_local_search_measures_hst_not_predicted_centroid(flc,options):
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:
        world=WCS(hd['SCI'].header,hd).all_pix2world([[128,128]],0)[0]
        hd['SCI'].header['CRPIX1']+=5
    result=i.associate(p,candidate(*world),options['source_identification'],options['relative']['detection'])
    assert result['status']=='valid_detection',result
    assert abs(result['x']-128)<.3 and abs(result['predicted_x']-result['x'])>4


def test_masked_core_is_not_rescued_by_external_expectation(flc,options):
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:
        world=WCS(hd['SCI'].header,hd).all_pix2world([[128,128]],0)[0];hd['DQ',7].data[125:132,125:132]=4096
    result=i.associate(p,candidate(*world),options['source_identification'],options['relative']['detection'])
    assert result['status']!='valid_detection' and result['status'] in ('detected_rejected','covered_undetected','ambiguous_association')


def test_duplicate_catalog_ids_merge_but_legacy_psf_is_not_declared_star(options):
    gaia=Table(dict(source_id=['g1'],ra=[10.],dec=[20.]))
    legacy=Table(dict(source_id=['l1','l2'],ra=[10.,10.001],dec=[20.,20.],type=['PSF','REX'],release=['DR9-North']*2))
    out=i.merge_candidates(gaia,legacy,options['source_identification'])
    assert len(out)==2 and out['legacy_id'][0]=='l1' and out['gaia_id'][0]=='g1'
    assert out['classification'][1]=='compact_galaxy_candidate' and 'DR9-North' in out['release'][0]


def test_multiple_plausible_hst_sources_are_ambiguous(flc,options):
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:
        y,x=np.indices((256,256));hd['SCI'].data+=100*np.exp(-((x-143)**2+(y-128)**2)/4)
        world=WCS(hd['SCI'].header,hd).all_pix2world([[128,128]],0)[0]
    result=i.associate(p,candidate(*world),options['source_identification'],options['relative']['detection'])
    assert result['status']=='ambiguous_association',result


def test_shared_identity_defeats_wrong_histogram_peak_and_ignores_external_ra(options):
    a=np.array([[0.,0],[100,0],[0,100],[100,100]])
    b=a+[3,2]
    ref=Table(dict(common_source_id=['a','b','c','d'],catalog_ra=[10.]*4))
    src=Table(dict(common_source_id=['d','c','b','a'],catalog_ra=[200.]*4))
    ri,ii=catalog_match(a,b[::-1],ref,src,options['relative'],True)
    assert np.allclose(a[ri]-b[::-1][ii],[-3,-2])


def test_network_failure_preserves_source_only_fallback(flc,tmp_path,monkeypatch,options):
    p=flc(size=128);paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    def fail(*a,**kw):
        from hstpipe.uvis.astrometry.catalogs import CatalogQueryError
        raise CatalogQueryError('offline')
    monkeypatch.setattr(i,'query_gaia_dr3',fail);monkeypatch.setattr(i,'query_legacy_dr10',fail)
    result=i.identify([p],paths,options)
    assert result['status']=='no_verified_common_sources' and all(q['status']=='query_failed' for q in result['queries'])
    assert not result['external_positions_used_for_relative_fit']
    # Reuse the independently cached failed-identification record offline.
    monkeypatch.setattr(i,'query_gaia_dr3',lambda *a,**kw:pytest.fail('Cached identification should resume'))
    i.identify([p],paths,options)


def test_external_positions_for_relative_fit_are_forbidden():
    with pytest.raises(ValueError,match='never relative'):
        settings({'astrometry':{'source_identification':{'use_external_positions_for_relative_fit':True}}})


def test_common_source_need_not_be_present_in_every_exposure(flc,tmp_path,monkeypatch,options):
    files=[flc(n,size=128) for n in ('a','b','c')];paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    legacy=Table(dict(source_id=['l'],ra=[10.],dec=[20.],type=['PSF'],release=['DR10-South']))
    monkeypatch.setattr(i,'query_gaia_dr3',lambda *a,**kw:Table())
    monkeypatch.setattr(i,'query_legacy_dr10',lambda *a,**kw:legacy)
    monkeypatch.setattr(i,'legacy_at_epoch',lambda t,e:t)
    def detection(p,c,*a,**kw):
        return dict(hst_source_id=str(a[0]['hst_source_id'][0]),source_id=c['source_id'],exposure=p.name,status='not_covered' if p.name.startswith('c') else 'valid_detection',x=64.,y=64.,extver=7,
            classification='legacy_point_candidate',gaia_id='',legacy_id='l',width_pix=1.3,ellipticity=.1,flux_rate=1.,snr=50.,centroid_error_pix=.1)
    monkeypatch.setattr(i,'associate_full_field',detection)
    report=i.identify(files,paths,options)
    assert len(report['catalogs'][0])==len(report['catalogs'][1])==1 and len(report['catalogs'][2])==0
    assert report['verified_common_ids']==['Ll']


def test_hst_only_native_fit_has_no_external_coordinate_constraints(flc,tmp_path,monkeypatch,options):
    from hstpipe.uvis.astrometry import relative_engine as e
    from hstpipe.uvis.astrometry.centroids import validate_centroids
    files=[flc(n,size=256) for n in ('a','b')];xy=np.array([[20,20],[100,20],[230,20],[20,100],[230,100],[20,230],[100,230],[230,230]],float)
    cats=[]
    for j,p in enumerate(files):
        with fits.open(p,mode='update') as hd:
            y,x=np.indices(hd['SCI'].shape)
            for xx,yy in xy:hd['SCI'].data+=200*np.exp(-((x-xx)**2+(y-yy)**2)/3.)
            hd['SCI'].header['CRPIX1']+=3*j
        c=Table(dict(x=xy[:,0],y=xy[:,1],extver=np.full(8,7),common_source_id=[str(k) for k in range(8)],catalog_ra=np.full(8,10+50*j),snr=np.full(8,100.),width_pix=np.full(8,1.3),minor_width_pix=np.full(8,1.3)))
        cats.append(validate_centroids(p,c,options['relative']['detection']))
    shifts=[];monkeypatch.setattr(e,'apply_native_fit',lambda image,reference,fit:shifts.append(fit['shift']))
    e.native_catalog_fit(files,cats,options['relative'])
    assert len(shifts)==1 and np.allclose(shifts[0],[3,0],atol=.001)


def test_source_selection_changes_invalidate_identification_receipt(flc,tmp_path,monkeypatch,options):
    files=[flc()];paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');calls=[]
    def missing(*a,**kw):calls.append(1);raise RuntimeError('offline')
    monkeypatch.setattr(i,'query_gaia_dr3',missing);monkeypatch.setattr(i,'query_legacy_dr10',missing)
    i.identify(files,paths,options);assert len(calls)==2
    i.identify(files,paths,options);assert len(calls)==2
    options['relative']['detection']['snr']=6
    i.identify(files,paths,options);assert len(calls)==4


def test_inconsistent_pairwise_closure_is_reported(flc,tmp_path,monkeypatch,options):
    from hstpipe.uvis.astrometry import relative_engine as e
    files=[flc(n) for n in ('a','b','c')];catalogs=[Table(dict(x=[20,40,60,80],y=[20,40,60,80],ra=[10.]*4,dec=[20.]*4,centroid_error_pix=[.1]*4,extver=[7]*4)) for p in files]
    monkeypatch.setattr(e,'match',lambda *a,**kw:(np.arange(4),np.arange(4)))
    shifts=iter(([0,0],[4,0],[0,0]))
    monkeypatch.setattr(e,'quality',lambda *a,**kw:dict(status='refinement_supported',matched=4,shift=next(shifts),transform_uncertainty_pix=.1))
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('intermediate','alignment',create=True)
    result=e.overlap_graph(files,catalogs,options['relative'],paths,'closure')
    assert result['closure']['status']=='inconsistent'


def test_verified_partial_raw_cache_retains_coverage_provenance(tmp_path):
    from hstpipe.uvis.astrometry import catalogs as c
    raw=Table(dict(ls_id=[1],ra=[10.],dec=[20.],type=['PSF'],brick_primary=[True],maskbits=[0],flux_r=[100.],flux_z=[100.],flux_ivar_r=[1.],flux_ivar_z=[1.],nobs_r=[1],release=[9011]))
    q=c.legacy_adql(10,20,.01,raw.colnames)
    c._cached('legacy_dr10_raw',dict(query=q,endpoint=c.TAP_URL),lambda:(raw,dict(query=q,endpoint=c.TAP_URL)),tmp_path,False,False)
    t=c.cached_footprint('legacy_dr10',10,20,.02,cache_dir=tmp_path)
    assert len(t)==1 and not t.meta['complete_query_coverage'] and t['release'][0]=='DR9-North'
    assert t.meta['cache_fallback'] and t.meta['original_query']==q


def test_masked_gaia_linked_psf_is_identification_prior_not_absolute_reference():
    from hstpipe.uvis.astrometry.catalogs import standardize_legacy
    raw=Table(dict(ls_id=[1],ra=[10.],dec=[20.],type=['PSF'],brick_primary=[True],maskbits=[51216],flux_r=[600.],flux_z=[800.],flux_ivar_r=[10.],flux_ivar_z=[3.],fracmasked_r=[.028],fracmasked_z=[.535],fracflux_r=[0.],fracflux_z=[0.],ref_cat=['GE'],ref_id=[123],nobs_r=[4],release=[10000]))
    assert len(standardize_legacy(raw))==0
    priors=standardize_legacy(raw,dict(allow_masked_stellar_priors=True))
    assert len(priors)==1 and priors['identification_prior_only'][0]


def test_saturated_star_is_reported_and_detached_halo_not_used(flc,options):
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:
        world=WCS(hd['SCI'].header,hd).all_pix2world([[128,128]],0)[0]
        hd['DQ',7].data[127:130,127:130]=256|4096
    result=i.associate(p,candidate(*world),options['source_identification'],options['relative']['detection'])
    assert result['status']=='saturated_source_in_window' and np.isnan(result['x'])
