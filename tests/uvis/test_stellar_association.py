"""Focused image-based counterpart resolution; no live network queries."""
import numpy as np
import pytest
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from hstpipe.uvis.astrometry.identification import associate_full_field
from hstpipe.uvis.config import settings
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity

def scene(flc,stars=((64.,64.,5000.),),cosmics=()):
    from photutils.psf import CircularGaussianPRF
    p=flc(size=128);yy,xx=np.indices((128,128));a=np.full((128,128),5.)
    for x,y,flux in stars:a+=CircularGaussianPRF(x_0=x,y_0=y,flux=flux,fwhm=2.2)(xx,yy)
    for x,y in cosmics:a+=1000*np.exp(-((xx-x)**2+(yy-y)**2)/.2)
    error=np.sqrt(a+4);a+=np.random.default_rng(17).normal(0,error)
    with fits.open(p,mode='update') as hd:hd['SCI',7].data=a;hd['ERR',7].data=error
    coords=list(stars)+[(x,y,1000.) for x,y in cosmics]
    cat=Table(dict(x=[c[0] for c in coords],y=[c[1] for c in coords],extver=[7]*len(coords),centroid_error_pix=[.1]*len(coords),hst_source_id=['HST'+str(i) for i in range(len(coords))]))
    with fits.open(p) as hd:
        from hstpipe.uvis.astrometry.wcs import detector_wcs
        ra,dec=detector_wcs(hd,7).all_pix2world([[64,64]],0)[0]
    prior=dict(source_id='Gaia-test',gaia_id='123',legacy_id='',classification='gaia_stellar_candidate',ra=ra,dec=dec,release='DR3')
    return p,cat,prior

def resolve(p,cat,prior):
    options=settings({})['astrometry']
    return associate_full_field(p,prior,cat,options['source_identification'],detection=options['relative']['detection'])

def test_star_among_cosmic_candidates_is_resolved_without_nearest_or_brightest(flc):
    p,c,g=scene(flc,stars=((72.,67.,5000.),),cosmics=((64.,64.),(45.,45.)))
    before=identity(p);r=resolve(p,c,g)
    assert r['status']=='valid_detection' and abs(r['x']-72)<.15 and abs(r['y']-67)<.15
    assert r['psf_audit']['raw_candidates']==3 and r['psf_audit']['rejected']==2
    assert identity(p)==before

def test_two_equally_credible_stars_remain_ambiguous(flc):
    p,c,g=scene(flc,stars=((55.,64.,5000.),(74.,64.,5000.)))
    r=resolve(p,c,g);assert r['status']=='ambiguous_stellar_profiles' and not r.get('psf_selected')

def test_crowded_window_is_fitted_instead_of_rejected_by_count(flc):
    stars=tuple((float(x),float(y),5000.) for y in (36,50,64,78,92) for x in (36,50,64,78,92))
    p,c,g=scene(flc,stars=stars);r=resolve(p,c,g)
    assert r['psf_audit']['raw_candidates']==25
    assert len(r['psf_audit']['fits'])==25 and r['status']=='ambiguous_stellar_profiles'

def test_saturated_neighbor_does_not_veto_unmasked_stellar_core(flc):
    p,c,g=scene(flc)
    with fits.open(p,mode='update') as hd:hd['DQ',7].data[40,40]=256
    r=resolve(p,c,g);assert r['status']=='valid_detection'

def test_manual_local_detector_uses_same_psf_resolution(flc):
    from hstpipe.uvis.astrometry.identification import associate
    p,c,g=scene(flc,stars=((72.,67.,5000.),),cosmics=((64.,64.),(45.,45.)))
    a=settings({})['astrometry'];r=associate(p,g,a['source_identification'],a['relative']['detection'])
    assert r['status']=='valid_detection' and abs(r['x']-72)<.2

@pytest.mark.parametrize('bit',[8192,256])
def test_masked_core_is_not_recovered_by_psf_interpolation(flc,bit):
    p,c,g=scene(flc)
    with fits.open(p,mode='update') as hd:hd['DQ',7].data[63:66,63:66]=bit
    r=resolve(p,c,g);assert r['status']=='no_reliable_stellar_profile'

def test_inherited_4096_reassessment_is_preserved(flc):
    p,c,g=scene(flc)
    with fits.open(p,mode='update') as hd:hd['DQ',7].data[64,64]=4096
    before=identity(p);r=resolve(p,c,g)
    assert r['status']=='valid_detection' and identity(p)==before

def test_failed_check_window_is_not_used_as_centroid_validation(flc,monkeypatch):
    from photutils.psf import PSFPhotometry
    original=PSFPhotometry.__call__;calls=[]
    def failed_check(self,*args,**kwargs):
        result=original(self,*args,**kwargs);calls.append(True)
        if len(calls)==2:result['flags']|=8
        return result
    monkeypatch.setattr(PSFPhotometry,'__call__',failed_check)
    p,c,g=scene(flc);r=resolve(p,c,g)
    assert r['status']=='no_reliable_stellar_profile'

def test_overlapping_stars_are_fitted_jointly_without_forced_identity(flc):
    p,c,g=scene(flc,stars=((62.,64.,5000.),(66.,64.,4000.)))
    r=resolve(p,c,g);assert any(q.get('group_size')==2 for q in r['psf_audit']['fits'])
    assert r['status']!='valid_detection'

def test_cross_exposure_profile_discrimination_does_not_use_coordinates():
    from hstpipe.uvis.astrometry.stellar import choose_profile
    anchor=dict(flux=5000.,flux_error=70.,flux_rate=50.,fwhm_pix=2.)
    a=dict(anchor,x=500.,y=-100.);b=dict(anchor,flux_rate=15.,x=0.,y=0.)
    assert choose_profile([a,b],anchor,{}) is a
    assert choose_profile([a,dict(a)],anchor,{}) is None

def test_duplicate_detector_rows_are_merged_before_psf_association():
    from hstpipe.uvis.astrometry.full_field import merge_detectors
    c=Table(dict(x=[10.],y=[20.],extver=[1],centroid_error_pix=[.1],source_id=['a'],exposure_rootname=['test']))
    second=c.copy();second['x']+=.1;second['source_id']=['b']
    assert len(merge_detectors([c,second]))==1

def test_field_star_fallback_is_local_and_does_not_update_working_wcs(flc,tmp_path):
    from hstpipe.uvis.astrometry.target_fallback import target_alignment
    from hstpipe.uvis.astrometry.centroids import QUALITY_VERSION
    inputs=[flc(name='a',size=256),flc(name='b',size=256)];cats=[]
    for i,p in enumerate(inputs):
        c=Table(dict(x=[25.+.2*i],y=[40.],extver=[7],centroid_error_pix=[.08],common_source_id=['field-star'],detection_method=['PSFPhotometry']))
        c.meta['centroid_quality']=dict(version=QUALITY_VERSION);cats.append(c)
    before=[identity(p) for p in inputs];opts=settings({})['astrometry']['relative'];paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');paths.area('intermediate','alignment',create=True)
    r,_=target_alignment(inputs,paths,opts,catalogs=cats)
    assert r['status']=='local_only' and r['result_state']=='local_shift_only'
    assert r['reference_source']['source_id']=='field-star' and not r['full_field_validated']
    assert r['independent_scatter_pix'] is None and [identity(p) for p in inputs]==before

def test_psf_policy_changes_invalidate_association_but_not_background_or_detection(flc,tmp_path):
    from hstpipe.uvis.astrometry.identification import identification_key
    from hstpipe.uvis.astrometry.full_field import detection_key
    from hstpipe.uvis.inventory import discover
    from hstpipe.uvis.storage import stage_keys
    p=flc();paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');cache=paths.cache()
    a=settings({});b=settings({'astrometry':{'source_identification':{'psf_resolution':{'max_qfit':.4}}}})
    assert identification_key(cache,[p],a['astrometry'])!=identification_key(cache,[p],b['astrometry'])
    assert detection_key(cache,[p],a['astrometry']['relative'])==detection_key(cache,[p],b['astrometry']['relative'])
    rows=discover(tmp_path)
    ak=stage_keys(cache,rows,a,[]);bk=stage_keys(cache,rows,b,[])
    assert ak[0]==bk[0] and ak[1]!=bk[1]


@pytest.mark.parametrize('field', ['flux_err', 'qfit', 'x_err'])
def test_rejected_nonfinite_psf_diagnostics_remain_json_serializable(flc,tmp_path,monkeypatch,field):
    import json
    from photutils.psf import PSFPhotometry
    from hstpipe.uvis.astrometry.stellar import diagnostic
    original=PSFPhotometry.__call__
    def undefined_measurement(self,*args,**kwargs):
        result=original(self,*args,**kwargs);result[field]=np.full(len(result),np.nan)
        return result
    monkeypatch.setattr(PSFPhotometry,'__call__',undefined_measurement)
    p,c,g=scene(flc);before=identity(p);r=resolve(p,c,g)
    assert r['status']=='no_reliable_stellar_profile' and not r.get('psf_selected')
    audit=r['psf_audit'];assert audit['rejected']==1
    row=audit['fits'][0]
    assert 'nonfinite_or_unconstrained_fit' in row['reason']
    assert row['nonfinite_fields'] and any(row[key] is None for key in row['nonfinite_fields'])
    json.dumps(audit,allow_nan=False)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    diagnostic(p,r,paths,'nonfinite_test')
    saved=json.loads(Path(r['psf_audit_file']).read_text())
    assert saved==audit and identity(p)==before
