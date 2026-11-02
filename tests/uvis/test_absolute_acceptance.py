"""Focused external-stage tests; no FLC relative registration or network."""
import numpy as np
from astropy.table import Table
from hstpipe.uvis.astrometry.absolute_matching import validate_absolute as validate,reference_audit
from hstpipe.uvis.config import settings

def sample(n=5,shift=(.1,-.05)):
    x=np.linspace(10,240,n);y=np.linspace(20,230,n)
    source=Table(dict(x=x,y=y,ra=10+x/3600,dec=20+y/3600,centroid_error_pix=np.full(n,.08)))
    ref=source.copy();ref['source_id']=np.arange(n).astype(str)
    ref['ra']+=shift[0]/3600/np.cos(np.deg2rad(20));ref['dec']+=shift[1]/3600
    return source,ref

def opts():return dict(min_matches=1,fitgeometry='shift',search_radius_arcsec=2,clip_sigma=3,max_rms_arcsec=.2,max_offset_arcsec=.1,systematic_floor_arcsec=.01,min_span_fraction=.1)

def test_shift_minimum_is_not_six():
    assert settings({'astrometry':{'absolute':{'min_matches':1}}})['astrometry']['absolute']['min_matches']==1
    a,b=sample(2,shift=(0,0));r,t=validate(a,b,opts(),(256,256))
    assert r['status']=='validated' and r['matched']==2 and r['independent_validation_sources']==0
    assert r['leave_one_out_max_arcsec']<1e-6

def test_one_reference_diagnostic_only():
    a,b=sample(1);r,t=validate(a,b,opts(),(256,256))
    assert r['status']=='diagnostic_single_reference' and r['postfit_rms_arcsec']<1e-8

def test_initial_offset_is_not_postfit_scatter():
    a,b=sample(5,shift=(.5,-.3));r,t=validate(a,b,opts(),(256,256))
    assert r['rms_arcsec']>.5 and r['postfit_rms_arcsec']<1e-3
    assert r['status']!='validated'

def test_false_legacy_association_does_not_inflate_clip_scale():
    a,b=sample(8,shift=(0,0));b['ra'][-1]+=.0003;b['dec'][-2]+=.0004
    r,t=validate(a,b,opts(),(256,256))
    assert r['rejected']==2 and r['matched']==6 and r['status']=='validated'
    assert r['initial_associations']==8

def test_reference_coverage_is_not_cone_membership(flc,tmp_path):
    p=flc(size=256);from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(p) as hd:w=WCS(hd['SCI'].header,hd)
    xy=np.array([[40.,40.],[500.,500.]]);world=w.all_pix2world(xy,0)
    ref=Table(dict(source_id=['inside','outside'],ra=world[:,0],dec=world[:,1]));sources=Table(dict(x=[40.],y=[40.],ra=[world[0,0]],dec=[world[0,1]]))
    audit=reference_audit(p,sources,ref,opts())
    assert list(audit['covered'])==[True,False]
    assert audit['association_status'][1]=='outside_coverage'

def test_native_stsci_two_reference_fit_and_wcs_preservation(flc,tmp_path,monkeypatch):
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry import alignment as a
    from hstpipe.uvis.inventory import discover
    from astropy.io import fits
    from astropy.wcs import WCS
    p=flc(size=256)
    with fits.open(p,mode='update') as hd:
        h=hd['SCI'].header
        h['CD1_1']=h['CDELT1'];h['CD1_2']=0.;h['CD2_1']=0.;h['CD2_2']=h['CDELT2']
        del h['CDELT1'];del h['CDELT2']
    with fits.open(p) as hd:w=WCS(hd['SCI'].header,hd)
    xy=np.array([[40.,40.],[210.,210.]]);sky=w.all_pix2world(xy,0)
    refs=Table(dict(source_id=['a','b'],ra=sky[:,0]+.02/3600/np.cos(np.deg2rad(20)),dec=sky[:,1]))
    opts=settings({})['astrometry'];paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    def detect(image,*args,**kwargs):
        with fits.open(image) as hd:world=WCS(hd['SCI'].header,hd).all_pix2world(xy,0)
        return Table(dict(x=xy[:,0],y=xy[:,1],ra=world[:,0],dec=world[:,1],flux=[100.,100.]))
    monkeypatch.setattr(a,'detect_sources',detect);monkeypatch.setattr(a,'query_gaia_dr3',lambda *args,**kwargs:refs)
    monkeypatch.setattr(a,'query_legacy_dr10',lambda *args,**kwargs:refs[:0])
    monkeypatch.setattr(a,'cached_footprint',lambda *args,**kwargs:(_ for _ in ()).throw(a.CatalogQueryError('no fixture cache')))
    solution,_=a.absolute_alignment(p,discover(tmp_path),paths,opts)
    assert solution['status']=='validated',solution
    assert solution['validation']['matched']==2
    assert solution['validation']['rms_arcsec']<.01
    assert solution['original_wcs_key'] and solution['transformations']

def test_absolute_report_distinguishes_failure_categories():
    from hstpipe.uvis.reporting.report import _absolute_summary
    statuses=['insufficient_verified_references','software_match_count_restriction','ambiguous_associations','excessive_postfit_residuals','external_catalog_unavailable','validated']
    record={'scientific_assessment':{'astrometry':{'attempts':[dict(reference='gaia_dr3',status=s) for s in statuses]}}}
    output=_absolute_summary(record)
    for phrase in ('Insufficient verified','Software match-count','ambiguous associations','Excessive post-fit','External catalog unavailable','Validated absolute'):
        assert phrase in output
    assert 'Native fit used' in output and 'Independent checks' in output


def test_worse_finite_absolute_correction_preserves_acceptable_wcs(flc,tmp_path,monkeypatch):
    from test_astrometry import preparation,grid,fake_tweakreg
    from hstpipe.uvis.astrometry import alignment as a
    import drizzlepac.tweakreg
    from astropy.io import fits
    p,rows,paths,refs,opts=preparation(flc,tmp_path)
    monkeypatch.setattr(a,'detect_sources',lambda image,*args,**kwargs:grid(image))
    monkeypatch.setattr(a,'query_gaia_dr3',lambda *args,**kwargs:refs)
    def worse(files,**kwargs):
        fake_tweakreg(files,**kwargs)
        with fits.open(files[0],mode='update') as hd:hd['SCI'].header['CRVAL1']+=.03/3600
    monkeypatch.setattr(drizzlepac.tweakreg,'TweakReg',worse)
    before=a.signature(p);solution,_=a.absolute_alignment(p,rows,paths,opts)
    assert solution['status']=='validated' and not solution['correction_applied']
    assert solution['attempts'][0]['rejected_transformation']
    assert a.signature(p)==before
