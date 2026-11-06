"""Focused cross-filter graph, WCS application, and transfer checks."""
import json
import shutil
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from hstpipe.common.paths import ReductionPaths
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry import cross_filter as c
from hstpipe.uvis.astrometry.wcs import detector_wcs


def catalog(image,options,paths=None):
    xy=np.array([[20,20],[90,20],[210,20],[20,100],[210,100],[20,220],[110,220],[220,220]],float)
    with fits.open(image) as hd:
        v=hd['SCI'].header.get('EXTVER',1);w=detector_wcs(hd,v);sky=w.all_pix2world(xy,0)
    t=Table(dict(x=xy[:,0],y=xy[:,1],ra=sky[:,0],dec=sky[:,1],flux=np.ones(len(xy)),extver=np.full(len(xy),v),centroid_error_pix=np.full(len(xy),.05)))
    if paths:
        p=paths.file('intermediate','alignment','cross_filter_sources','ecsv');p.parent.mkdir(parents=True,exist_ok=True);t.write(p,format='ascii.ecsv',overwrite=True)
    return t


def build(flc,tmp_path,monkeypatch,offset):
    records={};paths={};files=[]
    for key,shift in [('F467M',0),('F621M',offset)]:
        p=ReductionPaths(tmp_path/'redux','TEST','uvis',key);paths[key]=p
        image=p.file('intermediate','drizzle','provisional','fits');image.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(flc(key,size=256),image)
        with fits.open(image,mode='update') as hd:
            h=hd['SCI'].header;h['CRPIX1']+=shift;h['CD1_1']=h.pop('CDELT1');h['CD2_2']=h.pop('CDELT2');h['CD1_2']=0.;h['CD2_1']=0.
        files.append(image);records[key]=dict(filter=key,setup=None,scientific_assessment=dict(relative_astrometry=dict(status='validated')))
    (p.directory/'manifest.json').write_text(json.dumps(dict(filters=records)))
    monkeypatch.setattr(c,'cross_sources',catalog)
    return paths,files


def test_connections_allow_bridge_and_keep_disconnected_components():
    def edge(n):return dict(status='supported',before=dict(robust_inliers=n,transform_uncertainty_pix=.05))
    edges={('A','B'):edge(4),('B','A'):edge(4),('B','C'):edge(3),('C','B'):edge(3)}
    plan=c.connection_plan(['A','B','C','D'],edges,dict(A=20,B=10,C=8,D=100),requested='A')
    assert plan['C']['anchor_group']=='B' and plan['C']['depth']==2
    assert plan['C']['network_root']=='A' and plan['D']['connected_groups']==['D']


@pytest.mark.parametrize('offset',[0.,3.])
def test_native_correction_preserves_original_and_cache(flc,tmp_path,monkeypatch,offset):
    paths,files=build(flc,tmp_path,monkeypatch,offset);original=[p.read_bytes() for p in files]
    out=c.align_cross_filters(paths['F467M'].directory,options={'cross_filter':{'reference_filter':'F467M'}})
    assert out['F621M']['status']=='validated'
    assert out['F621M']['validation']['current_rms_pix']<1e-5
    assert [p.read_bytes() for p in files]==original
    monkeypatch.setattr(c,'assess_pair',lambda *a:pytest.fail('Unchanged internal pair should be cached'))
    c.align_cross_filters(paths['F467M'].directory,options={'cross_filter':{'reference_filter':'F467M'}})


def test_fixed_scale_rotation_uses_stsci_model_and_predictive_check(flc,tmp_path):
    a=flc('a',size=256);b=flc('b',size=256)
    ca=catalog(a,{});cb=ca.copy();theta=np.deg2rad(.08)
    xy=np.column_stack([cb['x'],cb['y']]);m=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    xy=(xy-128)@m.T+128+[2,-1];cb['x']=xy[:,0];cb['y']=xy[:,1]
    o=settings({})['astrometry']['relative'];ca['centroid_error_pix']=.005;cb['centroid_error_pix']=.005
    r,_,_=c.assess_pair(a,b,ca,cb,o)
    assert r['status']=='supported' and r['before']['geometry']=='rshift'
    assert abs(r['before']['rotation_deg']+.08)<1e-4
    assert r['before']['loo_rms_pix']<1e-4


def test_single_source_is_never_independent_cross_filter_validation(flc):
    a=flc('a');b=flc('b');cat=catalog(a,{})[:1]
    r,_,_=c.assess_pair(a,b,cat,cat,settings({})['astrometry']['relative'])
    assert r['status']=='unvalidated'


def test_saved_pairs_transfer_complete_wcs_and_archive_original(flc,tmp_path,monkeypatch):
    from hstpipe.common.inventory import Exposure
    from hstpipe.common.provenance import identity
    from hstpipe.uvis.astrometry.wcs import signature
    from stwcs.wcsutil import altwcs
    paths,files=build(flc,tmp_path,monkeypatch,3.)
    c.align_cross_filters(paths['F467M'].directory,options={'cross_filter':{'reference_filter':'F467M'}})
    ref,image=files;original_signature=signature(image)
    # Simulate the parent external correction without changing its archived
    # original drizzle frame. The child must inherit this correction as well.
    with fits.open(ref) as hd:v=hd['SCI'].header['EXTVER']
    k,_=altwcs.archive_wcs(str(ref),ext=[('SCI',v)],wcsname='HP-PROVISIONAL',mode=altwcs.ArchiveMode.AUTO_RENAME)
    with fits.open(ref,mode='update') as hd:hd['SCI'].header['HPORIGK']=k;hd['SCI'].header['CRVAL1']+=1e-5
    rows=[Exposure(image,dict(EXPSTART=58000.),'flc',identity(image))]
    options=settings({})['astrometry']
    anchor=(ref,dict(reference='gaia_dr3',provisional_reference=str(ref),reference_releases=['DR3'],final_validation=dict(status='validated',rms_arcsec=.01,uncertainty_arcsec=.01)))
    solution,refs=c.transfer_absolute_frame(image,rows,paths['F621M'],options,anchor)
    assert solution['status']=='validated' and solution['reference_mode']=='cross_filter'
    assert signature(image)==original_signature
    assert solution['validation']['rms_arcsec']<1e-5
    assert solution['validation']['fixed_pair_validation']['current_rms_pix']<1e-4
    assert solution['original_wcs_key'] and len(refs)==8


def test_stellar_fallback_keeps_measured_centroids_and_no_external_constraints(flc,tmp_path,monkeypatch):
    paths,files=build(flc,tmp_path,monkeypatch,3.)
    def sparse(image,options,paths=None):
        t=catalog(image,options)[:1]
        if paths:
            p=paths.file('intermediate','alignment','cross_filter_sources','ecsv');p.parent.mkdir(parents=True,exist_ok=True);t.write(p,format='ascii.ecsv',overwrite=True)
        return t
    monkeypatch.setattr(c,'cross_sources',sparse)
    monkeypatch.setattr(c,'detect_sources',lambda image,*a,**kw:catalog(image,{}))
    result=c.align_cross_filters(paths['F467M'].directory,options={'cross_filter':{'reference_filter':'F467M'}})
    r=result['F621M']
    assert r['status']=='validated' and r['source_selection']=='DAO_stellar_fallback'
    pairs=Table.read(r['matches'],format='ascii.ecsv')
    assert np.allclose(pairs['source_x'],catalog(files[1],{})['x'][np.asarray(pairs['source_index'],int)])
    assert r['compact_attempt']['status']=='unvalidated'


def test_transfer_can_use_bridge_and_its_new_provisional(tmp_path):
    from hstpipe.uvis.astrometry.alignment import _anchor
    paths=ReductionPaths(tmp_path/'redux','T','uvis','F275W')
    receipt=paths.file('intermediate','alignment','cross_filter');receipt.parent.mkdir(parents=True,exist_ok=True)
    receipt.write_text(json.dumps(dict(status='validated',anchor_group='F336W',network_root='F621M',matches='saved_pairs.ecsv')))
    product=tmp_path/'parent_final.fits';product.touch()
    new=tmp_path/'parent_corrected_provisional.fits';new.touch()
    record=dict(filter='F336W',products=[dict(path=str(product))],scientific_assessment=dict(
        relative_astrometry=dict(status='validated'),cross_filter_alignment=dict(status='validated',anchor_group='F621M'),
        provisional_image=dict(path=str(tmp_path/'obsolete_provisional.fits')),
        astrometry=dict(status='validated',observation_epoch=2020.,provisional_image_path=str(new))))
    (paths.directory/'manifest.json').write_text(json.dumps(dict(filters=dict(F336W=record))))
    options=settings({'astrometry':{'cross_filter':{'reference_filter':'F621M'}}})['astrometry']
    anchor=_anchor(paths,2020.,options)
    assert anchor[0]==product and anchor[1]['provisional_reference']==str(new)


def test_existing_final_retention_requires_measured_accuracy(flc,tmp_path,monkeypatch):
    from hstpipe.common.provenance import identity
    paths,files=build(flc,tmp_path,monkeypatch,.1)
    c.align_cross_filters(paths['F467M'].directory,options={'cross_filter':{'reference_filter':'F467M'}})
    ref,image=files
    with fits.open(image,mode='update') as hd:hd[0].header['HPQASTAT']='REVIEWED';hd['SCI'].header['ASTSTAT']='VALIDATED'
    anchor=(ref,dict(reference='gaia_dr3',provisional_reference=str(ref)))
    before=identity(image);opts=settings({})['astrometry']
    result=c.validate_existing_frame(image,paths['F621M'],opts,anchor)
    assert result['status']=='validated' and result['preserved_without_resampling']
    assert identity(image)==before
    with fits.open(image,mode='update') as hd:hd['SCI'].header['CRPIX1']+=10
    result=c.validate_existing_frame(image,paths['F621M'],opts,anchor)
    assert result['status']!='validated'
