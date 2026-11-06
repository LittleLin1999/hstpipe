"""Focused tests of the internal-registration/external-catalog barrier."""
import json
from pathlib import Path
import pytest
from hstpipe.common.config import ReductionConfig
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry import orchestration as o


def test_all_selected_provisionals_precede_any_absolute_stage(tmp_path,monkeypatch):
    root=tmp_path/'redux';directory=root/'T/uvis';directory.mkdir(parents=True)
    cfg=ReductionConfig(str(tmp_path/'raw'),str(root),'uvis',options=settings({}),resume=True)
    plans=[dict(configuration=dict(target='T',element=f,setup=None)) for f in ('F467M','F621M')]
    events=[]
    state=dict(filters={f:dict(scientific_assessment=dict(alignment=dict(status='validated'))) for f in ('F467M','F621M')})
    (directory/'manifest.json').write_text(json.dumps(state))
    def engine(c,*a,prepare_only=False):events.append(('prepare' if prepare_only else 'external',c.element))
    def cross(*a):
        assert events==[('prepare','F467M'),('prepare','F621M')]
        events.append(('cross',None));return {'F621M':dict(status='anchor'),'F467M':dict(status='validated')}
    def run(c,worker,dry_run=False,**kw):
        if dry_run:return plans
        worker(c,None,None)
        return [directory]
    from hstpipe.uvis import inventory,reporting
    monkeypatch.setattr(reporting,'generate_report',lambda *a:None)
    monkeypatch.setattr(inventory,'scan',lambda *a:([],[]))
    for p in plans:p['group']=['T',p['configuration']['element'],None]
    monkeypatch.setattr(o,'internal_cross_filter',cross);monkeypatch.setattr(o,'run_uvis',run)
    assert o.relative_then_absolute(cfg,engine)==[directory]
    assert events==[('prepare','F467M'),('prepare','F621M'),('cross',None),('external','F621M'),('external','F467M')]


def test_unvalidated_relative_groups_never_enter_absolute_stage(tmp_path,monkeypatch):
    root=tmp_path/'redux';d=root/'T/uvis';d.mkdir(parents=True)
    (d/'manifest.json').write_text(json.dumps(dict(filters={'F467M':dict(scientific_assessment=dict(alignment=dict(status='unvalidated')))})))
    cfg=ReductionConfig(str(tmp_path/'raw'),str(root),'uvis',options=settings({}))
    plan=[dict(configuration=dict(target='T',element='F467M',setup=None))]
    from hstpipe.uvis import inventory,reporting
    monkeypatch.setattr(reporting,'generate_report',lambda *a:None)
    monkeypatch.setattr(inventory,'scan',lambda *a:([],[]))
    plan[0]['group']=['T','F467M',None]
    def run(c,worker,dry_run=False,**kw):
        if dry_run:return plan
        assert worker is not None, 'Unvalidated filter entered absolute finalization'
        return [d]
    monkeypatch.setattr(o,'internal_cross_filter',lambda *a:{'F467M':dict(status='skipped')})
    monkeypatch.setattr(o,'run_uvis',run)
    assert o.relative_then_absolute(cfg,None)==[d]


def test_relative_external_catalog_configuration_is_rejected():
    with pytest.raises(ValueError,match='external'):
        settings({'astrometry':{'relative':{'external_catalog':True}}})
    with pytest.raises(ValueError,match='provisional'):
        settings({'astrometry':{'absolute':{'query_after_provisional_drizzle':False}}})


def _internal_catalog(image,*args,**kwargs):
    import numpy as np
    from astropy.table import Table
    from astropy.io import fits
    from astropy.wcs import WCS
    xy=np.array([[20,20],[100,20],[230,20],[20,100],[230,100],[20,230],[100,230],[230,230]],float)
    with fits.open(image) as hd:
        world=WCS(hd['SCI'].header,hd).all_pix2world(xy,0)
        version=hd['SCI'].header.get('EXTVER',1)
    from hstpipe.uvis.astrometry.centroids import QUALITY_VERSION
    c=Table(dict(x=xy[:,0],y=xy[:,1],ra=world[:,0],dec=world[:,1],flux=np.ones(8)*100,extver=np.full(8,version),centroid_error_pix=np.ones(8)*.1,width_pix=np.ones(8),centroid_valid=np.ones(8,bool)))
    c.meta['centroid_quality']=dict(version=QUALITY_VERSION)
    if kwargs.get('paths'):
        p=kwargs['paths'].file('intermediate','alignment','cross_filter_sources','ecsv');p.parent.mkdir(parents=True,exist_ok=True);c.write(p,format='ascii.ecsv',overwrite=True)
    return c



def test_overlap_graph_exposes_disconnected_exposure(flc,tmp_path):
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry.relative_engine import overlap_graph
    files=[flc('a',size=256),flc('b',size=256),flc('isolated',size=256,ra=11)]
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    result=overlap_graph(files,[_internal_catalog(p) for p in files],settings({})['astrometry']['relative'],paths,'test')
    assert result['status']=='disconnected' and len(result['validated_components'])==2


def test_tweakback_relative_qa_checks_internal_pairs_not_external_sources(flc,tmp_path):
    import shutil
    from astropy.io import fits
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry.workflow import validate_relative_transfer
    files=[flc('a',size=256),flc('b',size=256)];updated=[]
    for p in files:
        copy=tmp_path/('copy_'+p.name);shutil.copy2(p,copy);updated.append(copy)
        with fits.open(copy,mode='update') as hd:hd['SCI'].header['CRVAL1']+=.001
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W');entries=[]
    for p in files:
        out=paths.exposure_file('intermediate','alignment',p.stem,'relative_sources','ecsv')
        out.parent.mkdir(parents=True,exist_ok=True)
        _internal_catalog(p).write(out,format='ascii.ecsv');entries.append(dict(exposure=p.name,path=str(out)))
    relative=dict(reference=files[0].name,detector_catalogs=entries)
    options=settings({})['astrometry']['relative']
    result=validate_relative_transfer(files,updated,relative,options,paths)
    assert result['status']=='validated' and not result['external_catalog_used']
    with fits.open(updated[1],mode='update') as hd:hd['SCI'].header['CRPIX1']+=4
    assert validate_relative_transfer(files,updated,relative,options,paths)['status']=='unvalidated'


@pytest.mark.parametrize('changed,disconnected',[(False,False),(True,False),(False,True)])
def test_tweakback_preserves_accepted_connected_pairs(flc,tmp_path,changed,disconnected):
    import shutil
    from astropy.io import fits
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry.workflow import validate_relative_transfer
    files=[flc(n,size=256) for n in ('a','b','c')];updated=[];entries=[]
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    for i,p in enumerate(files):
        copy=tmp_path/('copy_'+p.name);shutil.copy2(p,copy);updated.append(copy)
        with fits.open(copy,mode='update') as hd:
            hd['SCI'].header['CRVAL1']+=.001
            if changed and i==2:hd['SCI'].header['CRPIX1']+=4
        cat=_internal_catalog(p);cat['common_source_id']=[str(k) for k in range(len(cat))]
        if i==0:cat=cat[:4]
        if i==2:cat=cat[4:]
        out=paths.exposure_file('intermediate','alignment',p.stem,'relative_sources','ecsv');out.parent.mkdir(parents=True,exist_ok=True)
        cat.write(out,format='ascii.ecsv');entries.append(dict(exposure=p.name,path=str(out)))
    edges=[dict(reference=files[j].name,exposure=files[i].name,status='validated') for j,i in [(0,1),(1,2)]]
    if disconnected:edges=edges[:1]
    relative=dict(reference=files[0].name,detector_catalogs=entries,network=dict(status='connected',closure=dict(status='not_constrained'),edges=edges))
    result=validate_relative_transfer(files,updated,relative,settings({})['astrometry']['relative'],paths)
    assert result['status']==('unvalidated' if changed or disconnected else 'validated')
    if not changed and not disconnected:
        assert result['validation_basis']=='accepted_connected_graph'
        assert [(p['reference'],p['exposure']) for p in result['pairs']]==[(files[0].name,files[1].name),(files[1].name,files[2].name)]


def test_catalog_failure_does_not_block_other_filter(flc,tmp_path,monkeypatch):
    import shutil
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry import cross_filter as c,align_cross_filters
    records={}
    for filt in ('F467M','F621M'):
        paths=ReductionPaths(tmp_path/'redux','TEST','uvis',filt)
        p=paths.file('intermediate','drizzle','provisional','fits');p.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(flc(filt,size=256),p)
        records[filt]=dict(filter=filt,setup=None,scientific_assessment=dict(alignment=dict(status='validated')))
    (paths.directory/'manifest.json').write_text(json.dumps(dict(filters=records)))
    def detect(image,*a,**kw):
        if 'F467M' in image.name:raise ValueError('Bad source catalog')
        return _internal_catalog(image)
    monkeypatch.setattr(c,'cross_sources',detect)
    results=align_cross_filters(paths.directory)
    assert results['F467M']['status']=='unvalidated'
    assert results['F621M']['status']=='unvalidated'
    assert results['F621M']['source_count']>0  # still eligible for its own absolute solution
    assert not results['F621M']['external_catalog_used']


def test_absolute_anchor_survives_running_manifest_entry(tmp_path):
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry.alignment import _anchor
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    product=tmp_path/'anchor.fits';product.touch()
    stage=paths.file('intermediate','alignment','cross_filter');stage.parent.mkdir(parents=True,exist_ok=True)
    stage.write_text(json.dumps(dict(status='validated',anchor_group='F621M')))
    own=dict(status='running',filter='F275W')
    anchor=dict(filter='F621M',products=[dict(path=str(product))],scientific_assessment=dict(
        relative_astrometry=dict(status='validated'),cross_filter_alignment=dict(status='anchor',anchor_group='F621M'),
        astrometry=dict(status='validated',observation_epoch=2020.,reference='gaia_dr3')))
    (paths.directory/'manifest.json').write_text(json.dumps(dict(filters=dict(F275W=own,F621M=anchor))))
    result=_anchor(paths,2020.,settings({})['astrometry'])
    assert result is not None and result[0]==product
    stage.write_text(json.dumps(dict(status='unvalidated',anchor_group='F621M')))
    assert _anchor(paths,2020.,settings({})['astrometry']) is None
