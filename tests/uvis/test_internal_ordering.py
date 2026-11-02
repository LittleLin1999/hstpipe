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
    def run(c,engine,dry_run=False):
        if dry_run:return plans
        for f in ('F467M','F621M'):engine(f,None,None)
        return [directory]
    def engine(c,*a,prepare_only=False):events.append(('prepare' if prepare_only else 'external',c))
    def cross(*a):
        assert events==[('prepare','F467M'),('prepare','F621M')]
        events.append(('cross',None));return {'F621M':dict(status='anchor'),'F467M':dict(status='validated')}
    def finish(c,engine,**kw):engine(c.element,None,None)
    monkeypatch.setattr(o,'run',run);monkeypatch.setattr(o,'internal_cross_filter',cross);monkeypatch.setattr(o,'run_uvis',finish)
    assert o.relative_then_absolute(cfg,engine)==[directory]
    assert events==[('prepare','F467M'),('prepare','F621M'),('cross',None),('external','F621M'),('external','F467M')]


def test_unvalidated_relative_groups_never_enter_absolute_stage(tmp_path,monkeypatch):
    root=tmp_path/'redux';d=root/'T/uvis';d.mkdir(parents=True)
    (d/'manifest.json').write_text(json.dumps(dict(filters={'F467M':dict(scientific_assessment=dict(alignment=dict(status='unvalidated')))})))
    cfg=ReductionConfig(str(tmp_path/'raw'),str(root),'uvis',options=settings({}))
    plan=[dict(configuration=dict(target='T',element='F467M',setup=None))]
    monkeypatch.setattr(o,'run',lambda c,e,dry_run=False:plan if dry_run else [d])
    monkeypatch.setattr(o,'internal_cross_filter',lambda *a:{'F467M':dict(status='skipped')})
    monkeypatch.setattr(o,'run_uvis',lambda *a,**k:pytest.fail('Absolute stage must be blocked'))
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
    return c


@pytest.mark.parametrize('native_backend',[False,True],ids=['dispatch','native'])
def test_internal_cross_filter_has_no_external_query_and_reuses_receipt(flc,tmp_path,monkeypatch,native_backend):
    import shutil
    from astropy.io import fits
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.astrometry import cross_filter as c,alignment
    opts=settings({})['astrometry'];root=tmp_path/'redux';records={};files=[]
    for filt,offset in [('F467M',0),('F621M',3)]:
        paths=ReductionPaths(root,'TEST','uvis',filt);paths.area('intermediate','drizzle',create=True)
        p=paths.file('intermediate','drizzle','provisional','fits');shutil.copy2(flc(filt,size=256),p)
        with fits.open(p,mode='update') as hd:
            header=hd['SCI'].header;header['CRPIX1']+=offset
            header['CD1_1']=header.pop('CDELT1');header['CD2_2']=header.pop('CDELT2');header['CD1_2']=0.;header['CD2_1']=0.
        files.append(p);records[filt]=dict(filter=filt,setup=None,scientific_assessment=dict(relative_astrometry=dict(status='validated')))
    statefile=paths.directory/'manifest.json';statefile.write_text(json.dumps(dict(filters=records)))
    monkeypatch.setattr(c,'sources',_internal_catalog)
    monkeypatch.setattr(alignment,'query_gaia_dr3',lambda *a,**k:pytest.fail('No external catalog in internal matching'))
    calls=[];real_native=c.native_catalog_fit
    def native(inputs,catalogs,options,reference_image=None,**kwargs):
        calls.append(1)
        if native_backend:return real_native(inputs,catalogs,options,reference_image=reference_image,**kwargs)
        with fits.open(reference_image) as hd:value=hd['SCI'].header['CRPIX1']
        with fits.open(inputs[0],mode='update') as hd:hd['SCI'].header['CRPIX1']=value
        return [dict(shift=[3,0])]
    monkeypatch.setattr(c,'native_catalog_fit',native)
    original=[p.read_bytes() for p in files]
    result=c.internal_cross_filter(paths.directory,list(records),opts)
    assert result['F467M']['status']=='anchor' and result['F621M']['status']=='validated'
    assert [p.read_bytes() for p in files]==original
    result=c.internal_cross_filter(paths.directory,list(records),opts)
    assert len(calls)==1 and result['F621M']['status']=='validated'
    opts['absolute']['search_radius_arcsec']=3
    c.internal_cross_filter(paths.directory,list(records),opts)
    assert len(calls)==1, 'External matching settings must not invalidate internal registration'


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
