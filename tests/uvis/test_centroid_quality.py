"""Focused pixel-mask/centroid regressions; no catalog queries or drizzle."""
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from hstpipe.uvis.astrometry.centroids import validate_centroids


def scene(tmp_path,kind='point'):
    y,x=np.indices((81,81));signal=200*np.exp(-((x-40)**2+(y-40)**2)/(2*1.3**2));dq=np.zeros((81,81),'u2')
    seg=(signal>3).astype(int);width=1.3;minor=1.3
    if kind in ('masked','annulus','fragmented'):
        dq[38:43,38:43]=8192;seg[38:43,38:43]=0
        if kind=='fragmented':seg[(x>40)&(y>40)]=0
    if kind=='irregular':
        signal=sum(70*np.exp(-((x-a)**2+(y-b)**2)/2) for a,b in [(36,36),(44,36),(35,44),(45,44)])
        seg=(signal>10).astype(int);width=4.5;minor=4.
    if kind=='compact':signal=100*np.exp(-((x-40)**2+(y-40)**2)/18);seg=(signal>3).astype(int);width=3.;minor=3.
    if kind=='cr':signal[:]=0;signal[40,40]=1e5;seg[:]=0;seg[40,40]=1;width=.1;minor=.1
    if kind=='saturated':dq[40,40]=256
    if kind=='blended':
        signal+=160*np.exp(-((x-43)**2+(y-40)**2)/3);seg[(x>41)&(signal>3)]=2
    p=tmp_path/(kind+'.fits');sci=fits.ImageHDU((5+signal).astype('f4'),name='SCI');sci.header['EXTVER']=7
    bad=fits.ImageHDU(dq,name='DQ');bad.header['EXTVER']=7
    err=fits.ImageHDU(np.ones((81,81),'f4'),name='ERR');err.header['EXTVER']=7
    primary=fits.PrimaryHDU();primary.header['ROOTNAME']='test'
    fits.HDUList([primary,sci,bad,err]).writeto(p)
    table=Table(dict(x=[40.],y=[40.],extver=[7],flux=[float(signal.sum())],snr=[100.],width_pix=[width],minor_width_pix=[minor],centroid_error_pix=[.1],segment_label=[1]))
    return p,table,{7:seg}


@pytest.mark.parametrize('method',['DAOStarFinder','segmentation','SEP','DrizzlePac','custom'])
def test_normal_source_gets_same_mandatory_quality_gate(tmp_path,method):
    p,t,seg=scene(tmp_path);accepted,audit=validate_centroids(p,t,{},method=method,segmentations=seg,return_all=True)
    assert len(accepted)==1 and audit['centroid_valid'][0]
    assert audit['detection_method'][0]==method and audit['local_masked_fraction'][0]==0


@pytest.mark.parametrize('kind,reason',[('masked','masked_core_centroid'),('annulus','ring_halo_unreliable_center'),('fragmented','ring_halo_unreliable_center'),('cr','cosmic_ray_like'),('saturated','saturated_or_invalid_core'),('blended','severe_blending')])
def test_unreliable_centroids_cannot_reach_matching(tmp_path,kind,reason):
    p,t,seg=scene(tmp_path,kind)
    accepted,audit=validate_centroids(p,t,{},segmentations=seg,return_all=True)
    assert not len(accepted) and reason in audit['quality_flags'][0]
    assert audit['rejection_reason'][0] and not audit['centroid_valid'][0]


@pytest.mark.parametrize('kind',['irregular','compact'])
def test_legitimate_compact_galaxies_survive(tmp_path,kind):
    p,t,seg=scene(tmp_path,kind)
    accepted,audit=validate_centroids(p,t,{},segmentations=seg,return_all=True)
    assert len(accepted)==1, audit['quality_flags']
    if kind=='irregular':assert not audit['centroid_in_segment'][0]


@pytest.mark.parametrize('bad',[np.nan,np.inf])
def test_nonfinite_centroids_retained_in_rejection_table(tmp_path,bad):
    p,t,seg=scene(tmp_path);t['x'][0]=bad
    accepted,audit=validate_centroids(p,t,{},return_all=True)
    assert not len(accepted) and 'nonfinite_centroid' in audit['quality_flags'][0]


def test_missing_centroid_is_rejected(tmp_path):
    p,t,seg=scene(tmp_path);t.remove_column('x')
    accepted,audit=validate_centroids(p,t,{},return_all=True)
    assert not len(accepted) and audit['rejection_reason'][0]=='nonfinite_centroid'


def test_imported_invalid_catalog_is_filtered_before_native_matching(tmp_path,monkeypatch):
    from hstpipe.uvis.astrometry.centroids import ensure_catalog
    from hstpipe.uvis.astrometry import relative_engine
    from hstpipe.uvis.config import settings
    from stwcs import wcsutil
    p,t,seg=scene(tmp_path,'masked')
    clean=ensure_catalog(p,t,{})
    assert len(clean)==0 and clean.meta['centroid_quality']['rejected']==1
    seen=[]
    def matcher(a,b,*args,**kw):
        assert len(a)==len(b)==0, 'Invalid centroid reached native matching'
        seen.append(True);return np.array([],int),np.array([],int)
    monkeypatch.setattr(relative_engine,'match',matcher)
    monkeypatch.setattr(relative_engine,'project',lambda c,*args:np.empty((len(c),2)))
    monkeypatch.setattr(wcsutil,'HSTWCS',lambda *args,**kw:None)
    with pytest.raises(ValueError,match='lacks validated source support'):
        relative_engine.native_catalog_fit([p],[t],settings({})['astrometry']['relative'],reference_image=p,reference_catalog=t)
    assert seen


def test_source_quality_change_invalidates_alignment_not_background(flc,tmp_path):
    from hstpipe.common.cache import StageCache
    from hstpipe.uvis.inventory import discover
    from hstpipe.uvis import storage
    from hstpipe.uvis.config import settings
    p=flc();rows=discover(tmp_path);options=settings({});real=StageCache(tmp_path/'cache')
    shadow=tmp_path/'quality_code.py';shadow.write_text('quality-v1')
    class MappedCache:
        def key(self,inputs,parameters):
            return real.key([shadow if q.name=='centroids.py' else q for q in inputs],parameters)
    cache=MappedCache();before=storage.stage_keys(cache,rows,options,[])
    matchfile=tmp_path/'matches.ecsv';matchfile.write_text('old unsafe matches')
    real.record('uvis_alignment',before[1],[matchfile]);assert real.valid('uvis_alignment',before[1])
    shadow.write_text('quality-v2');after=storage.stage_keys(cache,rows,options,[])
    assert before[0]==after[0] and before[1]!=after[1] and not real.valid('uvis_alignment',after[1])


def test_detector_edge_clipped_core_is_not_reported_as_full_coverage(tmp_path):
    p,t,seg=scene(tmp_path)
    y,x=np.indices((81,81))
    with fits.open(p,mode='update') as hd:hd['SCI',7].data[:]=5+200*np.exp(-((x-.4)**2+(y-40)**2)/(2*1.3**2))
    t['x'][0]=.4
    good,audit=validate_centroids(p,t,{},return_all=True)
    assert len(good)==0 and audit['core_valid_fraction'][0]<.75
    assert 'insufficient_core_coverage' in audit['quality_flags'][0]
