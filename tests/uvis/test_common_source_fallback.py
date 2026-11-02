"""Partial-coverage common objects and comparable no-CR source apertures."""
import json
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from hstpipe.common.paths import ReductionPaths
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry import relative_engine as engine
from hstpipe.uvis.astrometry.target_fallback import target_alignment


def catalogs():
    return [Table(dict(common_source_id=ids,centroid_error_pix=[.08]*len(ids),
                       detection_method=['PSFPhotometry']*len(ids)))
            for ids in (['A','B'],['B','C'],['A','C'])]


def test_partial_coverage_sources_reach_standard_joint_fitter(tmp_path,monkeypatch):
    cats=catalogs();seen=[]
    def joint(inputs,catalogs,options,paths):
        seen.extend([list(c['common_source_id']) for c in catalogs])
        return dict(status='unvalidated',may_combine=False,reason='Image check failed')
    monkeypatch.setattr(engine,'joint_registration',joint)
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F621M')
    r,c=target_alignment([tmp_path/f'{i}.fits' for i in range(3)],paths,
                        settings({})['astrometry']['relative'],catalogs=cats)
    assert seen==[['A','B'],['B','C'],['A','C']]
    assert r['available_common_sources']==['A','B','C']
    assert r['source_counts']==[2,2,2] and r['status']=='unvalidated'
    assert r['reason']=='Image check failed' and not r.get('may_combine')


def test_fallback_does_not_repeat_failed_joint_proposal(tmp_path,monkeypatch):
    def forbidden(*a,**kw):raise AssertionError('Repeated the same joint fit/coadds')
    monkeypatch.setattr(engine,'joint_registration',forbidden)
    previous=dict(status='unvalidated',may_combine=False,reason='Failed unchanged QA')
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F621M')
    r,_=target_alignment([tmp_path/f'{i}.fits' for i in range(3)],paths,
                        settings({})['astrometry']['relative'],catalogs=catalogs(),previous_joint=previous)
    assert r['reused_joint_assessment'] and not r.get('may_combine')
    assert r['validation']==previous


def test_no_cr_comparison_encloses_both_original_source_fragments(tmp_path):
    w=WCS(naxis=2);w.wcs.crpix=[64,64];w.wcs.crval=[10,20]
    w.wcs.cdelt=[-.04/3600,.04/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    yy,xx=np.indices((128,128))
    def source(x):return np.exp(-((xx-x)**2+(yy-64)**2)/(2*1.2**2))
    files=[]
    for name,data in [('before',.5*source(64)+.5*source(79)),('after',source(64))]:
        p=tmp_path/(name+'.fits')
        fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(data,header=w.to_header(),name='SCI'),
                      fits.ImageHDU(np.ones(data.shape),name='WHT')]).writeto(p);files.append(p)
    sky=w.all_pix2world([[64,64]],0)
    cat=Table(dict(ra=sky[:,0],dec=sky[:,1],source_id=['star'],width_pix=[1.2],measurement_extent_arcsec=[.6]))
    paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F621M')
    before_cat=cat.copy();after_cat=cat.copy()
    before_cat['measurement_world']=[json.dumps(w.all_pix2world([[64,64],[79,64]],0).tolist())]
    after_cat['measurement_world']=[json.dumps(w.all_pix2world([[64,64]],0).tolist())]
    r=engine.image_quality(*files,before_cat,settings({})['astrometry']['relative'],paths,'split',after_catalog=after_cat,minimum_sources=1)
    assert r['status']=='validated'
    assert .99<r['after']['flux']/r['before']['flux']<1.01
    assert r['after']['fwhm_pix']<r['before']['fwhm_pix']
    assert r['sources_by_phase']['before'][0]['radius_pix']>=22
    assert r['sources_by_phase']['before'][0]['radius_pix']==r['sources_by_phase']['after'][0]['radius_pix']
