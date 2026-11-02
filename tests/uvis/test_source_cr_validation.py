"""Only persistent-source/CR image QA; no registration or network tests."""
import numpy as np
from hstpipe.uvis.source_qa import aperture_assessment

def scene():
    yy,xx=np.indices((77,77));star=100*np.exp(-((xx-38.2)**2+(yy-37.9)**2)/(2*1.3**2));rng=np.random.default_rng(14)
    stack=np.array([star+rng.normal(0,.05,star.shape) for _ in range(4)]);weights=np.ones_like(stack)*100
    return star,stack,weights

def test_correctly_removed_cosmic_rays_do_not_fail_stellar_qa():
    star,stack,weights=scene()
    for i in range(4):stack[i,32+i*2,35+i*3]+=10000
    result,_=aperture_assessment(stack,weights,star,np.ones_like(star)*400,(38.2,37.9))
    assert result['status']=='validated',result

def test_rejected_real_source_flux_is_not_excused_as_cosmic_ray():
    star,stack,weights=scene();result,_=aperture_assessment(stack,weights,star*.5,np.ones_like(star)*400,(38.2,37.9))
    assert result['status']=='review_required' and 'significant_persistent_source_flux_loss' in result['failures']

def test_source_shaped_weight_hole_is_detected():
    star,stack,weights=scene();yy,xx=np.indices(star.shape);wht=np.ones_like(star)*400;wht[(xx-38.2)**2+(yy-37.9)**2<3**2]=100
    result,_=aperture_assessment(stack,weights,star,wht,(38.2,37.9))
    assert 'source_shaped_weight_deficit' in result['failures']

def test_real_source_broadening_is_detected():
    star,stack,weights=scene();yy,xx=np.indices(star.shape);broad=100*(1.3/2.)**2*np.exp(-((xx-38.2)**2+(yy-37.9)**2)/(2*2.**2))
    result,_=aperture_assessment(stack,weights,broad,np.ones_like(star)*400,(38.2,37.9))
    assert result['status']=='review_required' and 'persistent_source_broadening' in result['failures']

def test_missing_coverage_is_not_a_success():
    star,stack,weights=scene();weights[:,:,35:40]=0
    result,_=aperture_assessment(stack,weights,star,np.ones_like(star)*400,(38.2,37.9))
    assert result['status']=='review_required' and 'insufficient_source_coverage' in result['failures']

def test_partial_background_cutout_does_not_discard_a_covered_star():
    star,stack,weights=scene()
    stack[:,:12]=np.nan;weights[:,:12]=0
    final=star.copy();wht=np.ones_like(star)*400;final[:12]=np.nan;wht[:12]=0
    result,_=aperture_assessment(stack,weights,final,wht,(38.2,37.9))
    assert result['status']=='validated',result

def test_publication_rejects_unvalidated_assessment(tmp_path):
    import pytest
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.source_qa import publish_validated_image
    paths=ReductionPaths(tmp_path,'T','uvis','F275W')
    with pytest.raises(ValueError,match='must pass'):
        publish_validated_image(tmp_path/'does_not_exist.fits',paths,{'status':'review_required'})
    assert not paths.product('drc').exists()

def test_publication_preserves_arrays_review_and_previous_product(flc,tmp_path):
    from astropy.io import fits
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.common.provenance import identity
    from hstpipe.uvis.source_qa import publish_validated_image
    source=flc(size=64)
    with fits.open(source) as hd:
        header=hd['SCI'].header.copy();header['BUNIT']='ELECTRONS/S';header['ASTSTAT']='VALIDATED';header['ASTREL']='validated'
        data=hd['SCI'].data.copy()
    review=tmp_path/'review.fits';fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(data,header=header,name='SCI'),fits.ImageHDU(np.ones_like(data),name='WHT'),fits.ImageHDU(np.ones_like(data,dtype='i4'),name='CTX')]).writeto(review)
    before=identity(review);paths=ReductionPaths(tmp_path/'redux','T','uvis','F275W');paths.copy(review,paths.product('drc'))
    with fits.open(paths.product('drc'),mode='update') as hd:hd['SCI'].data*=.9
    prior=identity(paths.product('drc'))
    product=publish_validated_image(review,paths,dict(status='validated',validated_sources=2))
    with fits.open(product) as hd:
        assert np.array_equal(hd['SCI'].data,data) and hd[0].header['HPQASTAT']=='VALIDATED'
    assert identity(review)==before
    backups=list(paths.directory.parent.glob('.uvis-science-before-*/*_drc.fits'))
    assert len(backups)==1 and identity(backups[0])['sha256']==prior['sha256']

def test_no_input_coverage_is_explicit_and_json_safe():
    import json
    star,stack,weights=scene();weights[:]=0
    result,_=aperture_assessment(stack,weights,star,np.zeros_like(star),(38.2,37.9))
    assert result['status']=='review_required'
    json.dumps(result,allow_nan=False)

import pytest
@pytest.mark.parametrize('source_status',['validated','review_required'])
def test_final_pipeline_uses_persistent_source_gate_not_raw_peaks(flc,tmp_path,monkeypatch,source_status):
    from test_astrometry import preparation,grid
    from hstpipe.uvis.astrometry import workflow as w
    from hstpipe.uvis import source_qa
    from hstpipe.uvis.inventory import discover
    from hstpipe.common.provenance import identity
    from astropy.io import fits
    from pathlib import Path
    p,_,paths,refs,opts=preparation(flc,tmp_path);second=flc('other',size=256);rows=discover(tmp_path);before=[identity(x) for x in (p,second)]
    provisional=paths.file('intermediate','drizzle','provisional','fits');paths.copy(p,provisional)
    solution=dict(status='validated',reference='gaia_dr3',reference_mode='gaia_dr3',observation_epoch=2017.7,exposure_epoch_span_years=0,
        correction_applied=True,corrected_wcs_name='HP-GAIA-DR3',thresholds=dict(opts['absolute'],**opts['catalogs']['gaia_dr3']['validation']))
    monkeypatch.setattr(w,'absolute_alignment',lambda *args:(dict(solution),refs));monkeypatch.setattr(w,'gaia_at_epoch',lambda r,e:r)
    def tb(provisional,relative,working,*args):
        for path in working:
            with fits.open(path,mode='update') as hd:hd['SCI'].header['WCSNAME']='HP-GAIA-DR3'
        h=paths.file('intermediate','alignment','qa_fixture_headerlet','fits');fits.PrimaryHDU().writeto(h,overwrite=True)
        return dict(headerlets=[identity(h)],status='applied')
    monkeypatch.setattr(w,'tweakback_inputs',tb);monkeypatch.setattr(w,'validate_relative_transfer',lambda *args:dict(status='validated'))
    def drizzle(inputs,prefix,options):
        with fits.open(inputs[0]) as hd:
            header=hd['SCI'].header.copy();header['BUNIT']='ELECTRONS/S';shape=hd['SCI'].shape
            fits.HDUList([fits.PrimaryHDU(header=hd[0].header.copy()),fits.ImageHDU(hd['SCI'].data/100,header=header,name='SCI'),fits.ImageHDU(np.ones(shape),name='WHT'),fits.ImageHDU(np.ones(shape,'i4'),name='CTX')]).writeto(str(prefix)+'_drc.fits')
    monkeypatch.setattr(w,'drizzle',drizzle);monkeypatch.setattr(w,'apply_cr_masks',lambda *args:None)
    def cr(original,working,mask,*args,**kwargs):
        fits.PrimaryHDU().writeto(mask,overwrite=True)
        return [dict(extver=1,current_cr_fraction=.03,target_cr_fraction=.05 if source_status=='validated' else 0.,bright_source_cr_fraction=.9 if source_status=='validated' else 0.)]
    monkeypatch.setattr(w,'cr_report',cr);monkeypatch.setattr(w,'final_report',lambda *args,**kwargs:dict(warnings=[],limitations=[]))
    monkeypatch.setattr(w,'detect_sources',lambda image,*args,**kwargs:grid(image))
    monkeypatch.setattr(source_qa,'validate_source_rejection',lambda *args,**kwargs:dict(status=source_status,validated_sources=2))
    from hstpipe.uvis.config import settings
    products,result=w.finalize(provisional,[p,second],rows,paths,opts,settings({})['drizzle'],dict(status='validated'))
    assert bool(products)==(source_status=='validated')
    assert [identity(x) for x in (p,second)]==before
    if products:
        with fits.open(products[0]) as hd:assert hd[0].header['HPQASTAT']=='VALIDATED'
    else:assert result['status']=='image_quality_review_required' and result['unvalidated_final']
