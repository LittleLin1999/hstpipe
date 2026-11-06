"""Only persistent-source/CR image QA; no registration or network tests."""
import numpy as np
from hstpipe.uvis.source_qa import aperture_assessment

def scene():
    yy,xx=np.indices((77,77));star=100*np.exp(-((xx-38.2)**2+(yy-37.9)**2)/(2*1.3**2));rng=np.random.default_rng(14)
    stack=np.array([star+rng.normal(0,.05,star.shape) for _ in range(4)]);weights=np.ones_like(stack)*100
    return star,stack,weights

def test_nonpositive_reference_flux_has_missing_ratio_not_nan():
    import json
    star,stack,weights=scene()
    stack=-np.abs(stack);science=-star
    result,_=aperture_assessment(stack,weights,science,np.ones_like(star)*400,(38.2,37.9))
    assert all(row['reference_flux']<0 and row['loss_fraction'] is None and row['status']=='low_snr' for row in result['apertures'])
    json.dumps(result,allow_nan=False)


def test_no_cr_reference_clears_only_4096_on_copies(flc,tmp_path):
    from astropy.io import fits
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.uvis.source_qa import _diagnostic_copy
    source=flc(chips=(1,2));destination=tmp_path/'qa_flc.fits'
    with fits.open(source,mode='update') as hd:
        hd[0].header['HPCRRUN']=True
        for v in (1,2):hd['DQ',v].data[5,5]=4096|8192|256
    before=source.read_bytes()
    _diagnostic_copy(source,destination,ReductionPaths(tmp_path/'redux','T','uvis','F275W'))
    with fits.open(source) as raw,fits.open(destination) as qa:
        for v in (1,2):
            assert qa['DQ',v].data[5,5]==8192|256
            assert np.array_equal(raw['SCI',v].data,qa['SCI',v].data)
        assert not qa[0].header['HPCRRUN'] and qa[0].header['QAONLY']
    assert source.read_bytes()==before

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
    backups=list(paths.directory.parent.glob('.uvis-T-science-before-*/*_drc.fits'))
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


@pytest.fixture
def pending_human_review(flc,tmp_path):
    from astropy.io import fits
    from hstpipe.common.paths import ReductionPaths
    from hstpipe.common.provenance import identity
    from hstpipe.common.target_runner import save
    from hstpipe.uvis.inventory import validate_flc
    raw=flc(size=64);paths=ReductionPaths(tmp_path/'redux','TEST','uvis','F275W')
    with fits.open(raw) as hd:
        header=hd['SCI'].header.copy();header['BUNIT']='ELECTRONS/S';header['ASTSTAT']='VALIDATED';header['ASTREL']='validated';data=hd['SCI'].data.copy()
    review=paths.file('intermediate','drizzle','unvalidated_final','fits');review.parent.mkdir(parents=True)
    fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(data,header=header,name='SCI'),fits.ImageHDU(np.ones_like(data),name='WHT'),fits.ImageHDU(np.ones_like(data,dtype='i4'),name='CTX')]).writeto(review)
    source=dict(source_id='Gfield',status='review_required',output_x=55.,output_y=55.,failures=['significant_persistent_source_flux_loss'],observations=[],apertures=[dict(radius_pix=2.5,status='review_required',loss_fraction=.03,loss_significance=4.)])
    assessment=dict(status='image_quality_review_required',relative_astrometry=dict(status='validated'),relative_after_tweakback=dict(status='validated'),astrometry=dict(status='validated',final_validation=dict(status='validated')),final=dict(warnings=[],validation_status='review_required'),source_validation=dict(status='review_required',version='persistent-source-3',validated_sources=2,failed_sources=1,sources=[source]),unvalidated_final=str(review),effective_parameters=dict(background=dict(target_radius=.3)))
    record=dict(filter='F275W',setup=None,status='image_quality_review_required',products=[],inputs=[identity(raw)],observing_configurations=[validate_flc(raw)],scientific_assessment=assessment)
    state=dict(schema_version=2,instrument='uvis',target='TEST',filters={'F275W':record,'F621M':dict(filter='F621M',status='relative_astrometry_unvalidated',products=[],scientific_assessment={})})
    save(paths.directory/'manifest.json',state)
    return paths,review,raw,state


def test_human_review_publishes_exact_pixels_with_warning_provenance(pending_human_review):
    import json
    from astropy.io import fits
    from hstpipe.common.provenance import identity
    from hstpipe.uvis import approve_source_review
    paths,review,raw,state=pending_human_review;before=identity(review);raw_before=identity(raw)
    result=approve_source_review(paths.directory,'F275W',reviewer='Human reviewer',reason='Accept field-source aperture discrepancy for intended target analysis',source_ids=['Gfield'],expected_sha256=before['sha256'])
    updated=json.loads((paths.directory/'manifest.json').read_text());record=updated['filters']['F275W'];a=record['scientific_assessment']
    assert record['status']=='complete' and a['status']=='science_ready'
    assert a['source_validation']==state['filters']['F275W']['scientific_assessment']['source_validation']
    assert updated['filters']['F621M']==state['filters']['F621M']
    assert identity(review)==before and identity(raw)==raw_before
    with fits.open(result) as hd,fits.open(review) as original:
        assert hd[0].header['HPQASTAT']=='REVIEWED' and hd[0].header['HPRVAUTO']=='REVIEW'
        for name in ('SCI','WHT','CTX'):assert np.array_equal(hd[name].data,original[name].data)
        assert hd['SCI'].header['ASTSTAT']=='VALIDATED'
    product_before=identity(result)
    assert approve_source_review(paths.directory,'F275W',reviewer='Human reviewer',reason='Same approval',source_ids=['Gfield'],expected_sha256=before['sha256'])==result
    assert identity(result)==product_before
    html=paths.report().read_text();assert 'Human review accepted' in html and 'Ready — human reviewed' in html
    assert paths.report().stat().st_size<1_000_000


@pytest.mark.parametrize('problem',['astrometry','checksum','target','morphology'])
def test_human_review_cannot_waive_unrelated_or_stale_failures(pending_human_review,problem):
    from hstpipe.common.provenance import identity
    from hstpipe.common.target_runner import save
    from hstpipe.uvis import approve_source_review
    paths,review,raw,state=pending_human_review;a=state['filters']['F275W']['scientific_assessment'];sha=identity(review)['sha256']
    if problem=='astrometry':a['astrometry']['final_validation']['status']='unvalidated'
    if problem=='checksum':sha='0'*64
    if problem=='target':a['source_validation']['sources'][0].update(output_x=31.,output_y=31.)
    if problem=='morphology':a['source_validation']['sources'][0]['failures']=['source_shaped_weight_deficit']
    save(paths.directory/'manifest.json',state)
    with pytest.raises(ValueError):approve_source_review(paths.directory,'F275W',reviewer='Human reviewer',reason='Accept field-source flux only',source_ids=['Gfield'],expected_sha256=sha)
    assert not paths.product('drc').exists()


@pytest.mark.parametrize('problem',[None,'changed','missing_identity','astrometry'])
def test_user_review_command_keeps_checksum_internal(pending_human_review,problem,capsys):
    from astropy.io import fits
    from hstpipe.common.provenance import identity
    from hstpipe.common.target_runner import save
    from hstpipe.cli import main
    paths,review,raw,state=pending_human_review
    record=state['filters']['F275W'];record['artifacts']=[identity(review)]
    if problem=='missing_identity':record['artifacts']=[]
    if problem=='astrometry':record['scientific_assessment']['astrometry']['final_validation']['status']='unvalidated'
    save(paths.directory/'manifest.json',state)
    if problem=='changed':
        with fits.open(review,mode='update') as hd:hd['SCI'].data[0,0]+=1
    args=['uvis','review','--directory',str(paths.directory),'--filter','F275W']
    if problem in ('changed','missing_identity'):
        with pytest.raises(ValueError):main(args+['--approve','--reviewer','User','--reason','Accept field-source flux discrepancy'])
        assert not paths.product('drc').exists()
        return
    main(args)
    assert 'Preview only' in capsys.readouterr().out and not paths.product('drc').exists()
    if problem=='astrometry':
        with pytest.raises(ValueError):main(args+['--approve','--reviewer','User','--reason','Accept field-source flux discrepancy'])
        assert not paths.product('drc').exists()
    else:
        before=identity(raw)
        main(args+['--approve','--reviewer','User','--reason','Accept field-source aperture discrepancy'])
        assert 'Published:' in capsys.readouterr().out
        assert identity(raw)==before
        assert fits.getheader(paths.product('drc'))['HPQASTAT']=='REVIEWED'
