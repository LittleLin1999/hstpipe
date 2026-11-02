"""Offline scientific guardrails, using synthetic UVIS data and official API mocks."""
import json,sys,types
from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from hstpipe.uvis.inventory import validate_flc,discover
from hstpipe.uvis.grouping import group_exposures
from hstpipe.uvis.config import settings
from hstpipe.uvis.background import estimate,prepare_background
from hstpipe.uvis.alignment import acceptable_fit
from hstpipe.uvis.qa import cr_report,validate_drizzle
from hstpipe.common.provenance import identity


def test_inventory_calibration_geometry(flc,tmp_path):
    a=flc(chips=(7,9),sub=False);m=validate_flc(a)
    assert m['imsets']==[7,9] and m['geometry'][0]['sci']==2
    with fits.open(a,mode='update') as h:h[0].header['PCTECORR']='OMIT'
    assert any('PCTECORR' in s for s in validate_flc(a)['warnings'])
    with fits.open(a,mode='update') as h:h['ERR',7].data=np.ones((5,5))
    with pytest.raises(ValueError,match='Mismatched'):validate_flc(a)

def test_groups_across_geometry_and_filters(flc,tmp_path):
    flc('a',sub=False,chips=(1,2));flc('b');flc('c',filt='F336W');flc('d',ra=11)
    groups=group_exposures(discover(tmp_path))
    assert sorted(map(len,groups.values()))==[1,1,2]
    assert sum(k[2] is not None for k in groups)==2

def test_background_protects_target_and_fallback(flc):
    p=flc(size=256);opts=settings({})['background'];opts.update(box_size=32,target_radius=0.4,min_pixels=100)
    with fits.open(p) as h:
        a=h['SCI',7];model,mask,r=estimate(a.data,h['ERR',7].data,h['DQ',7].data,WCS(a.header),[10,20],opts)
        assert mask[128,128] and r['reliable'] and abs(r['applied']-5)<0.2
        opts['min_blank']=1.;model,mask,r=estimate(a.data,h['ERR',7].data,h['DQ',7].data,WCS(a.header),[10,20],opts)
        assert not model.any() and r['warnings']

def test_sky_once_and_immutability(flc,tmp_path):
    p=flc();before=identity(p);directory=tmp_path/'intermediate/background';directory.mkdir(parents=True)
    dest=tmp_path/'work.fits';opts=settings({})['background'];opts.update(box_size=16,target_radius=0.1,min_pixels=100)
    prepare_background(p,dest,directory,opts)
    assert identity(p)==before and (directory/'test_mask.fits').exists()
    with fits.open(dest) as h:assert h['SCI',7].header['HPSKY']>0 and h[0].header['HPSKYSUB']
    with pytest.raises(ValueError,match='already subtracted'):prepare_background(dest,tmp_path/'again.fits',directory,opts)

def test_alignment_acceptance():
    o=settings({})['alignment']
    assert acceptable_fit(20,.1,.2,1,2,o)
    assert not acceptable_fit(2,.1,.1,1,1,o)
    assert not acceptable_fit(20,1,.1,1,1,o)
    assert not acceptable_fit(20,np.nan,.1,1,1,o)

def test_cr_core_flag_and_detector_bits(flc,tmp_path):
    import shutil
    p=flc();work=tmp_path/'work.fits';shutil.copy2(p,work)
    with fits.open(work,mode='update') as h:h['DQ',7].data[62:68,62:68]|=4096
    report=cr_report(p,work,tmp_path/'mask.fits',tmp_path,.001)
    assert report[0]['warnings'] and report[0]['bright_source_cr_fraction']>0
    with fits.open(tmp_path/'mask.fits') as h:assert h['SCI',7].data.sum()==36

def test_drizzle_configuration(monkeypatch,tmp_path):
    from hstpipe.uvis.drizzle import drizzle
    calls=[];monkeypatch.setitem(sys.modules,'drizzlepac',types.SimpleNamespace(astrodrizzle=types.SimpleNamespace(AstroDrizzle=lambda *a,**kw:calls.append(kw))))
    drizzle(['a','b'],tmp_path/'out',settings({})['drizzle']);o=calls[0]
    assert o['resetbits']==0 and not o['skysub'] and o['skyuser']=='HPSKY'
    assert o['final_wcs'] and o['final_scale']==.04 and o['final_wht_type']=='EXP'
    assert o['build'] and o['context'] and o['driz_cr'] and not o['driz_cr_corr']

def test_mock_full_reduction_and_stage_reuse(flc,tmp_path,monkeypatch):
    from hstpipe.uvis import reduce_uvis
    from hstpipe.common.config import ReductionConfig
    from hstpipe.uvis import pipeline
    raw=tmp_path/'raw';raw.mkdir();p=flc();p.rename(raw/p.name)
    before=identity(raw/p.name);calls=[]
    def fake(inputs,prefix,options):
        calls.append(options)
        with fits.open(inputs[0]) as h:
            header=h['SCI',7].header.copy();header['BUNIT']='ELECTRONS/S'
            hd=[fits.PrimaryHDU(header=h[0].header.copy()),fits.ImageHDU(h['SCI',7].data/100,header=header,name='SCI'),
                fits.ImageHDU(np.full((128,128),100.),name='WHT'),fits.ImageHDU(np.ones((128,128),'i4'),name='CTX')]
            fits.HDUList(hd).writeto(str(prefix)+'_drc.fits',overwrite=True)
    monkeypatch.setattr(pipeline,'drizzle',fake)
    c=ReductionConfig(str(raw),str(tmp_path/'redux'),'uvis',resume=True,options={'alignment':{'mode':'original'},'astrometry':{'method':'relative_only'}})
    [out]=reduce_uvis(c);product=out/'products/TEST_F275W_drc.fits';validate_drizzle(product)
    with fits.open(product) as h:
        assert [x.name for x in h]==['PRIMARY','SCI','WHT','CTX'];assert h['WHT'].header['WHTTYPE']=='EXP';assert not h[0].header['HPCRCLN']
    receipt=(out/'.cache/TEST_F275W/uvis_alignment.json').read_bytes();cr=(out/'.cache/TEST_F275W/uvis_cr.json').read_bytes()
    c.options['drizzle']={'final_scale':.03};reduce_uvis(c)
    assert (out/'.cache/TEST_F275W/uvis_alignment.json').read_bytes()==receipt and (out/'.cache/TEST_F275W/uvis_cr.json').read_bytes()==cr
    assert identity(raw/p.name)==before
    assert json.loads((out/'manifest.json').read_text())['filters']['F275W']['scientific_assessment']['status']=='review_required'
    assert len(calls)==2

def test_amp_detector_geometry():
    from hstpipe.uvis.amplifiers import detector_x
    h=fits.Header(dict(LTV1=-1900.,LTM1_1=1.))
    x=detector_x(h,(128,300));assert x[0]==1901 and (x<=2048).sum()==148

def test_four_real_dithered_exposures(tmp_path):
    """Opt-in local-FLC integration; cached catalogs may be supplied offline."""
    import os,shutil
    directory=os.environ.get('HSTPIPE_UVIS_INTEGRATION_INPUT')
    if not directory:pytest.skip('Set HSTPIPE_UVIS_INTEGRATION_INPUT to four local dithered FLCs')
    pytest.importorskip('drizzlepac')
    files=list(Path(directory).glob('*_flc.fits'))
    assert len(files)==4,'Integration directory must contain exactly four exposures'
    rows=discover(directory);groups=group_exposures(rows);assert len(groups)==1
    before=[identity(p) for p in files]
    from hstpipe.uvis import reduce_uvis
    config=dict(input=str(directory),output=str(tmp_path/'redux'),instrument='uvis',
        options={'astrometry':{'catalogs':{'offline':os.environ.get('HSTPIPE_UVIS_INTEGRATION_OFFLINE')=='1'}}})
    [out]=reduce_uvis(config)
    state=next(iter(json.loads((out/'manifest.json').read_text())['filters'].values()))
    assert list((out/'intermediate/drizzle').glob('*_provisional.fits'))
    products=list((out/'products').glob('*_drc.fits'))
    if products:
        validate_drizzle(products[0])
        assert state['scientific_assessment']['astrometry']['status'] in ('validated','final_validation_failed')
        assert len(list((out/'intermediate/drizzle').glob('*_final_crmask.fits')))==4
    else:
        assert state['status']=='astrometry_unvalidated'
        assert state['scientific_assessment']['astrometry']['status']=='unvalidated'
    assert len(list((out/'intermediate/drizzle').glob('*_crmask.fits'))) in (4,8)
    assert [identity(p) for p in files]==before

def test_none_and_user_sky(flc,tmp_path):
    p=flc();dest=tmp_path/'copy.fits';d=tmp_path/'intermediate/background';d.mkdir(parents=True)
    opts=settings({'background':{'mode':'none'}})['background']
    prepare_background(p,dest,d,opts)
    with fits.open(p) as a,fits.open(dest) as b:assert np.array_equal(a['SCI',7].data,b['SCI',7].data)
    opts.update(mode='user',user={'test':{'7':[2.0,.1]}})
    prepare_background(p,dest,d,opts)
    with fits.open(p) as a,fits.open(dest) as b:assert np.allclose(a['SCI',7].data-2,b['SCI',7].data)

def test_amp_single_subarray_has_no_relative_correction(flc):
    from hstpipe.uvis.amplifiers import equalize
    p=flc();opts=settings({})['background']
    with fits.open(p) as hd:
        corrections,report=equalize(hd,opts)
    assert not corrections and not report['applied'] and 'one amplifier' in report['warnings'][0]

def test_existing_detector_dq_cannot_reset(flc,tmp_path):
    import shutil
    p=flc()
    with fits.open(p,mode='update') as hd:hd['DQ',7].data[0,0]=256
    work=tmp_path/'work.fits';shutil.copy2(p,work)
    with fits.open(work,mode='update') as hd:hd['DQ',7].data[0,0]=0
    with pytest.raises(ValueError,match='DQ bits were reset'):cr_report(p,work,tmp_path/'mask.fits',tmp_path,.02)

@pytest.mark.parametrize('rms,accepted',[(0.1,True),(1.0,False)])
def test_tweakreg_candidate_acceptance(flc,tmp_path,monkeypatch,rms,accepted):
    from hstpipe.uvis import alignment
    import photutils.detection
    from astropy.table import Table
    catalog=Table(dict(xcentroid=np.arange(20.),ycentroid=np.arange(20.),flux=np.ones(20)))
    monkeypatch.setattr(photutils.detection,'DAOStarFinder',lambda *a,**k:lambda *a,**k:catalog)
    inputs=[flc('one'),flc('two')];before=[identity(p) for p in inputs]
    def fake_align(paths,options):
        lines=[]
        for p in paths:
            lines.append(f'{p.name} 0.1 0.1 0 1 {rms} {rms}')
            Path(p.stem+'_catalog_fit.match').write_text('\n'.join(['1 2 3 4']*20))
            with fits.open(p,mode='update') as hd:hd['SCI',7].header['WCSNAME']='CANDIDATE'
        Path('shifts.txt').write_text('\n'.join(lines))
    monkeypatch.setattr(alignment,'align',fake_align)
    d=tmp_path/'alignment';d.mkdir()
    report=alignment.register(inputs,d,tmp_path/'plots',settings({})['alignment'])
    assert (report['status']=='accepted')==accepted
    if not accepted:assert before==[identity(p) for p in inputs]
    else:
        with fits.open(inputs[0]) as hd:assert hd['SCI',7].header['WCSNAME']=='CANDIDATE'

def test_amp_equalization_with_mapped_flat(flc,tmp_path):
    from hstpipe.uvis.amplifiers import equalize
    p=flc(size=300);flat=tmp_path/'flat.fits'
    h=fits.ImageHDU(np.ones((300,300),dtype='f4'),name='SCI');h.header['CCDCHIP']=2;h.header['LTV1']=-1900.
    fits.HDUList([fits.PrimaryHDU(),h]).writeto(flat)
    with fits.open(p,mode='update') as hd:
        hd[0].header['PFLTFILE']=str(flat);hd[0].header['CCDAMP']='CD';a=hd['SCI',7];a.header['LTV1']=-1900.
        a.data[:,:148]+=3
    opts=settings({'background':{'amplifier_correction':True,'target_radius':.1,'box_size':32}})['background']
    with fits.open(p) as hd:
        original=hd['SCI',7].data.copy();models,report=equalize(hd,opts)
        assert report['applied'] and report['flat']['sha256']
        assert np.allclose(hd['SCI',7].data,original-models[7])
        assert np.median(models[7][:,:148])>1 and np.median(models[7][:,148:])<-1

def test_cr_flags_cannot_be_permitted_in_final_bits():
    with pytest.raises(ValueError,match='CR bit'):settings({'drizzle':{'final_bits':4096}})


def test_official_cr_mask_polarity_and_sci_preservation(flc,tmp_path):
    from hstpipe.uvis.drizzle import apply_cr_masks
    p=flc();before=fits.getdata(p,('SCI',7)).copy()
    good=np.ones(before.shape,dtype='u1');good[64,64]=0
    fits.PrimaryHDU(good).writeto(tmp_path/'test_sci7_crmask.fits')
    apply_cr_masks([p])
    with fits.open(p) as hd:
        assert hd['DQ',7].data[64,64]&4096
        assert not hd['DQ',7].data[10,10]&4096
        assert np.array_equal(hd['SCI',7].data,before)
