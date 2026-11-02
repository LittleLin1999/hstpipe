"""Reference-assisted *HST-only* counterpart recovery and held-out checks."""
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry.common_sources import build_common,resolve_repeated_counterparts
from hstpipe.common.provenance import identity

def example(flc,count=4,wrong_label=False):
    inputs=[flc(name=chr(97+i)) for i in range(count)];cats=[];records=[]
    for i,p in enumerate(inputs):
        points=[(20.,20.,'PSFPhotometry')]
        if wrong_label:points.append((70.+20*i,70.,'segmentation'))
        elif i<3:points.append((90.,90.,'segmentation'))
        points.extend([(110.-7*i,35.+4*i,'segmentation'),(35.+6*i,110.-5*i,'segmentation')])
        c=Table(dict(x=[x for x,y,m in points],y=[y for x,y,m in points],extver=[7]*len(points),detection_method=[m for x,y,m in points],
            centroid_error_pix=[.1]*len(points),flux=[1000.]*len(points),snr=[30.]*len(points),width_pix=[1.3]*len(points),ellipticity=[.1]*len(points),
            hst_source_id=[p.stem+'_'+str(j) for j in range(len(points))],source_id=[p.stem+'_'+str(j) for j in range(len(points))],exposure_rootname=[p.stem]*len(points)))
        cats.append(c)
        records.append(dict(source_id='Gstar',exposure=p.name,hst_source_id=str(c['hst_source_id'][0]),status='confirmed_common'))
        if wrong_label:records.append(dict(source_id='Lfalse',exposure=p.name,hst_source_id=str(c['hst_source_id'][1]),status='confirmed_common'))
    return inputs,cats,Table(rows=records)

def test_partial_three_exposure_object_survives_failed_full_catalog_fit(flc):
    inputs,cats,external=example(flc);r=settings({})['astrometry']['relative']
    selected,common,report=build_common(inputs,cats,external,r)
    assert [len(c) for c in selected]==[2,2,2,1]
    assert report['HST_guide_used'] and report['HST_repeat_components']>=1
    assert any(p['status']!='validated' for p in report['pairs'])
    assert report['status']=='candidates_established'

def test_repeated_component_resolves_galaxy_without_filling_absent_exposure(flc):
    inputs,cats,external=example(flc);opts=settings({})['astrometry'];selected,common,report=build_common(inputs,cats,external,opts['relative'])
    rows=[dict(source_id='Lgalaxy',exposure=p.name,hst_source_id='',status='ambiguous_association',reason='many candidates',centroid_method='detector',extver=7,
               predicted_x=90.,predicted_y=90.,x=np.nan,y=np.nan,snr=np.nan,width_pix=np.nan,ellipticity=np.nan,centroid_error_pix=np.nan,flux_rate=np.nan) for p in inputs]
    associations=Table(rows=rows);decisions=resolve_repeated_counterparts(inputs,cats,common,associations,opts['source_identification'])
    assert decisions and list(associations['status'])==['confirmed_common']*3+['ambiguous_association']
    assert np.isnan(associations['x'][-1]) and np.allclose(associations['x'][:3],90)
    assert associations['hst_source_id'][0]==cats[0]['hst_source_id'][1]

def test_conflicting_external_labels_do_not_override_HST_position_evidence(flc):
    inputs,cats,external=example(flc,count=2,wrong_label=True)
    selected,common,report=build_common(inputs,cats,external,settings({})['astrometry']['relative'])
    assert [len(c) for c in selected]==[1,1]

def test_additional_object_is_a_held_out_check_not_forced_to_reference(flc):
    from hstpipe.uvis.astrometry.target_fallback import _additional_checks
    inputs,cats,external=example(flc,count=2);cats=[c[:2] for c in cats]
    for c in cats:c['common_source_id']=['Gstar','HSTnucleus']
    cats[1]['x'][1]+=.2
    record=dict(transformations=[dict(exposure=inputs[1].name,shift=[0.,0.],uncertainty_pix=.14)])
    checks=_additional_checks(inputs,cats,'Gstar',record,settings({})['astrometry']['relative'])
    assert len(checks)==1 and checks[0]['held_out_from_primary_fit'] and checks[0]['consistent']
    assert np.isclose(checks[0]['radial_residual_pix'],.2,atol=1e-6)

def test_partial_object_without_anchor_exposure_still_contributes(flc):
    from hstpipe.uvis.astrometry.target_fallback import _additional_checks
    inputs,cats,external=example(flc,count=3);cats=[cats[0][:1],cats[1][:2],cats[2][:2]]
    cats[0]['common_source_id']=['Gstar']
    for c in cats[1:]:c['common_source_id']=['Gstar','HSTnucleus']
    record=dict(transformations=[dict(exposure=p.name,shift=[0.,0.],uncertainty_pix=.14) for p in inputs[1:]])
    checks=_additional_checks(inputs,cats,'Gstar',record,settings({})['astrometry']['relative'])
    assert len(checks)==1 and checks[0]['reference']==inputs[1].name and checks[0]['exposure']==inputs[2].name

@pytest.mark.parametrize('masked',[False,True])
def test_missing_centroid_is_measured_from_pixels_or_explicitly_rejected(flc,masked):
    from hstpipe.uvis.astrometry.source_recovery import recover_missing
    inputs,cats,external=example(flc);options=settings({})['astrometry'];selected,common,_=build_common(inputs,cats,external,options['relative'])
    with fits.open(inputs[-1],mode='update') as hd:
        yy,xx=np.indices(hd['SCI',7].shape);hd['SCI',7].data+=200*np.exp(-((xx-90)**2+(yy-90)**2)/5)
        if masked:hd['DQ',7].data[89:92,89:92]=8192
    before=[identity(p) for p in inputs];records=recover_missing(inputs,cats,common,options['relative']['detection'],options['relative'])
    assert records
    if masked:assert records[0]['status']!='measured_candidate' and len(cats[-1])==3
    else:
        assert records[0]['status']=='measured_candidate' and len(cats[-1])==4
        assert abs(records[0]['x']-90)<.1 and abs(records[0]['y']-90)<.1
    assert [identity(p) for p in inputs]==before


def test_aperture_vetted_guide_keeps_repeated_sources_when_blind_fit_fails(flc):
    inputs,cats,external=example(flc)
    for c in cats:
        c['detection_method']=np.asarray(c['detection_method']).astype('U32')
        c['detection_method'][0]='aperture_profile'
    selected,common,report=build_common(inputs,cats,external,settings({})['astrometry']['relative'])
    assert report['HST_guide_used']
    assert [len(c) for c in selected]==[2,2,2,1]
    assert report['HST_repeat_components']>=1
    assert any(p['status']!='validated' for p in report['pairs'])
    assert report['status']=='candidates_established'


def test_hst_only_repeat_candidates_do_not_require_external_or_psf_labels(flc):
    inputs,cats,external=example(flc)
    for c in cats:c['detection_method']=np.full(len(c),'segmentation')
    selected,common,report=build_common(inputs,cats,external[:0],settings({})['astrometry']['relative'])
    assert report['initial_wcs_repeat_pairs']>0
    assert [len(c) for c in selected]==[2,2,2,1]
    assert not report['HST_guide_used'] and not report['external_positions_used_for_fit']
    assert report['status']=='candidates_established'

@pytest.mark.parametrize('core_bit',[0,4096,8192])
def test_repeated_primary_source_can_be_recovered_with_recorded_cr_policy(flc,core_bit):
    from hstpipe.uvis.astrometry.source_recovery import recover_missing
    inputs,cats,external=example(flc);options=settings({})['astrometry'];opts=options['relative']
    # One object repeats in three exposures; no second guide in the fourth.
    cats=[c[:1] if i<3 else c[:0] for i,c in enumerate(cats)]
    external=external[:3]
    _,common,_=build_common(inputs,cats,external,opts)
    with fits.open(inputs[-1],mode='update') as hd:
        yy,xx=np.indices(hd['SCI',7].shape)
        hd['SCI',7].data+=200*np.exp(-((xx-20)**2+(yy-20)**2)/3)
        hd['SCI',7].data[17:20,23:25]+=1800
        hd['DQ',7].data[17:20,23:25]=4096
        if core_bit:hd['DQ',7].data[19:22,19:22]=core_bit
    before=[identity(p) for p in inputs];records=recover_missing(inputs,cats,common,opts['detection'],opts)
    assert records and records[0]['search_prior']=='existing_WCS_of_repeated_HST_source'
    if core_bit:
        assert records[0]['status']!='measured_candidate'
        assert len(cats[-1])==0
    else:
        r=records[0];assert r['status']=='measured_candidate' and r['centroid_ignore_bits']==0
        assert r['mask_attempts'][0]['status']=='rejected' and r['mask_attempts'][-1]['status']=='measured_candidate'
        assert abs(r['x']-20)<.1 and abs(r['y']-20)<.1
        assert r['core_valid_fraction']>=opts['detection']['min_core_valid_fraction']
        assert r['external_coordinates_used_for_fit'] is False and r['relative_validation_required']
        assert len(cats[-1])==1 and cats[-1]['centroid_ignore_bits'][0]==0
    assert [identity(p) for p in inputs]==before

def test_recovery_keeps_valid_source_with_false_inherited_core_flag(flc):
    from hstpipe.uvis.astrometry.source_recovery import recover_missing
    inputs,cats,external=example(flc);opts=settings({})['astrometry']['relative']
    cats=[c[:1] if i<3 else c[:0] for i,c in enumerate(cats)]
    _,common,_=build_common(inputs,cats,external[:3],opts)
    with fits.open(inputs[-1],mode='update') as hd:
        yy,xx=np.indices(hd['SCI',7].shape);hd['SCI',7].data+=200*np.exp(-((xx-20)**2+(yy-20)**2)/3)
        hd['DQ',7].data[19:22,19:22]=4096
    records=recover_missing(inputs,cats,common,opts['detection'],opts)
    r=records[0];assert r['status']=='measured_candidate' and r['centroid_ignore_bits']==4096
    assert len(r['mask_attempts'])==1 and abs(r['x']-20)<.1
