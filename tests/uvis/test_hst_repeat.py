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
    assert all(p['status']!='validated' for p in report['pairs'])

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
