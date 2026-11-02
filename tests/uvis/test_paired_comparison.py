"""Same measurements on both sides of an astrometric improvement decision."""
import numpy as np
from astropy.table import Table
from hstpipe.uvis.astrometry.relative_engine import paired_comparison
from hstpipe.uvis.config import settings


def case():
 a=np.array([[0.,0.],[100.,0.]])
 xy=[a,a+[-1.,0.],a+[1.,0.]]
 cats=[Table(dict(common_source_id=['A','B'])) for _ in xy]
 return xy,cats,settings({})['astrometry']['relative']


def test_unchanged_wcs_cannot_gain_from_different_exposure_pair_sets():
 xy,cats,options=case();q=paired_comparison(xy,xy,cats,options)
 assert q['exposure_pairs']==3 and q['matched_pair_measurements']==6
 assert np.isclose(q['before_rms_pix'],np.sqrt(2.))
 assert q['improvement_pix']==0.
 # Old comparison used all three pairs before, reference pairs only after.
 assert np.sqrt(2.)-1.>.1


def test_missing_sources_use_identical_available_pairs_before_after():
 xy,cats,options=case();xy[2]=xy[2][:1];cats[2]=cats[2][:1]
 after=[xy[0],xy[1]+[1.,0.],xy[2]+[-1.,0.]]
 q=paired_comparison(xy,after,cats,options)
 assert q['matched_pair_measurements']==4
 assert np.isclose(q['before_rms_pix'],np.sqrt(7/4)) and q['after_rms_pix']==0
 assert [p['matched'] for p in q['pairs']]==[2,1,1]


def test_rejected_observation_is_excluded_equally_before_and_after():
 xy,cats,options=case();fits=[dict(source_ids=['A','B']) for c in cats]
 q=paired_comparison(xy,xy,cats,options,fits_report=fits,inlier_mask=[True,True,True,False,True,True])
 assert q['matched_pair_measurements']==4 and q['improvement_pix']==0
 assert q['pairs'][0]['source_ids']==['A'] and q['pairs'][2]['source_ids']==['A']


def test_report_exposes_same_pair_improvement_and_criterion():
 from hstpipe.uvis.reporting.astrometry import transformations
 q=dict(status='measured',before_rms_pix=.2684,after_rms_pix=.2026,improvement_pix=.0658,exposure_pairs=6,matched_pair_measurements=28,basis='Identical measured source pairs')
 ctx=dict(relative=dict(status='validated',method='existing_wcs'),joint=dict(fits=[],paired_comparison=q,required_improvement_pix=.1))
 # An empty proposal is not a real fit; provide a reference record for rendering.
 ctx['joint']['fits']=[dict(exposure='ref_flc.fits',shift=[0.,0.],geometry='shift')]
 text=transformations(ctx)
 assert '0.2684' in text and '0.2026' in text and '0.0658' in text and '0.1000' in text
 assert 'identical matched sources' in text
