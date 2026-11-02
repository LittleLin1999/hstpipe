"""Sparse relative models use mathematical rank, measured errors and HST IDs."""
import numpy as np
from astropy.table import Table
from hstpipe.uvis.config import settings
from hstpipe.uvis.astrometry.relative_fit import quality,select_model,joint_models


def opts():return settings({})['astrometry']['relative']

def rotate(a,angle):
 t=np.deg2rad(angle);return a@np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]]).T


def test_three_source_fixed_scale_rotation():
 a=np.array([[-500.,-200.],[500.,-200.],[0.,500.]])
 b=rotate(a,.05)+[.7,-.3]
 fit=select_model(b,a,np.full(3,.02),opts())
 assert fit['selected']['geometry']=='rshift'
 q=fit['selected'];assert q['robust_inliers']==3 and abs(q['rotation_deg']-.05)<1e-6
 assert q['scale']==1 and q['loo_rms_pix']<1e-5 and q['rotation_uncertainty_deg']>0


def test_two_source_rotation_is_fit_with_limited_independent_validation():
 a=np.array([[-500.,0.],[500.,0.]])
 q=quality(rotate(a,.05)+[1.,-.5],a,np.full(2,.02),opts(),'rshift')
 assert 'matrix' in q and abs(q['rotation_deg']-.05)<1e-6
 assert q['loo_rms_pix'] is None and q['validation_scope']=='limited_redundancy'
 assert q['rotation_uncertainty_deg']>0


def test_single_source_translation_retains_calibrated_rotation():
 a=np.array([[30.,80.]])
 q=select_model(a+[.2,-.1],a,np.array([.08]),opts())['selected']
 assert q['geometry']=='shift' and np.allclose(q['shift'],[.2,-.1])
 assert not q['rotation_measured'] and q['rotation_uncertainty_deg'] is None
 assert q['validation_scope']=='local_shift_only' and q['transform_uncertainty_pix']>0


def test_clustered_sources_do_not_support_rotation():
 a=np.array([[0.,0.],[.2,0.],[0.,.2]])
 fit=select_model(rotate(a,.05),a,np.full(3,.2),opts())
 assert fit['selected']['geometry']=='shift'
 assert not fit['candidates']['rshift']['rotation_constrained']


def test_joint_missing_measurements_and_HST_only_frame():
 points=np.array([[-500.,0.],[500.,0.],[0.,600.]])
 xy=[points,points+np.array([.4,-.2]),points[[0]]+np.array([-.2,.1])]
 cats=[]
 for i,a in enumerate(xy):
  c=Table(dict(common_source_id=['A','B','C'][:len(a)],centroid_error_pix=[.03]*len(a)))
  cats.append(c)
 r=joint_models(xy,cats,opts())
 assert r['status']=='fit_supported' and r['external_coordinates_used_for_fit'] is False
 assert np.allclose(r['fits'][1]['shift'],[-.4,.2],atol=1e-5)
 assert r['fits'][2]['geometry']=='shift' and not r['fits'][2]['rotation_measured']
 assert r['fits'][2]['validation_scope']=='local_shift_only'
 assert r['distinct_sources']==3 and r['reference_index']==0


def test_duplicate_or_nonfinite_points_cannot_support_rotation():
 a=np.array([[1.,2.],[1.,2.]])
 assert quality(a,a,np.ones(2)*.1,opts(),'rshift')['status']=='degenerate_distribution'
 a[0,0]=np.nan
 assert quality(a,a,np.ones(2)*.1,opts(),'rshift')['status']=='nonfinite'


def test_three_source_native_rotation_updates_all_sci_without_changing_data(flc):
 from astropy.io import fits
 from hstpipe.uvis.astrometry.relative_engine import apply_native_fit
 from hstpipe.uvis.astrometry.wcs import detector_wcs
 p=flc(chips=(1,2),size=128,sub=False)
 # Real calibrated HST WCS use CD matrices; the generic fixture uses PC/CDELT.
 with fits.open(p,mode='update') as hd:
  for v in (1,2):
   h=hd['SCI',v].header;h['CD1_1']=h['CDELT1'];h['CD2_2']=h['CDELT2'];h['CD1_2']=0.;h['CD2_1']=0.
   for key in list(h):
    if key.startswith(('CDELT','PC1_','PC2_')):del h[key]
 with fits.open(p) as hd:
  ref=detector_wcs(hd,1);before=[hd['SCI',v].data.copy() for v in (1,2)]
  positions=np.array([[20.,20.],[100.,20.],[60.,100.]])
  original=[detector_wcs(hd,v).all_pix2world(positions,0) for v in (1,2)]
 f=quality(rotate(positions-60,.05)+60+[.3,-.2],positions,np.full(3,.02),opts(),'rshift')
 apply_native_fit(p,ref,f)
 with fits.open(p) as hd:
  for v,data,old in zip((1,2),before,original):
   assert np.array_equal(hd['SCI',v].data,data)
   assert hd['SCI',v].header['WCSNAME']=='HP-RELATIVE'
   world=detector_wcs(hd,v).all_pix2world(positions,0)
   mapped=ref.wcs_world2pix(world,0)
   expected=(positions-f['center'])@np.asarray(f['matrix']).T+f['center']+f['shift']
   assert np.allclose(mapped,expected,atol=.01)


def test_two_oblique_unequal_weight_sources_use_proper_rotation():
 for angle in (.005,.025,.05):
  a=np.array([[512.07,512.66],[888.77,1009.05]])
  b=rotate(a,angle)+[.03,-.02]
  q=quality(b,a,np.array([.23,.08]),opts(),'rshift')
  assert np.linalg.det(q['matrix'])>0
  assert abs(q['rotation_deg']-angle)<1e-6
  assert q['postfit_scatter_pix']<1e-7


def test_joint_adjustment_can_select_fixed_scale_rotation():
 a=np.array([[-500.,-200.],[500.,-200.],[0.,500.]])
 b=rotate(a,-.05)+[.4,-.2]
 cats=[Table(dict(common_source_id=['A','B','C'],centroid_error_pix=[.02]*3)) for _ in range(2)]
 r=joint_models([a,b],cats,opts())
 assert r['status']=='fit_supported'
 assert r['fits'][1]['geometry']=='rshift'
 assert abs(r['fits'][1]['rotation_deg']-.05)<1e-6
 assert r['fits'][1]['scale']==1 and r['residual_rms_pix']<1e-6


def test_uniform_scale_is_considered_only_when_significantly_measured():
 a=np.array([[-500.,-200.],[500.,-200.],[0.,500.]])
 model=select_model(a*1.0008+[.2,-.1],a,np.full(3,.01),opts())
 assert model['selected']['geometry']=='rscale'
 assert abs(model['selected']['scale']-1.0008)<1e-9
 assert model['selected']['scale_uncertainty']>0
