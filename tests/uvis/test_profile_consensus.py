from copy import deepcopy
from hstpipe.uvis.astrometry.stellar import choose_profile,resolve_profile_group

def profile(sid,flux=10.,concentration=.45,scale=1.):
    return dict(source_id=sid,flux=flux*100,flux_error=20.,flux_rate=flux,fwhm_pix=2.4,
        concentration=concentration,concentration_error=.02,aperture_radii_pix=[1.25*scale,5*scale],
        aperture_definition=dict(unit='arcsec',radii=[.05,.2]),status='credible_stellar_fit',
        detection_method='aperture_profile',extver=1,x=50.,y=50.,snr=50.,width_pix=1.1,ellipticity=.1,centroid_error_pix=.1)

def record(name,choices):
    return dict(exposure=name,psf_audit=dict(fits=choices),psf_selected=choices[0],status='valid_detection')

def test_equal_angular_apertures_allow_different_detector_scale():
    a=profile('a');b=profile('b',scale=1.000024);bad=profile('other',flux=100.)
    assert choose_profile([b,bad],a,{}) is b
    b['aperture_definition']['radii']=[.1,.4]
    assert choose_profile([b,bad],a,{}) is None

def test_unique_artifact_is_not_accepted_against_inconsistent_anchor():
    a=profile('a');b=profile('b',flux=100.)
    assert choose_profile([b],a,{}) is None

def test_bad_exposure_does_not_veto_three_repeatable_profiles():
    group=[record('a',[profile('a'),profile('artifact',flux=100.)]),record('b',[profile('b',flux=1000.)]),
           record('c',[profile('c',scale=1.000024)]),record('d',[profile('d')])]
    resolve_profile_group(group,{})
    assert [r.get('psf_selected',{}).get('source_id') for r in group]==['a',None,'c','d']
    assert group[1]['status']=='inconsistent_profile'
    assert all(r['profile_consensus']['supporting_exposures']==3 for r in group)

def test_competing_repeated_objects_are_not_selected_by_brightness():
    group=[record(str(i),[profile(str(i)+'A'),profile(str(i)+'B',flux=100.)]) for i in range(3)]
    resolve_profile_group(group,{})
    assert all(r['status']=='ambiguous_stellar_profiles' and not r.get('psf_selected') for r in group)

def test_ready_report_keeps_single_source_pair_limit_visible():
    from hstpipe.uvis.report_decisions import decision
    record=dict(products=[{'path':'final.fits'}],scientific_assessment=dict(status='science_ready',alignment=dict(status='validated',pairs=[dict(validation_scope='local_shift_only')])) )
    state,why,action=decision(record)
    assert state=='Ready with sparse-source limits' and 'not independently validated' in why

def test_blocked_report_names_exposure_without_common_centroid():
    from hstpipe.uvis.report_decisions import decision
    record=dict(scientific_assessment=dict(alignment=dict(status='unvalidated',original=dict(pairs=[dict(exposure='missing_flc.fits',matched=0)]))))
    state,why,action=decision(record)
    assert state=='Blocked' and 'missing' in why and 'No reliable common centroid' in why
