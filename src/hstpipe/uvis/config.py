"""Validated UVIS stage settings; detector calibration is never rerun."""
from .astrometry.config import DEFAULTS as ASTROMETRY_DEFAULTS
from ..common.calibration import DEFAULTS as CALIBRATION_DEFAULTS
DEFAULTS=dict(calibration=CALIBRATION_DEFAULTS,execution=dict(target_workers=1,filter_workers=1),astrometry=ASTROMETRY_DEFAULTS,background=dict(mode='auto',threshold=2.5,dilate=8,target_radius=3.0,min_blank=0.2,
    min_pixels=1000,significance=3.0,box_size=128,user={},amplifier_correction=False,amplifier_assess=False,dq_ignore_bits=0,source_fills_field=False),
    alignment=dict(mode='auto',minobj=10,max_rms=0.3,max_shift=10.0,fitgeometry='shift',max_sources=1000),
    drizzle=dict(final_scale=0.04,final_pixfrac=0.8,final_kernel='square',final_wht_type='EXP',
    driz_cr_snr='3.5 3.0',driz_cr_scale='1.2 0.7',driz_sep_pixfrac=1.0,driz_sep_scale=None,
    driz_sep_bits=80,final_bits=80),exclude=[],allow_incomplete_calibration=False,
    target_cr_fraction=0.02,keep_drizzle_intermediates=False)

def settings(options):
    import copy
    out=copy.deepcopy(DEFAULTS)
    unknown=set(options)-set(out)
    if unknown:raise ValueError(f'Unknown UVIS options: {sorted(unknown)}')
    for k,v in options.items():
        if k=='astrometry':
            from ..common.target_runner import merge
            out[k]=merge(out[k],v)
        elif isinstance(out[k],dict):out[k].update(v)
        else:out[k]=v
    for section,extra in (('background',set()),('alignment',{'refcat','transformations'})):
        invalid=set(out[section])-set(DEFAULTS[section])-extra
        if invalid:raise ValueError(f'Unknown {section} settings: {invalid}')
    calibration=out['calibration']
    if set(calibration)-set(CALIBRATION_DEFAULTS):raise ValueError('Unknown calibration settings')
    if any(type(calibration[k]) is not bool for k in ('auto_download','offline')):raise ValueError('Calibration switches must be boolean')
    if not isinstance(calibration['timeout'],(int,float)) or not 0<calibration['timeout']<=600:raise ValueError('Calibration timeout must be in (0,600] seconds')
    if type(calibration['retries']) is not int or not 0<=calibration['retries']<=3:raise ValueError('Calibration retries must be an integer from 0 to 3')
    if out['alignment']['minobj']<3:raise ValueError('At least three matches are required')
    if out['background']['box_size']<8 or out['background']['dilate']<1 or out['background']['target_radius']<0:raise ValueError('Invalid sky mask geometry')
    if out['background']['mode'] not in ('auto','constant','2d','none','user'):raise ValueError('Invalid background mode')
    if out['alignment']['mode'] not in ('auto','relative','original','external','user'):raise ValueError('Invalid alignment mode')
    forbidden=set(out['drizzle'])-set(DEFAULTS['drizzle'])-{'final_rot','final_ra','final_dec','final_outnx','final_outny'}
    if forbidden:raise ValueError(f'Unsafe or unsupported AstroDrizzle overrides: {forbidden}')
    if out['drizzle']['final_wht_type'] not in ('EXP','IVM','ERR'):raise ValueError('Invalid weight type')
    if not 0<out['drizzle']['final_pixfrac']<=1:raise ValueError('pixfrac must be in (0,1]')
    from astropy.nddata.bitmask import bitfield_to_boolean_mask
    import numpy as np
    if bitfield_to_boolean_mask(np.array([4096],dtype='u2'),ignore_flags=out['drizzle']['final_bits'],good_mask_value=True)[0]:
        raise ValueError('final_bits must reject the CR bit 4096')
    if any(type(out['execution'].get(k)) is not int or not 1<=out['execution'][k]<=16 for k in ('target_workers','filter_workers')):raise ValueError('Worker counts must be integers from 1 to 16')
    astro=out['astrometry']
    detection=astro['detection']
    if detection.get('method','dao') not in ('dao','compact_stable'):raise ValueError('Unknown absolute source detection method')
    if detection.get('method')=='compact_stable':
        widths=np.asarray(detection['compact_kernel_widths'],float)
        if widths.shape!=(2,) or not np.isfinite(widths).all() or not 0<widths[0]<widths[1]:raise ValueError('Compact kernels must be finite, positive and increasing')
        for name in ('compact_min_snr','compact_max_centroid_change_pix'):
            value=detection[name]
            if not np.isfinite(value) or value<=0:raise ValueError(f'{name} must be positive and finite')
    error_limit=astro['catalogs']['legacy_dr10']['quality'].get('max_position_error_mas')
    if error_limit is not None and (not np.isfinite(error_limit) or error_limit<=0):raise ValueError('Legacy position-error limit must be positive and finite')
    if astro['workflow']!='relative_then_absolute':raise ValueError('Supported UVIS workflow is relative_then_absolute')
    ident=astro['source_identification']
    if ident['use_external_positions_for_relative_fit']:raise ValueError('External source identities may guide associations, never relative-fit coordinates')
    if set(ident['catalogs'])-{'gaia_dr3','legacy_dr10'}:raise ValueError('Unknown source-identification catalog')
    policies=ident.get('catalog_quality',{})
    if not isinstance(policies,dict) or set(policies)-{'gaia_dr3','legacy_dr10'} or any(not isinstance(v,dict) for v in policies.values()):raise ValueError('Source-identification catalog_quality must map supported catalogs to quality policies')
    if ident['search_radius_arcsec']<=0 or ident['max_search_radius_pix']<=0:raise ValueError('Source search windows must be bounded and positive')
    rel=astro['relative']
    if rel.get('external_catalog',False) or not rel.get('enabled',True) or not rel.get('require_validated_solution',True):raise ValueError('FLC relative registration is mandatory and cannot use external catalogs')
    if not astro['absolute'].get('query_after_provisional_drizzle',True):raise ValueError('External queries require a validated provisional image')
    if isinstance(rel['detection'],str):
        if rel['detection']!='full_field':raise ValueError('Only full_field astrometric detection is supported')
        rel['detection']=copy.deepcopy(DEFAULTS['astrometry']['relative']['detection'])
    if not rel['validate']:raise ValueError('Relative WCS validation is mandatory')
    if rel['fit_geometry'] not in ('auto','shift','rshift','rscale'):raise ValueError('Unsupported relative fit geometry')
    if not 1<=rel['max_iterations']<=5:raise ValueError('Relative iterations must be bounded (1..5)')
    if set(rel['methods'])-{'existing_wcs','tweakreg','improved_catalogs','offset_histogram','iterative_mosaic','image_cross_correlation'}:raise ValueError('Unknown relative matching strategy')
    if astro['absolute']['mode']!='auto':raise ValueError('Absolute reference selection supports mode=auto; use relative_only for comparison reductions')
    if astro['catalogs']['timeout']<=0 or astro['catalogs']['retries']<0:raise ValueError('Catalog timeout must be positive and retries nonnegative')
    if astro['method']=='two_pass' and not astro['final_validation']['enabled']:raise ValueError('two_pass requires final astrometric validation')
    if astro['method']=='two_pass' and not astro['absolute']['require_validated_solution']:raise ValueError('Use relative_only explicitly for unvalidated comparison products')
    if astro['method'] not in ('two_pass','relative_only'):raise ValueError('Unknown astrometry method')
    if astro['absolute']['preferred']!='gaia_dr3' or astro['absolute']['fallback']!='legacy_dr10':raise ValueError('Reference priority must be explicit Gaia DR3 then Legacy DR10')
    if astro['catalogs']['gaia_dr3']['provider']!='astroquery.gaia' or astro['catalogs']['legacy_dr10']['provider']!='pyvo':raise ValueError('Only direct Astroquery/pyvo providers are supported')
    if astro['catalogs']['legacy_dr10']['table']!='ls_dr10.tractor':raise ValueError('Legacy table must be ls_dr10.tractor')
    if astro['absolute']['fitgeometry'] not in ('shift','rscale'):raise ValueError('Sparse-field astrometry supports shift or rscale fits')
    if astro['absolute']['min_matches']<(1 if astro['absolute']['fitgeometry']=='shift' else 2):raise ValueError('Absolute fit is mathematically underconstrained')
    return out
