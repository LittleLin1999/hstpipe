"""Compact offline reduction reports, derived from persisted target manifests."""
from __future__ import annotations
import base64,html,io,json,os
from pathlib import Path
from ..common.paths import ReductionPaths
from ..common.target_runner import save

LIMIT=1_000_000
STAGES=('source_detection','source_identification','relative_alignment','target_specific_alignment',
        'provisional_drizzle','cross_filter_alignment','absolute_alignment','tweakback','final_drizzle','final_validation')

def _escape(value):return html.escape(str(value),quote=True)

def _brief(value,limit=600):
    text=json.dumps(value,ensure_ascii=False,default=str) if isinstance(value,(dict,list)) else str(value)
    return _escape(text if len(text)<=limit else text[:limit]+'…')

def _badge(status):
    if status in ('validated','validated_existing','validated_refined','complete','applied','not_needed_validated_existing','identified'):return 'PASS'
    if status in ('blocked',):return 'BLOCKED'
    if status in ('not_run','not_attempted','not_applicable','disabled','skipped','pending','anchor'):return 'NOT_RUN'
    if status in ('single_exposure','local_only','local_shift_only','coverage_not_verified'):return 'WARNING'
    return 'FAIL'

def _chips(record):
    a=record.get('scientific_assessment',{});relative=a.get('relative_astrometry',a.get('alignment',{}))
    detected={r['exposure']:r.get('chips',[]) for r in a.get('source_identification',{}).get('full_field',{}).get('exposures',[])}
    pairs=relative.get('overlap_graph',{}).get('edges',[])
    rows=[]
    for exposure in a.get('multi_chip',[]):
        name=exposure['exposure'];root=exposure['rootname'];original=root+'_flc.fits'
        for chip in exposure['chips']:
            v=chip['extver'];counts=next((c for c in detected.get(original,[]) if c['extver']==v),{})
            matches=[]
            for pair in pairs:
                side='reference' if pair.get('reference')==name else ('exposure' if pair.get('exposure')==name else None)
                if side:
                    matches.extend(c['robust_inliers'] for c in pair.get('chip_matches',[]) if c['side']==side and c['extver']==v)
            tb=[c for e in a.get('tweakback',{}).get('chip_wcs_after',[]) if e['exposure']==name for c in e['chips'] if c['extver']==v]
            provisional=next((c['status'] for c in a.get('provisional_chips',[]) if c['exposure']==name and c['extver']==v),'NOT_RUN')
            values=(root,str(v),str(chip.get('ccdchip')),str(counts.get('detected_method_rows','unknown')),str(counts.get('valid_unique',chip.get('relative_catalog_sources','unknown'))),
                str(max(matches,default=0))+' candidate inliers/pair',provisional,'updated '+str(tb[0].get('wcsname')) if tb else 'NOT_RUN', 'excluded' if chip.get('excluded') else 'included')
            rows.append('<tr>'+''.join('<td>'+_escape(x)+'</td>' for x in values)+'</tr>')
    if not rows:
        # Validation can fail before detection; discovered imsets still remain explicit.
        for inp in record.get('observing_configurations',[]):
            for c in inp.get('geometry',[]):
                rows.append('<tr><td>'+_escape(inp.get('ROOTNAME',''))+'</td><td>'+str(c['extver'])+'</td><td>'+_escape(c.get('CCDCHIP'))+'</td><td colspan="6">Discovered; later processing NOT_RUN</td></tr>')
    return '<h3>CCD processing</h3><p>Detection rows overlap across methods. QC-passing detections and DQ-recovered bright candidates are not confirmed stars. Match counts are candidate mathematical inliers, not accepted registration. A normal single-SCI subarray requires no missing-chip warning.</p><table><tr>'+''.join('<th>'+h+'</th>' for h in ('Exposure','SCI','CCDCHIP','Detected rows','QC-passing unique','Matches','Provisional coverage','TweakBack','Selection'))+'</tr>'+''.join(rows)+'</table>'

def _absolute_summary(record):
    """Separate catalog eligibility, counterpart selection and WCS accuracy."""
    ast=record.get('scientific_assessment',{}).get('astrometry',{})
    descriptions={
        'insufficient_verified_references':'Insufficient verified references',
        'insufficient_reference_sources':'Insufficient references under the previous software count policy',
        'software_match_count_restriction':'Software match-count restriction',
        'diagnostic_single_reference':'One verified reference: diagnostic shift only; absolute accuracy is not independently validated',
        'ambiguous_associations':'Incorrect or ambiguous associations',
        'excessive_postfit_residuals':'Excessive post-fit residuals or unstable leave-one-out prediction',
        'excessive_postfit_residuals_or_distribution':'Excessive corrected residuals or inadequate source distribution',
        'query_failed':'External catalog unavailable',
        'external_catalog_unavailable':'External catalog unavailable',
        'validated':'Validated absolute solution'}
    rows=[];reference_details=[]
    def number(v):
        if v is None:return '—'
        if isinstance(v,(list,tuple)):return ', '.join(number(x) for x in v)
        return f'{v:.4f}' if isinstance(v,float) else str(v)
    attempts=list(ast.get('attempts',[]))
    compared=record.get('scientific_assessment',{}).get('absolute_reference_comparison',{})
    used={a.get('reference') for a in attempts}
    attempts.extend(dict(a,comparison_only=True) for a in compared.get('attempts',[]) if a.get('reference') not in used)
    for attempt in attempts:
        initial=attempt.get('initial_validation',{});final=attempt.get('validation',{})
        decision=descriptions.get(attempt.get('status'),attempt.get('status'))
        if attempt.get('rejected_transformation'):decision+='; original WCS retained; proposed correction rejected because it worsened the residuals'
        values=(str(attempt.get('reference'))+(' (comparison only)' if attempt.get('comparison_only') else ''),decision,
            attempt.get('eligible_catalog_objects',attempt.get('reference_sources')),attempt.get('covered_references'),
            initial.get('initial_associations'),initial.get('fit_inliers',initial.get('matched')),attempt.get('fit_source_count','not fitted'),initial.get('rejected'),
            initial.get('fitted_shift_arcsec',initial.get('offset_arcsec')),initial.get('rms_arcsec'),initial.get('postfit_rms_arcsec'),
            final.get('rms_arcsec'),initial.get('leave_one_out_max_arcsec'),attempt.get('independent_validation_sources',0))
        rows.append('<tr>'+''.join('<td>'+_escape(number(v))+'</td>' for v in values)+'</tr>')
        if attempt.get('reference')=='gaia_dr3' and attempt.get('reference_audit'):
            try:
                from astropy.table import Table
                audit=Table.read(attempt['reference_audit'],format='ascii.ecsv');details=[]
                for obj in audit[:20]:
                    cells=(obj['source_id'],'covered' if obj['covered'] else 'outside image coverage',obj['association_status'],
                        str(obj['pm_applied']) if 'pm_applied' in audit.colnames else 'not supplied',
                        number(float(obj['position_epoch'])) if 'position_epoch' in audit.colnames else 'not supplied')
                    details.append('<tr>'+''.join('<td>'+_escape(c)+'</td>' for c in cells)+'</tr>')
                reference_details.append('<details><summary>Gaia source identities and epoch checks</summary><table><tr><th>Gaia ID</th><th>Coverage</th><th>HST association</th><th>Proper motion applied</th><th>Position epoch</th></tr>'+''.join(details)+'</table></details>')
            except (OSError,ValueError,KeyError):pass
    headers=('Reference','Decision','Eligible','Covered','Associations','Consensus inliers','Native fit used','Rejected','Proposed east/north shift (arcsec)',
        'Initial WCS RMS (arcsec)','Post-fit scatter (arcsec)','Corrected WCS RMS (arcsec)','LOO maximum (arcsec)','Independent checks')
    return '<h3>Absolute-reference assessment</h3><p>Accepted WCS: '+_escape(ast.get('corrected_wcs_name','none'))+'; archived provisional WCS: '+_escape(ast.get('original_wcs_name','none'))+' (key '+_escape(ast.get('original_wcs_key','none'))+').</p>' + ''.join(reference_details) + '<div class="scroll"><table><tr>'+''.join('<th>'+h+'</th>' for h in headers)+'</tr>'+''.join(rows)+'</table></div><p>Eligible cone objects can lie outside detector coverage. Proposed shifts are not applied until validation passes. Leave-one-out uses alternate fitted subsets and is not a separate withheld reference. A single external reference cannot certify absolute accuracy. Gaia and Legacy retain their catalog-specific accuracy limits; the 0.20-arcsec goal is not relaxed.</p>'


def _source_validation_summary(record):
    q=record.get('scientific_assessment',{}).get('source_validation',{})
    if not q:return ''
    rows=[]
    for source in q.get('sources',[]):
        aperture=next((a for a in source.get('apertures',[]) if a.get('radius_pix')==5),{})
        before=source.get('before') or {};after=source.get('after') or {}
        def num(x):return '—' if x is None else f'{x:.4f}'
        values=(source.get('source_id'),len(source.get('observations',[])),source.get('status'),
            num(1-aperture['loss_fraction']) if aperture.get('loss_fraction') is not None else '—',num(aperture.get('loss_significance')),
            num(before.get('fwhm_pix'))+' / '+num(after.get('fwhm_pix')),num(source.get('core_weight_ratio')))
        rows.append('<tr>'+''.join('<td>'+_escape(v)+'</td>' for v in values)+'</tr>')
    return '<h3>Persistent-source and CR validation: '+_escape(q.get('status'))+'</h3><p>Independent no-CR input drizzles share the exact final grid. Aperture flux, morphology and WHT are compared using empirical input scatter. Raw high pixels and background-area CR fractions mix real transients with astronomical emission and do not establish lost stellar flux. These checks do not relax the CR parameters or manufacture source pixels.</p><div class="scroll"><table><tr>'+''.join('<th>'+x+'</th>' for x in ('HST source','Measured exposures','Decision','Final / reference flux (5-pixel radius)','Deficit significance','FWHM reference / final (pixels)','Core / local WHT'))+'</tr>'+''.join(rows)+'</table></div>'


def _summary(key,record,workflow_error=None,astrometry_context=None):
    if workflow_error:
        # A pre-stage/migration failure must not present old products as new.
        record=dict(record,status='blocked_workflow',error=workflow_error,
                    previous_products=record.get('products',[])+record.get('previous_products',[]),
                    products=[],scientific_assessment={})
    a=record.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}))
    local=rel.get('local_registration',{});stages=a.get('stages',{})
    calibration=all(e.get(k)=='COMPLETE' for e in record.get('observing_configurations',[]) for k in ('BIASCORR','BLEVCORR','DARKCORR','FLATCORR','PCTECORR'))
    rows=['<tr><td>Input calibration / SCI discovery</td><td>'+('PASS' if calibration and record.get('inputs') and not record.get('input_validation_errors') else 'FAIL')+'</td><td>'+_brief(record.get('input_validation_errors',[]),1500)+'</td></tr>']
    sky=[c for chips in a.get('background',{}).values() for c in chips]
    sky_warnings=list(dict.fromkeys(str(w) for c in sky for w in c.get('warnings',[])))
    sky_status='NOT_RUN' if not sky else ('WARNING' if sky_warnings or any(not c.get('reliable',False) for c in sky) else 'PASS')
    rows.append('<tr><td>Background assessment</td><td>'+sky_status+'</td><td>'+_brief(sky_warnings,600)+'</td></tr>')
    for name in STAGES:
        value=stages.get(name,{})
        if name=='target_specific_alignment':value=local or {'status':'not_attempted'}
        if name=='relative_alignment' and not value:value=rel
        if name=='absolute_alignment' and not value:value=a.get('astrometry',{})
        if name=='tweakback' and not value:value=a.get('tweakback',{})
        if name=='final_drizzle' and not value:value={'status':'complete' if record.get('products') else 'not_produced'}
        status=value.get('status','not_run')
        if name=='cross_filter_alignment' and status=='anchor':value=dict(value,reason='Selected reference filter; no cross-filter correction is applied to the anchor')
        badge='WARNING' if name=='cross_filter_alignment' and status=='unvalidated' and a.get('astrometry',{}).get('status')=='validated' else _badge(status)
        rows.append('<tr><td>'+name.replace('_',' ')+'</td><td>'+badge+' — '+_brief(status)+'</td><td>'+_brief(value.get('reason',value.get('error','')),260)+'</td></tr>')
    counts=a.get('source_identification',{}).get('full_field',{}).get('exposures',[])
    inventory=record.get('observing_configurations',[])
    metrics={'exposures':len(record.get('inputs',[])),
             'input_calibration_warnings':[w for e in inventory for w in e.get('warnings',[])][:16],
             'background':{n:[{k:r.get(k) for k in ('extver','measured','applied','uncertainty','reliable','warnings')} for r in rs] for n,rs in a.get('background',{}).items()},
             'residual_sky_range_electrons_per_second':a.get('background_discontinuity_range'),
             'amplifier_assessment':{n:[r['amplifiers'] for r in rs if r.get('amplifiers')] for n,rs in a.get('background',{}).items() if any(r.get('amplifiers') for r in rs)},
             'passing_HST_detections':{e['exposure']:e.get('unique_accepted') for e in counts},
             'relative_method':rel.get('method'),'relative_status':rel.get('status'),
             'local_reference_result':rel.get('local_registration'),
             'HST_repeat_recovery':a.get('source_identification',{}).get('HST_repeat_recovery'),
             'stellar_association':[dict(exposure=s.get('exposure'),source_id=s.get('source_id'),status=s.get('status'),candidates=s.get('raw_candidates'),duplicates=s.get('duplicates'),rejected=s.get('rejected'),selected=s.get('selected'),reason=s.get('reason')) for s in a.get('source_identification',{}).get('stellar_resolution',[])],
             'original_pairs':[{k:r.get(k) for k in ('exposure','status','matched','robust_inliers','rejected','current_rms_pix','postfit_scatter_pix','loo_rms_pix','transform_uncertainty_pix')} for r in rel.get('original',{}).get('pairs',[])],
             'final_absolute_validation':a.get('astrometry',{}).get('final_validation',{}),
             'current_DRC_products':[Path(p['path']).name for p in record.get('products',[])],
             'previous_products_preserved':len(record.get('previous_products',[]))}
    joint=rel.get('joint_models',{})
    metrics['joint_relative_models']={k:joint.get(k) for k in ('status','scope','distinct_sources','measurement_count','model_selection_reason','residual_rms_pix','degrees_of_freedom','corrections_applied','reason')}
    metrics['joint_relative_models']['exposures']=[{k:f.get(k) for k in ('exposure','matched','robust_inliers','rejected','geometry','shift','dx_uncertainty_pix','dy_uncertainty_pix','rotation_measured','rotation_deg','rotation_uncertainty_deg','rotation_fixed_reason','postfit_scatter_pix','robust_scatter_pix','max_correction_uncertainty_over_footprint_pix','validation_scope')} for f in joint.get('fits',[])]
    metrics['joint_relative_models']['candidate_models']={i:{g:{k:f.get(k) for k in ('status','geometry','rotation_deg','rotation_uncertainty_deg','rotation_constrained','postfit_scatter_pix','loo_rms_pix','scale','scale_uncertainty')} for g,f in comparison.get('candidates',{}).items()} for i,comparison in joint.get('model_comparisons',{}).items()}
    metrics['joint_relative_models']['model_rejections']={i:comparison.get('reasons') for i,comparison in joint.get('model_comparisons',{}).items()}
    metrics['dq_policy']=a.get('dq_policy')
    metrics['mdriztab']=a.get('mdriztab')
    metrics['DQ_by_chip']={name:[dict(extver=c.get('extver'),astrometric_mask_fraction=c.get('effective_astrometric_mask_fraction'),bits={bit:round(v['fraction'],5) for bit,v in c.get('dq_bits',{}).items()}) for c in chips] for name,chips in a.get('exposures',{}).items()}
    metrics['recovered_bright_sources']={e['exposure']:{str(c['extver']):c.get('recovered_valid_bright_sources') for c in e.get('chips',[])} for e in counts}
    metrics['CR_reassessment']={name:[{k:r.get(k) for k in ('extver','inherited_4096_fraction','cleared_inherited_4096_fraction','new_cr_fraction','calwf3_8192_fraction','raw_bright_pixel_cr_fraction','raw_bright_pixel_metric_caveat','persistent_source_validation','warnings')} for r in chips] for name,chips in a.get('cr',{}).items()}
    metrics['drizzle']=a.get('effective_parameters',record.get('configuration',{}).get('options',{})).get('drizzle')
    metrics['final_image_QA']={k:a.get('final',{}).get(k) for k in ('coverage_fraction','target_core_weight_ratio','adjacent_pixel_correlation','warnings')}
    metrics['final_source_masked_sky']={k:a.get('final_background_assessment',{}).get(k) for k in ('measured','units','pixel_rms','blank_fraction','source_masked_tile_scatter','structure_peak_to_peak','original_science_modified')}
    notes=[]
    if record.get('error'):notes.append('Processing error: '+record['error'])
    if rel.get('status') not in ('validated','single_exposure'):notes.append('Relative registration is not validated: absolute propagation and new science-ready DRC are blocked. Inspect common-source identities, masked cores and no-CR diagnostics.')
    if local.get('status')=='local_only':notes.append('Target-local translation only; full-field accuracy, rotation and scale are not independently established.')
    if rel.get('status')=='single_exposure':notes.append('One exposure: relative registration cannot be measured; reliable multi-exposure CR cleaning is unavailable.')
    notes.extend(rel.get('warnings',[]))
    if joint:
        notes.append('Relative solution scope: '+str(joint.get('scope','unvalidated'))+'. Fitted-source residuals are not independent full-field validation.')
        local_exposures=[f.get('exposure') for f in joint.get('fits',[]) if f.get('validation_scope')=='local_shift_only']
        if local_exposures:notes.append('Translation only, with calibrated rotation retained and not measured: '+', '.join(local_exposures))
    notes.extend(a.get('final',{}).get('limitations',[]))
    for q in a.get('source_identification',{}).get('queries',[]):
        if q.get('error'):notes.append(q.get('catalog','Catalog')+': '+q['error'])
    notes=list(dict.fromkeys(map(str,notes)))[:24]
    attempts=rel.get('attempts',[])
    alternatives=[{'method':x.get('method'),'status':x.get('status',x.get('assessment',{}).get('status')),'error':x.get('error')} for x in attempts]
    relative_label=rel.get('status','not_run')
    absolute_label=stages.get('absolute_alignment',{}).get('status',a.get('astrometry',{}).get('status','not_run'))
    content='<section id="filter-'+_escape(key)+'"><h2>'+_escape(key)+' — relative '+_escape(relative_label)+'; absolute '+_escape(absolute_label)+'</h2>'
    content+='<p>Processing state: <code>'+_escape(record.get('status','unknown'))+'</code>. Last reduction run: '+_escape(record.get('started_utc','not recorded'))+'.</p>'
    content+='<table><thead><tr><th>Stage</th><th>Status</th><th>Reason</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table>'
    if astrometry_context is not None:
        from .report_astrometry import render
        content+=render(astrometry_context,record)
    content+=_spatial_summary(record)
    content+=_absolute_summary(record)
    content+=_source_validation_summary(record)
    content+=_chips(record)
    content+='<details><summary>Additional QA measurements</summary><pre>'+_brief(metrics,18000)+'</pre></details>'
    content+='<details><summary>Relative alternatives</summary><pre>'+_brief(alternatives,5000)+'</pre></details>'
    content+='<h3>Warnings and actions</h3><ul>'+''.join('<li>'+_escape(n[:1500])+'</li>' for n in notes)+'</ul>'
    if not notes:content+='<p>Inspect final morphology, CR masks, background and coverage before science use.</p>'
    return content+'</section>'

def _figures(paths,filters,budget=650_000,descriptions=None):
    """Bounded JPEG thumbnails: prioritize alignment failures, then final QA."""
    from PIL import Image,ImageOps
    files=list(paths.area('diagnostics').glob('*.png'))+list(paths.area('diagnostics').glob('*.jpg'))
    def rank(p):
        text=p.name.lower()
        if 'centroid_recovery' in text:return (-4,0,p.name)
        return (-3 if 'source_validation' in text else (-2 if 'absolute_references' in text else (-1 if 'stellar_association' in text else (0 if any(x in text for x in ('identity','rejected','coadd','relative_sources')) else 1))),
                0 if any(x in text for x in ('target','weight','coverage')) else 1,p.name)
    chosen=[];used=0;per_filter={k:0 for k in filters}
    for p in sorted(files,key=rank):
        group=next((k for k in filters if '_'+k+'_' in p.name),None)
        if group and per_filter[group]>=2:continue
        try:
            with Image.open(p) as original:
                image=ImageOps.exif_transpose(original).convert('RGB');image.thumbnail((900,620))
                stream=io.BytesIO();image.save(stream,format='JPEG',quality=65,optimize=True)
                data=base64.b64encode(stream.getvalue()).decode('ascii')
        except (OSError,ValueError):continue
        caption=next((desc for (filt,sid),desc in (descriptions or {}).items() if filt==group and sid in p.name),p.name)
        fragment='<figure><img alt="'+_escape(p.stem)+'" src="data:image/jpeg;base64,'+data+'"><figcaption>'+_escape(caption)+'</figcaption></figure>'
        size=len(fragment.encode('utf-8'))
        if used+size>budget:continue
        chosen.append(fragment);used+=size
        if group:per_filter[group]+=1
        if len(chosen)>=3:break
    return chosen

def generate_report(directory):
    directory=Path(directory);statefile=directory/'manifest.json'
    state=json.loads(statefile.read_text())
    paths=ReductionPaths(directory.parent.parent,state['target'],'uvis','report')
    filters=state.get('filters',{})
    from .report_astrometry import context,table,_kind
    contexts={}
    for key,record in filters.items():
        group_paths=ReductionPaths(directory.parent.parent,state['target'],'uvis',record.get('filter',key),record.get('setup'))
        if not state.get('workflow_error'):contexts[key]=context(group_paths,record)
    ready=bool(filters) and all(r.get('products') and r.get('scientific_assessment',{}).get('status')=='science_ready' for r in filters.values())
    title=state['target']+' UVIS reduction'
    from .report_decisions import ready as filter_ready,overview_row,concise_filter,decision
    count=sum(filter_ready(r) for r in filters.values())
    css="""body{font:16px -apple-system,BlinkMacSystemFont,Segoe UI,system-ui,sans-serif;max-width:1050px;margin:36px auto;padding:0 20px;color:oklch(25% .015 250);background:oklch(99% .005 250)}h1{font-size:1.7rem;line-height:1.25}h2{font-size:1.3rem;margin-top:2rem}h3{font-size:1.1rem}p,li{max-width:74ch;line-height:1.55}a{color:oklch(42% .13 250)}table{border-collapse:collapse;width:100%;font-size:.9rem}td,th{padding:10px 12px;border-bottom:1px solid oklch(87% .01 250);text-align:left;vertical-align:top}th{background:oklch(95% .008 250);font-weight:600}td{overflow-wrap:anywhere}td>a{white-space:nowrap}.scroll{overflow:auto}summary{cursor:pointer;font-weight:600;line-height:1.5}summary:focus-visible,a:focus-visible{outline:2px solid oklch(45% .13 250);outline-offset:3px}details{margin:14px 0}.filter{padding:14px 0;border-top:1px solid oklch(82% .01 250)}.filter>summary{font-size:1.15rem}.technical{background:oklch(96% .005 250);padding:12px}.muted,figcaption{color:oklch(45% .015 250);font-size:.85rem}.pass{color:oklch(38% .09 155)}.attention{color:oklch(42% .12 45)}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.8rem;padding:12px}img{max-width:100%;height:auto}figure{margin:20px 0}.offset-plot{max-width:640px}.offset-plot svg{width:100%;height:auto}code{overflow-wrap:anywhere}@media(max-width:650px){body{margin:20px auto;padding:0 12px}td,th{padding:8px 6px}h1{font-size:1.4rem}}"""
    head='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>'+_escape(title)+'</title><style>'+css+'</style><body><h1>'+_escape(title)+'</h1>'
    head+='<p><strong>'+str(count)+' of '+str(len(filters))+' filters science-ready.</strong> '+('All recorded filters passed.' if ready else 'Use ready filters; the others need attention.')+'</p>'
    head+='<h2>Results</h2>'+table(('Filter','Science status','Relative alignment','Absolute check','What this means'),[overview_row(k,r) for k,r in sorted(filters.items())])
    attention=[]
    for k,r in sorted(filters.items()):
        if not filter_ready(r):
            status,why,action=decision(r)
            attention.append('<li><strong>'+_escape(k)+':</strong> '+_escape(why)+' '+_escape(action)+'</li>')
    if attention:head+='<h2>Needs attention</h2><ul>'+''.join(attention)+'</ul>'
    if state.get('workflow_error'):head+='<p><strong>Workflow error:</strong> '+_escape(state['workflow_error'])+'</p>'
    detail_heading='<h2>Filter details</h2><p class="muted">Expand a filter to see its sources, applied offsets and WCS. Full diagnostic tables are inside Technical details.</p>'
    body=''.join(concise_filter(k,r,contexts.get(k,{}),_summary(k,r,state.get('workflow_error'),contexts.get(k)).replace(' id="filter-'+_escape(k)+'"','')) for k,r in sorted(filters.items()))
    tail='<footer><p class="muted">Offline report from saved pipeline results. Refreshing it does not rerun a reduction. Fit-source scatter is not an independent absolute-accuracy measurement; drizzled noise is correlated.</p></footer></body></html>'
    descriptions={}
    for key,ctx in contexts.items():
        for row in ctx['associations']:
            sid=str(row.get('source_id',''));descriptions[(key,sid)]=key+' — '+_kind(row,filters[key])+' — '+sid+'. Native HST cutouts: green selected; red rejected; orange unresolved. Selection identifies a counterpart; WCS acceptance is shown in the tables above.'
    try:
        from .report_alignment import previews
        comparisons,comparison_count=previews(paths,filters,min(300_000,max(0,LIMIT-len((head+body+tail).encode('utf-8'))-8192))) if not state.get('workflow_error') else ([],0)
    except ImportError:comparisons,comparison_count=[],0
    comparison_html=('<section><h2>Before vs after alignment</h2><p class="muted">Open a filter to inspect each exposure separately. After alignment, each source should lie on the same cyan reference cross.</p>'+''.join(comparisons)+'</section>') if comparisons else ''
    body=comparison_html+detail_heading+body
    try:figures=_figures(paths,list(filters),descriptions=descriptions,budget=min(350_000,max(0,LIMIT-len((head+body+tail).encode('utf-8'))-4096)))
    except ImportError:figures=[]
    content=head+body+('<details><summary>Diagnostic previews ('+str(len(figures))+')</summary>'+''.join(figures)+'</details>' if figures else '')+tail
    while len(content.encode('utf-8'))>=LIMIT and figures:
        figures.pop();content=head+body+('<details><summary>Diagnostic previews ('+str(len(figures))+')</summary>'+''.join(figures)+'</details>' if figures else '')+tail
    if len(content.encode('utf-8'))>=LIMIT:
        # Essential per-filter status survives even pathological large manifests.
        comparison_count=0
        body=''.join('<h2>'+_escape(k)+'</h2><p>'+_escape(r.get('status','unknown'))+'</p><p>'+_brief(r.get('error','Inspect saved diagnostic measurements before science use.'),1200)+'</p>' for k,r in sorted(filters.items()))
        content=head+body+tail
    data=content.encode('utf-8')
    if len(data)>=LIMIT:raise ValueError('Report essential text exceeds strict 1 MB limit')
    output=paths.report();output.parent.mkdir(parents=True,exist_ok=True)
    temp=output.with_suffix('.html.tmp');temp.write_bytes(data);os.replace(temp,output)
    state['report']=dict(path=str(output),bytes=len(data),limit_bytes=LIMIT,self_contained=True,embedded_figures=len(figures)+comparison_count,alignment_source_pairs=comparison_count,astrometric_details_embedded=True,science_readiness='ready' if ready else 'review_required')
    save(statefile,state)
    return output


def _spatial_summary(record):
    a=record.get('scientific_assessment',{});transfer=a.get('relative_after_tweakback',{})
    pairs=transfer.get('pairs',[])
    if not pairs:pairs=a.get('relative_astrometry',{}).get('original',{}).get('pairs',[])
    def fmt(value):return '—' if value is None else f'{value:.4g}'
    rows=[]
    for p in pairs:
        r=p.get('after',p);t=r.get('spatial_systematic',{})
        if not t:continue
        values=[p.get('exposure',''),r.get('geometry',''),str(r.get('robust_inliers',r.get('matched',''))),
            fmt(r.get('postfit_scatter_pix')),fmt(t.get('max_abs_correlation')),
            fmt(t.get('amplitude_pix'))+' ± '+fmt(t.get('amplitude_uncertainty_pix')),
            fmt(t.get('amplitude_mas')),t.get('status',''),r.get('status','')]
        rows.append('<tr>'+''.join('<td>'+_escape(v)+'</td>' for v in values)+'</tr>')
    if not rows:return ''
    headers=['Exposure','Model','Inliers','RMS (pix)','|correlation|','Systematic amplitude ± uncertainty (pix)','Amplitude (mas)','Spatial assessment','WCS validation']
    return '<h3>Spatial residual assessment</h3><p>Correlation is diagnostic only. Weighted first-order trends are evaluated over matched-source coverage, with centroid-error and leave-one-out uncertainty. Underconstrained trends do not independently fail alignment.</p><table><tr>'+''.join('<th>'+h+'</th>' for h in headers)+'</tr>'+''.join(rows)+'</table>'
