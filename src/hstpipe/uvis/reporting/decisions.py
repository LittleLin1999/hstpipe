"""Brief human decisions from saved QA; no scientific state changes."""
from __future__ import annotations
from .astrometry import esc,number,table


def ready(record):
    return bool(record.get('products')) and record.get('scientific_assessment',{}).get('status')=='science_ready'


def decision(record):
    a=record.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}));ab=a.get('astrometry',{})
    if ready(record):
        if a.get('source_review',{}).get('decision')=='approved':
            return 'Ready — human reviewed','Astrometry passed. The reviewer accepted the recorded field-source flux flags; automatic warnings remain.','Use the reviewed DRC with its documented photometric limitations.'
        pairs=rel.get('pairs',[]) or rel.get('original',{}).get('pairs',[])
        if rel.get('scope')=='sparse_source_constrained' or any(p.get('validation_scope')=='local_shift_only' for p in pairs):
            return 'Ready with sparse-source limits','Final image checks passed. Relative registration is constrained locally; field rotation is not independently validated.','Use the DRC with the recorded sparse-source astrometric limitation.'
        return 'Ready','Final image passed astrometric and image-quality checks.','Use the DRC; retain its recorded limitations.'
    if record.get('error'):return 'Processing error',str(record['error'])[:240],'Resolve the recorded processing error, then resume.'
    if record.get('status')=='running':return 'Running','Reduction is still in progress.','Wait for final validation.'
    transfer=a.get('relative_after_tweakback',{})
    if transfer.get('status')=='unvalidated':return 'Blocked','Relative-registration check after WCS propagation did not pass.','Inspect the saved before/after source residuals.'
    if rel.get('status') not in ('validated','single_exposure'):
        pairs=rel.get('original',{}).get('pairs',[])
        failed=[p for p in pairs if p.get('postfit_scatter_pix') is not None and p.get('scatter_limit_pix') is not None and p['postfit_scatter_pix']>p['scatter_limit_pix']]
        if failed:
            p=max(failed,key=lambda p:p['postfit_scatter_pix']/p['scatter_limit_pix'])
            why=f"Existing-WCS pair check: {p['postfit_scatter_pix']:.2f} pixel scatter, above its {p['scatter_limit_pix']:.2f} pixel limit. No tested alternative was accepted."
        else:
            missing=[str(p.get('exposure','')).removesuffix('_flc.fits') for p in pairs if p.get('matched')==0]
            why=('No reliable common centroid in '+', '.join(missing)+'. Other exposures do not establish this exposure’s registration.') if missing else 'No reliable exposure-to-exposure transformation was accepted.'
        return 'Blocked',why,'Verify common-source identities and centroid consistency; do not force a fit.'
    if ab.get('status') not in ('validated',):return 'Blocked','Relative registration passed; absolute astrometry is not validated.','Review external matches, coverage and residuals.'
    if a.get('status') in ('final_drizzle_pending','provisional_ready'):return 'Pending','Astrometric prerequisites passed; final image validation is pending.','Resume final drizzle and scientific QA.'
    if a.get('source_validation',{}).get('status') not in (None,'validated'):
        return 'Review image QA','Astrometry passed; source-preservation checks did not pass.','Inspect source flux, shape and CR rejection before science use.'
    return 'Review required','No science-ready product was recorded for this run.','Inspect the failed or incomplete stage below.'



def image_review_section(filters,figures=None):
    """Visible review instructions from actual failed source checks, not fit status."""
    sections=[];pending=False
    for key,record in sorted(filters.items()):
        a=record.get('scientific_assessment',{});q=a.get('source_validation',{})
        failed=[source for source in q.get('sources',[]) if source.get('status')=='review_required']
        if q.get('status')!='review_required' and not failed:continue
        approval=a.get('source_review') or {};accepted=approval.get('decision')=='approved' and ready(record);pending|=not accepted
        image=(record['products'][0]['path'] if accepted else a.get('unvalidated_final')) or ''
        if isinstance(image,dict):image=image.get('path','')
        text='<h3>'+esc(key)+(' — review accepted' if accepted else '')+'</h3>'
        if accepted:text+='<p><strong>Approved by '+esc(approval.get('reviewer'))+' on '+esc(approval.get('reviewed_utc'))+'.</strong> '+esc(approval.get('reason'))+' Automated measurements below are preserved; approval does not demonstrate zero field-source flux error.</p>'
        if image:text+='<p><strong>File to inspect:</strong> <code>'+esc(str(image).split('/')[-1])+'</code> — SCI and WHT extensions.</p>'
        rows=[]
        for source in failed:
            x,y=source.get('output_x'),source.get('output_y')
            position=number(x+1,2)+', '+number(y+1,2) if x is not None and y is not None else 'Not recorded'
            apertures=[v for v in source.get('apertures',[]) if v.get('status')=='review_required']
            measurements=[]
            for v in apertures:
                loss=v.get('loss_fraction');delta=number(100*loss,2)+'% lower flux' if loss is not None else 'Flux difference not measured'
                measurements.append(number(v.get('radius_pix'),1)+'-pixel radius: '+delta+'; '+number(v.get('loss_significance'),2)+'σ')
            reasons={'significant_persistent_source_flux_loss':'Source flux deficit',
                     'persistent_source_broadening':'Source broadening',
                     'source_centroid_inconsistency':'Centroid disagreement'}
            cause='; '.join(reasons.get(v,str(v).replace('_',' ')) for v in source.get('failures',[]))
            if not measurements:measurements=[cause or 'See source QA details']
            agreeing=[number(v.get('radius_pix'),1) for v in source.get('apertures',[]) if v.get('status')=='consistent']
            note=('Other aperture radii agreeing within uncertainty: '+', '.join(agreeing)+' pixels.') if agreeing else 'No agreeing aperture recorded.'
            rows.append([esc(source.get('source_id','Unknown')),esc(position),esc(cause),esc('; '.join(measurements)),esc(note)])
        text+=table(('Flagged source','SCI X, Y (1-based)','Reason','Failed measurement','Other apertures'),rows) if rows else '<p>No per-source failure details were recorded; inspect the saved QA diagnostics.</p>'
        for source in failed:text+=(figures or {}).get((key,str(source.get('source_id',''))),'')
        text+='<p><strong>What to inspect:</strong> Locate the source at the SCI position above. Compare its core in the individual exposures, independent no-CR reference, final SCI and WHT. Inspect CR masks for rejection of real source pixels, residual cosmic rays or source broadening.</p>'
        text+='<p><strong>Decision:</strong> '+('This image is published as human-reviewed science-ready data, with the field-source caveat retained. ' if accepted else 'These are review images, not validated science products. ')+ ' A flux discrepancy alone does not prove CR overrejection: sampling and median-versus-weighted-combination differences also matter. Assess precision-photometry suitability for the intended source and aperture; do not change validation flags by hand.</p>'
        sections.append(text)
    return '<section class="review-required"><h2>'+('Human review required — inspect these sources' if pending else 'Human review accepted — retained QA evidence')+'</h2>'+''.join(sections)+'</section>' if sections else ''


def overview_row(key,record):
    a=record.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}));ab=a.get('astrometry',{});v=ab.get('final_validation') or {}
    state,why,action=decision(record)
    relative='Passed' if rel.get('status')=='validated' else ('Single exposure' if rel.get('status')=='single_exposure' else 'Unvalidated')
    absolute=({'gaia_dr3':'Gaia DR3','legacy_dr10':'Legacy Surveys','cross_filter':'Cross-filter anchor'}.get(ab.get('reference'),str(ab.get('reference','')).replace('_',' '))+'; '+number(v.get('rms_arcsec'),4)+' arcsec RMS') if ab.get('status')=='validated' else 'Not validated'
    return ['<a href="#filter-'+esc(key)+'">'+esc(key)+'</a>','<strong class="'+('pass' if ready(record) else 'attention')+'">'+esc(state)+'</strong>',esc(relative),esc(absolute),esc(why)]


def concise_filter(key,record,ctx,technical):
    a=record.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}));ab=a.get('astrometry',{});joint=rel.get('joint_models',{})
    state,why,action=decision(record);out='<details class="filter" id="filter-'+esc(key)+'"><summary>'+esc(key)+' · '+esc(state)+'</summary><p>'+esc(why)+'</p>'
    out+='<p><strong>Next step:</strong> '+esc(action)+'</p>'
    products=[p['path'].split('/')[-1] for p in record.get('products',[])]
    if ready(record):out+='<p><strong>Science product:</strong> <code>'+esc(', '.join(products))+'</code></p>'
    # Keep identity and acceptance separate. A detected star is not a solved WCS.
    used=ctx.get('used',{});accepted=rel.get('status')=='validated' and ctx.get('membership_basis','').startswith('accepted_');groups={}
    for r in ctx.get('associations',[]):
        if r.get('status')=='confirmed_common':groups.setdefault(str(r['source_id']),[]).append(r)
    selected=[(sid,rs) for sid,rs in groups.items() if len({r.get('exposure') for r in rs})>=2]
    selected.sort(key=lambda x:(not any(s==x[0] and keep is not False for (_,s),keep in used.items()),not x[0].startswith('G'),x[0]))
    stars=sum(any(r.get('classification')=='gaia_stellar_candidate' for r in rs) for _,rs in selected)
    out+='<h3>Sources</h3><p>'+str(stars)+' Gaia stellar candidates and '+str(len(selected)-stars)+' other common-source candidates identified. '
    if rel.get('status')!='validated':out+='These detections have not produced an accepted relative solution.'
    elif not accepted:out+='The relative solution passed, but exact accepted source membership is not fully recorded.'
    else:out+='Accepted-fit membership is listed separately from identification.'
    out+='</p>'
    recovered=[r for r in a.get('source_identification',{}).get('hst_association',{}).get('centroid_recovery',[]) if r.get('status')=='measured_candidate' and len(r.get('mask_attempts',[]))>1]
    if recovered:
        out+='<p><strong>Flag-aware centroid recovery:</strong> '+esc(', '.join(r['exposure'].removesuffix('_flc.fits') for r in recovered))+'. The original local measurement failed; excluding inherited 4096 pixels recovered a stable centroid with a valid core. Raw pixels and WCS were not changed by recovery. Registration status is assessed separately.</p>'
    rows=[];total=len(record.get('inputs',[]))
    for sid,rs in selected[:6]:
        n=len({r.get('exposure') for r in rs});fit=sum(s==sid and keep is not False for (_,s),keep in used.items())
        role='Accepted assessment' if accepted and fit else 'Identified; fit use not established'
        rows.append([esc(sid),'Gaia stellar candidate' if any(r.get('classification')=='gaia_stellar_candidate' for r in rs) else 'Compact / other source',f'{n}/{total} exposures',esc(role)])
    if rows:out+=table(('Source ID','Type','Measured in','Use'),rows)
    if len(selected)>6:out+='<p>'+str(len(selected)-6)+' additional identities are in Technical details.</p>'
    out+='<h3>Applied corrections and current WCS</h3>'
    fits=joint.get('fits',[]);applied=bool(rel.get('status')=='validated' and joint.get('may_combine') and joint.get('corrections_applied'));rows=[]
    for i,f in enumerate(fits):
        shift=f.get('shift',[None,None]);change=(number(shift[0],3)+', '+number(shift[1],3)) if applied and i else ('Reference' if i==0 else ('0, 0 (retained)' if rel.get('status')=='validated' else 'Not applied'))
        rotation=number(f.get('rotation_deg'),5)+'°' if applied and i and f.get('rotation_measured') else 'Fixed / not applied'
        rows.append([esc(str(f.get('exposure','')).removesuffix('_flc.fits')),esc(change),esc(rotation)])
    if rows:out+=table(('Exposure','Applied dx, dy (pixels)','Rotation'),rows)
    else:out+='<p>'+('Accepted '+esc(rel.get('method','relative WCS'))+'. Exposure-level corrections are in Technical details.' if rel.get('status')=='validated' else 'No relative correction accepted; original WCS retained.')+'</p>'
    out+='<p><strong>Absolute WCS:</strong> '+esc(ab.get('corrected_wcs_name','No accepted absolute solution'))+'. '
    out+='Corrections above use the internal tangent-plane pixels, not commanded dithers.</p>'
    detection=record.get('configuration',{}).get('options',{}).get('astrometry',{}).get('detection',{})
    selection=ab.get('validation',{}).get('source_selection') or {}
    if selection.get('method')=='compact_stable' or detection.get('method')=='compact_stable':
        snr=selection.get('min_hst_snr',detection.get('compact_min_snr',10.))
        change=selection.get('max_centroid_change_pix',detection.get('compact_max_centroid_change_pix',.5))
        out+='<p><strong>Absolute reference selection:</strong> Compact sources passing pixel/morphology checks in both segmentation kernels; HST S/N ≥'+number(snr,1)+', centroid change ≤'+number(change,2)+' pixel. Selection is independent of external fit residuals.</p>'
    elif a.get('compact_absolute_selection'):
        # Read historical resume manifests without requiring their ad-hoc field
        # in newly generated reports.
        out+='<p><strong>Absolute reference selection:</strong> '+esc(a['compact_absolute_selection'].get('description','Compact, stable HST sources only.'))+'</p>'
    if a.get('source_validation'):
        q=a['source_validation'];out+='<p><strong>Automated image checks:</strong> '+esc(q.get('status','not run'))+'; '+str(q.get('validated_sources',0))+' sources checked for flux, shape and CR rejection.</p>'
    out+='<p class="muted">Last reduction: '+esc(record.get('started_utc','not recorded'))+'. Refreshing this report does not rerun the reduction.</p>'
    # All numerical provenance remains available without dominating the page.
    out+='<details class="technical"><summary>Technical details: stages, all sources, offsets, WCS and QA</summary>'+technical.replace('<details open>','<details>')+'</details></details>'
    return out
