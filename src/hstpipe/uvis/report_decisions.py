"""Brief human decisions from saved QA; no scientific state changes."""
from __future__ import annotations
from .report_astrometry import esc,number,table


def ready(record):
    return bool(record.get('products')) and record.get('scientific_assessment',{}).get('status')=='science_ready'


def decision(record):
    a=record.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}));ab=a.get('astrometry',{})
    if ready(record):return 'Ready','Final image passed astrometric and image-quality checks.','Use the DRC; retain its recorded limitations.'
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
        else:why='No reliable exposure-to-exposure transformation was accepted.'
        return 'Blocked',why,'Verify common-source identities and centroid consistency; do not force a fit.'
    if ab.get('status') not in ('validated',):return 'Blocked','Relative registration passed; absolute astrometry is not validated.','Review external matches, coverage and residuals.'
    if a.get('status') in ('final_drizzle_pending','provisional_ready'):return 'Pending','Astrometric prerequisites passed; final image validation is pending.','Resume final drizzle and scientific QA.'
    if a.get('source_validation',{}).get('status') not in (None,'validated'):
        return 'Review image QA','Astrometry passed; source-preservation checks did not pass.','Inspect source flux, shape and CR rejection before science use.'
    return 'Review required','No science-ready product was recorded for this run.','Inspect the failed or incomplete stage below.'


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
    if a.get('source_validation'):
        q=a['source_validation'];out+='<p><strong>Image checks:</strong> '+esc(q.get('status','not run'))+'; '+str(q.get('validated_sources',0))+' sources checked for flux, shape and CR rejection.</p>'
    out+='<p class="muted">Last reduction: '+esc(record.get('started_utc','not recorded'))+'. Refreshing this report does not rerun the reduction.</p>'
    # All numerical provenance remains available without dominating the page.
    out+='<details class="technical"><summary>Technical details: stages, all sources, offsets, WCS and QA</summary>'+technical.replace('<details open>','<details>')+'</details></details>'
    return out
