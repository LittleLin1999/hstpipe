"""Readable astrometric provenance for offline reports; no fitting or FITS edits."""
from __future__ import annotations
import html,math
from pathlib import Path


def esc(value):return html.escape(str(value),quote=True)


def number(value,digits=4):
    try:return f'{float(value):.{digits}f}' if math.isfinite(float(value)) else '—'
    except (TypeError,ValueError):return '—'


def table(headers,rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(cell)+'</td>' for cell in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def context(paths,record):
    """Read saved measured catalogs and header-only current WCS snapshots."""
    from astropy.table import Table
    from astropy.io import fits
    import numpy as np
    a=record.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}));joint=rel.get('joint_models',{})
    result=dict(assessment=a,relative=rel,joint=joint,associations=[],wcs=[],warnings=[])
    sourcefile=paths.file('intermediate','alignment','common_sources','ecsv')
    if sourcefile.exists():
        try:
            for row in Table.read(sourcefile,format='ascii.ecsv'):
                clean={}
                for key in row.colnames:
                    v=row[key];clean[key]=None if np.ma.is_masked(v) else (v.item() if hasattr(v,'item') else v)
                result['associations'].append(clean)
        except Exception as exc:result['warnings'].append('Saved source associations could not be read: '+str(exc))
    # Preserve the active *pre-fit* identity separately from the current working
    # header; input_wcs_solutions may have been inspected after relative updates.
    before={r['exposure']:r for r in rel.get('wcs_solutions',a.get('input_wcs_solutions',[]))}
    current={r['exposure']:r for r in a.get('multi_chip',[])}
    configs=record.get('observing_configurations',[])
    for inp,config in zip(record.get('inputs',[]),configs):
        name=Path(inp['path']).name;root=config.get('ROOTNAME',name.removesuffix('_flc.fits'))
        working=paths.exposure_file('intermediate','alignment',root,'flc')
        old={int(c['extver']):c for c in before.get(name,{}).get('chips',[])}
        expected={int(c['extver']):c for c in current.get(name,{}).get('chips',[])}
        snapshot=[]
        if working.exists():
            try:
                with fits.open(working,memmap=True) as hd:
                    for h in hd:
                        if h.name!='SCI':continue
                        v=int(h.header.get('EXTVER',1));header=h.header
                        mismatch=v in expected and expected[v].get('wcsname')!=header.get('WCSNAME')
                        snapshot.append(dict(exposure=name,extver=v,ccdchip=header.get('CCDCHIP'),original=old.get(v,{}).get('primary'),
                            current=header.get('WCSNAME'),snapshot_status='Header differs from recorded accepted WCS' if mismatch else 'Working header read for this report',
                            header={k:header.get(k) for k in ('WCSNAME','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','PC1_1','PC1_2','PC2_1','PC2_2','CDELT1','CDELT2','A_ORDER','B_ORDER')},
                            distortion={k:header.get(k,hd[0].header.get(k)) for k in ('IDCTAB','NPOLFILE','D2IMFILE')},alternates=old.get(v,{}).get('solutions',{})))
            except Exception as exc:result['warnings'].append(name+': working WCS snapshot unavailable: '+str(exc))
        if not snapshot:
            for v,c in expected.items() or old.items():
                snapshot.append(dict(exposure=name,extver=v,ccdchip=c.get('ccdchip'),original=old.get(v,{}).get('primary'),current=c.get('wcsname',c.get('primary')),snapshot_status='Recorded WCS identity only; numerical working header unavailable',header={},distortion=old.get(v,{}).get('fit_metadata',{}),alternates=old.get(v,{}).get('solutions',{})))
        result['wcs'].extend(snapshot)
    measurements=paths.file('intermediate','alignment','joint_measurements','ecsv')
    if measurements.exists():
        try:
            known={(r.get('exposure'),r.get('source_id')) for r in result['associations']}
            for row in Table.read(measurements,format='ascii.ecsv'):
                key=(str(row['exposure']),str(row['source_id']))
                if key not in known:
                    result['associations'].append(dict(exposure=key[0],source_id=key[1],extver=int(row['extver']),x=float(row['x']),y=float(row['y']),centroid_error_pix=float(row['centroid_error_pix']),classification='HST_only',status='measured_HST_source',reason='Saved native HST centroid; external identification unavailable'))
        except Exception as exc:result['warnings'].append('Saved HST measurement table could not be read: '+str(exc))
    result['used']={};offset=0;mask=joint.get('measurement_inlier_mask',[])
    for fit in joint.get('fits',[]):
        ids=fit.get('source_ids',[])
        for j,sid in enumerate(ids):result['used'][(fit.get('exposure'),sid)]=bool(mask[offset+j]) if offset+j<len(mask) else None
        offset+=len(ids)
    result['membership_basis']='accepted_joint_assessment' if joint.get('may_combine') else 'candidate_joint_proposal'
    if rel.get('status')=='validated' and rel.get('method')=='existing_wcs' and not joint.get('may_combine'):
        # An existing WCS may pass after a proposed joint correction was rejected.
        # Use the accepted pair inlier records, not the rejected proposal mask.
        result['used']={};result['membership_basis']='accepted_existing_wcs_checks'
        catalog_paths={r['exposure']:r['path'] for r in rel.get('detector_catalogs',[])}
        try:
            catalogs={name:Table.read(path,format='ascii.ecsv') for name,path in catalog_paths.items()}
            pairs=rel.get('pairs',rel.get('original',{}).get('pairs',[]))
            for pair in pairs:
                if pair.get('status')!='validated':continue
                exposure=pair['exposure'];reference=pair['reference']
                matched=Table.read(paths.file('intermediate','alignment','original_'+Path(exposure).stem+'_matches','ecsv'),format='ascii.ecsv')
                keep=pair.get('fit_inlier_mask',[True]*len(matched))
                for n,row in enumerate(matched):
                    for name,index in ((reference,int(row['reference_index'])),(exposure,int(row['source_index']))):
                        cat=catalogs[name];column='common_source_id' if 'common_source_id' in cat.colnames else 'source_id'
                        key=(name,str(cat[column][index]));result['used'][key]=result['used'].get(key,False) or bool(keep[n])
        except Exception as exc:result['warnings'].append('Accepted existing-WCS source membership is incomplete: '+str(exc))
    return result


def _kind(row,record):
    classification=str(row.get('classification',''))
    kind='Gaia stellar candidate' if classification=='gaia_stellar_candidate' else ('Legacy PSF-class object; not necessarily a star' if classification=='legacy_point_candidate' else 'Compact-source candidate')
    if classification=='HST_only':kind='HST-only measured source; no external identification'
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        config=record['observing_configurations'][0]
        sep=SkyCoord(float(row['catalog_ra'])*u.deg,float(row['catalog_dec'])*u.deg).separation(SkyCoord(float(config['RA_TARG'])*u.deg,float(config['DEC_TARG'])*u.deg)).arcsec
        if sep<1.:kind+='; near science target'
    except (KeyError,ValueError,TypeError,IndexError):pass
    return kind


def source_table(ctx,record):
    rel=ctx['relative'];joint=ctx['joint'];ident=ctx['assessment'].get('source_identification',{});used=ctx['used'];groups={}
    for row in ctx['associations']:
        if row.get('status')=='not_covered':continue
        groups.setdefault(str(row.get('source_id')),[]).append(row)
    selected=set(sid for _,sid in used)|set(ident.get('verified_common_ids',[]))
    selected|={sid for sid,rows in groups.items() if any(r.get('classification')=='gaia_stellar_candidate' for r in rows)}
    ids=sorted(selected,key=lambda sid:(not any(k[1]==sid for k in used),not sid.startswith('G'),sid))
    rows=[]
    for sid in ids[:32]:
        group=groups.get(sid,[]);sample=group[0] if group else {};measurements=sum(k[1]==sid and v is not False for k,v in used.items())
        accepted=ctx.get('membership_basis','').startswith('accepted_') and measurements>0
        role=('Used in accepted relative fit / frame validation' if accepted else ('Candidate in unaccepted fit' if measurements else 'Identification candidate only; no accepted fit use recorded'))
        native=[r for r in group if r.get('status')=='confirmed_common']
        rows.append([esc(sid),esc(_kind(sample,record)),number(sample.get('catalog_ra'),7)+' / '+number(sample.get('catalog_dec'),7),str(measurements) if used else str(len(native))+' associated',esc(role)])
    out='<h3>Objects used for astrometry</h3><p>Catalog coordinates identify the object only. Relative fits use independently measured HST centroids. A PSF-class catalog object is not automatically a star.</p>'
    out+=table(('Source / catalog ID','Object type','Catalog RA / Dec (degrees)','HST measurements','Role in this filter'),rows) if rows else '<p>No saved object-level fit identities are available for this run.</p>'
    if len(ids)>32:out+='<p>'+str(len(ids)-32)+' further candidate identities omitted from this compact view; accepted-fit objects are listed first.</p>'
    if not joint and rel.get('status')!='validated':out+='<p><strong>No source has an accepted relative-WCS solution in this saved run.</strong> Associated objects below remain candidates. No per-exposure joint-model record is available; check the last evaluated date. Associated candidates are not accepted astrometric references.</p>'
    details=[]
    for sid in ids:
        for r in groups.get(sid,[]):
            key=(r.get('exposure'),sid);fit_use=key in used
            status=('Used in accepted assessment' if ctx.get('membership_basis','').startswith('accepted_') else 'Candidate fit input') if fit_use and used[key] is not False else ('Robust-fit rejection' if fit_use else str(r.get('status','unknown')).replace('_',' '))
            reason=r.get('reason','')
            for recovery in ident.get('hst_association',{}).get('centroid_recovery',[]):
                if recovery.get('source_id')==sid and recovery.get('exposure')==r.get('exposure') and recovery.get('status')!='measured_candidate':
                    status=recovery['status'].replace('_',' ');reason=recovery.get('reason',reason)
            details.append([esc(sid),esc(str(r.get('exposure','')).removesuffix('_flc.fits')),'SCI '+str(r.get('extver','?')),number(r.get('x'),3)+' / '+number(r.get('y'),3),number(r.get('centroid_error_pix'),3),esc(status),esc(str(reason)[:240])])
    if details:out+='<details open><summary>Measured HST centroids and exposure-by-exposure use</summary><p>X/Y are zero-based native detector pixels, before WCS correction; the SCI extension identifies the CCD coordinate system.</p>'+table(('Object','Exposure','Extension','Measured X / Y','Centroid error (pixel)','Status','Reason'),details[:128])+'</details>'
    if len(details)>128:out+='<p>Only the first 128 object/exposure records are shown.</p>'
    return out


def offset_plot(joint):
    fits=[f for f in joint.get('fits',[])[1:] if len(f.get('shift',[]))==2 and all(math.isfinite(float(v)) for v in f['shift'])]
    if not fits:return ''
    extent=max(.05,max(abs(float(v)) for f in fits for v in f['shift']))*1.35;scale=110/extent
    svg='<svg viewBox="0 0 460 290" role="img" aria-label="Proposed relative WCS correction vectors in tangent-plane pixels"><rect width="460" height="290" fill="#f5f7fa"/><path d="M60 140H330 M195 20V260" stroke="#aab4bf"/><text x="338" y="144">dx</text><text x="201" y="20">dy</text><text x="64" y="279" font-size="12">Limits ±'+number(extent,3)+' pixel; arrows point toward reference frame</text>'
    palette=['#145da0','#b44d12','#6b4597','#1b7f5a']
    for i,f in enumerate(fits):
        x=195+float(f['shift'][0])*scale;y=140-float(f['shift'][1])*scale;color=palette[i%len(palette)];dash='' if joint.get('corrections_applied') else ' stroke-dasharray="4 3"'
        svg+=f'<path d="M195 140L{x:.2f} {y:.2f}" stroke="{color}" stroke-width="2"{dash}/><circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="{color}"/><text x="340" y="{45+22*i}" font-size="12" fill="{color}">{esc(str(f.get("exposure","")).removesuffix("_flc.fits"))}</text>'
    return '<figure class="offset-plot">'+svg+'</svg><figcaption>Proposed corrections. '+('Solid arrows: corrections applied to working copies.' if joint.get('corrections_applied') else 'Dashed arrows: fitted proposals only; original WCS retained.')+'</figcaption></figure>'


def transformations(ctx):
    rel=ctx['relative'];j=ctx['joint'];out='<h3>Offsets and corrections actually applied</h3>'
    if not j:
        out+='<p>No joint-model fit is recorded. Relative status: <strong>'+esc(rel.get('status','not run'))+'</strong>. No accepted correction is implied by a failed candidate fit.</p>'
        rows=[[esc(str(p.get('exposure','')).removesuffix('_flc.fits')),str(p.get('matched','—')),esc(p.get('status','unknown')),number(p.get('initial_rms_pix')),number(p.get('postfit_scatter_pix'))] for p in rel.get('original',{}).get('pairs',[])]
        return out+(table(('Exposure','Candidate matches','Original-WCS assessment','Initial RMS (pixel)','Post-fit RMS (pixel)'),rows) if rows else '')
    accepted_existing=rel.get('status')=='validated' and rel.get('method')=='existing_wcs' and not j.get('may_combine')
    if accepted_existing:out+='<p><strong>The original WCS passed the relative-registration checks. The joint correction below was proposed and rejected; no correction was applied.</strong></p>'
    out+='<p>Reference exposure: <strong>'+esc(j.get('fits',[{}])[0].get('exposure','unavailable'))+'</strong>. '+esc(j.get('reason',j.get('model_selection_reason','')))+'.</p>'
    out+='<p>dx/dy map measured image positions into the internal reference TAN plane. They are <strong>not commanded dithers or RA/Dec offsets</strong>. Fixed rotation means retained calibration, not a measured zero.</p>'
    rows=[]
    for i,f in enumerate(j.get('fits',[])):
        shift=f.get('shift',[None,None]);applied=rel.get('status')=='validated' and j.get('may_combine') and j.get('corrections_applied') and i>0
        rotation=(number(f.get('rotation_deg'),6)+' ± '+number(f.get('rotation_uncertainty_deg'),6)) if f.get('rotation_measured') else 'Fixed; not measured'
        state='Reference frame anchor' if not i else ('Applied' if applied else ('Not applied; existing WCS retained' if rel.get('status')=='validated' else 'Not accepted'))
        actual=(number(shift[0])+' / '+number(shift[1])) if applied else ('0 / 0 (retained)' if rel.get('status')=='validated' else '—')
        rows.append([esc(str(f.get('exposure','')).removesuffix('_flc.fits')),str(f.get('matched','—')),esc(f.get('geometry','?')),number(shift[0])+' ± '+number(f.get('dx_uncertainty_pix')),number(shift[1])+' ± '+number(f.get('dy_uncertainty_pix')),esc(rotation),esc(actual),esc(state),esc(f.get('validation_scope',''))])
    out+=table(('Exposure','Objects','Model','Fitted dx ± σ (pixel)','Fitted dy ± σ (pixel)','Rotation ± σ (degree)','Applied dx / dy','Decision','Constraint scope'),rows)
    out+='<p>'+('Rejected proposal training RMS: ' if accepted_existing else 'Joint training RMS: ')+'<strong>'+number(j.get('residual_rms_pix'))+' pixel</strong>; residual degrees of freedom: '+esc(j.get('degrees_of_freedom','—'))+'. Training residuals are not independent astrometric precision.</p>'
    comparison=j.get('paired_comparison',{})
    if comparison.get('status')=='measured':
        out+='<h4>Before/after comparison on identical matched sources</h4>'+table(('Before RMS (pixel)','After RMS (pixel)','Improvement (pixel)','Required improvement (pixel)','Exposure pairs','Paired measurements'),[[number(comparison.get('before_rms_pix')),number(comparison.get('after_rms_pix')),number(comparison.get('improvement_pix')),number(j.get('required_improvement_pix')),str(comparison.get('exposure_pairs','—')),str(comparison.get('matched_pair_measurements','—'))]])
        out+='<p>'+esc(comparison.get('basis',''))+'. '+esc(comparison.get('caveat',''))+'.</p>'
    out+=offset_plot(j)
    alternatives=[]
    for i,comparison in j.get('model_comparisons',{}).items():
        name=j.get('fits',[])[int(i)].get('exposure','') if int(i)<len(j.get('fits',[])) else str(i)
        for model,f in comparison.get('candidates',{}).items():
            if model=='shift':continue
            alternatives.append([esc(name.removesuffix('_flc.fits')),esc(model),number(f.get('rotation_deg'),6)+' ± '+number(f.get('rotation_uncertainty_deg'),6),number(f.get('scale'),7),number(f.get('postfit_scatter_pix')),esc(comparison.get('reasons',{}).get(model,f.get('status','unknown')))])
    if alternatives:out+='<details><summary>Why rotation or scale was accepted or rejected</summary>'+table(('Exposure','Candidate model','Rotation ± σ (degree)','Uniform scale','Fit RMS (pixel)','Decision'),alternatives)+'</details>'
    return out


def wcs_table(ctx):
    a=ctx['assessment'];rel=ctx['relative'];absolute=a.get('astrometry',{});tb=a.get('tweakback',{});j=ctx['joint']
    out='<h3>Which WCS is in use?</h3>'
    if absolute.get('status')!='validated':out+='<p class="warning"><strong>No validated absolute WCS.</strong> The working WCS below is relative/internal only, or the preserved input solution. It must not be described as Gaia-aligned or Legacy-aligned. TweakBack: '+esc(tb.get('status',a.get('stages',{}).get('tweakback',{}).get('status','NOT_RUN')))+'.</p>'
    else:out+='<p>Validated absolute reference: '+esc(absolute.get('reference','recorded catalog'))+'; corrected solution '+esc(absolute.get('corrected_wcs_name','unavailable'))+'.</p>'
    applied=j.get('may_combine') and j.get('corrections_applied');anchor=j.get('fits',[{}])[0].get('exposure')
    rows=[];numeric=[]
    for c in ctx['wcs']:
        name=c['exposure'];state='Relative correction applied' if applied and name!=anchor else ('Original WCS retained' if rel.get('status')=='validated' else 'Unvalidated; no accepted correction')
        rows.append([esc(name.removesuffix('_flc.fits')),'SCI '+str(c['extver'])+' / CCD '+esc(c.get('ccdchip','?')),esc(c.get('original') or 'Not recorded'),esc(c.get('current') or 'Unavailable'),esc(state),esc(c['snapshot_status'])])
        h=c['header'];refs=c['distortion']
        if h:
            matrix='<br>'.join(' / '.join(number(h.get('CD'+str(i)+'_'+str(k)),12) for k in (1,2)) for i in (1,2)) if h.get('CD1_1') is not None else 'PC/CDELT representation: '+esc(str({k:v for k,v in h.items() if k.startswith(('PC','CDELT'))}))
            numeric.append([esc(name.removesuffix('_flc.fits')),'SCI '+str(c['extver']),number(h.get('CRVAL1'),9)+' / '+number(h.get('CRVAL2'),9),number(h.get('CRPIX1'),5)+' / '+number(h.get('CRPIX2'),5),matrix,esc(str(h.get('A_ORDER')))+' / '+esc(str(h.get('B_ORDER'))),'<br>'.join(esc(k)+': '+esc(refs.get(k) or 'not recorded') for k in ('IDCTAB','NPOLFILE','D2IMFILE'))])
    out+=table(('Exposure','SCI / physical CCD','Pre-fit WCSNAME','Current working WCSNAME','Relative action','Snapshot provenance'),rows) if rows else '<p>No per-extension WCS identity is available.</p>'
    if numeric:out+='<details><summary>Numerical working WCS and calibrated distortion references</summary><p>CRVAL is the WCS reference sky coordinate, not necessarily the target or image center. CRPIX uses FITS one-based pixels. CD is degrees per pixel. These linear parameters do not replace the retained SIP and lookup-table distortion model.</p>'+table(('Exposure','SCI','CRVAL1 / CRVAL2 (degrees)','CRPIX1 / CRPIX2','CD matrix, rows 1 / 2','SIP A / B order','Distortion references'),numeric)+'</details>'
    return out


def render(ctx,record):
    return '<div class="astrometric-detail">'+source_table(ctx,record)+transformations(ctx)+wcs_table(ctx)+''.join('<p class="warning">'+esc(w)+'</p>' for w in ctx['warnings'])+'</div>'


def cross_filter_summary(filters, workflow_error=None):
    """Visible internal-fit and final-frame provenance; never infer final QA."""
    rows=[];anchors=set()
    for key,record in sorted(filters.items()):
        a=record.get('scientific_assessment',{}) if not workflow_error else {}
        c=a.get('cross_filter_alignment',a.get('stages',{}).get('cross_filter_alignment',{}))
        status=c.get('status','not_recorded');anchor=c.get('anchor_group')
        if anchor:anchors.add(c.get('network_root',anchor))
        before=c.get('before',{});fit=c.get('validation') or before
        ast=a.get('astrometry',{});mode=ast.get('reference_mode')
        if ast.get('status')=='validated' and mode=='cross_filter':frame='Shared frame transferred; inspect final validation'
        elif ast.get('status')=='validated':frame='Independent absolute fit; final cross-filter offsets not measured'
        else:frame='No validated absolute frame'
        refresh=a.get('cross_filter_reprocessing',{})
        existing=c.get('final_validation',{})
        if existing.get('status')=='validated':frame='Existing published DRC passes the cross-filter check; retained unchanged'
        elif refresh:frame='New candidate: '+str(refresh.get('status','failed'))+'; previous published image retained'
        if status=='anchor':reason='Reference chosen for this attempt; this alone does not align other filters'
        elif status=='validated':reason='Internal common-source fit passed'
        else:reason=c.get('reason',c.get('error','No saved cross-filter result'))
        if before.get('status')=='candidate_fit_only':reason+='; association has insufficient independent validation'
        rms=fit.get('postfit_scatter_pix');limit=fit.get('scatter_limit_pix')
        if rms is not None and limit is not None and rms>limit:reason+='; scatter exceeds '+number(limit,3)+' pix'
        uncertainty=fit.get('transform_uncertainty_pix')
        threshold=record.get('configuration',{}).get('options',{}).get('astrometry',{}).get('relative',{}).get('max_transform_error_pix')
        if uncertainty is not None and threshold is not None and uncertainty>threshold:reason+='; shift uncertainty exceeds '+number(threshold,3)+' pix'
        badge={'validated':'PASS (internal)','anchor':'REFERENCE','unvalidated':'FAIL','not_applicable':'NOT APPLICABLE','skipped':'SKIPPED'}.get(status,'NOT RECORDED')
        model=c.get('geometry')
        if model:reason+='; '+model+(' ('+number(before.get('rotation_deg'),4)+' deg)' if before.get('rotation_measured') else ', rotation fixed')
        final=existing.get('relative_to_anchor_rms_arcsec') if existing.get('status')=='validated' else (ast.get('final_validation',{}).get('relative_to_anchor_rms_arcsec') if mode=='cross_filter' and not refresh else None)
        cells=[key,anchor or '—',badge,str(before.get('matched',fit.get('matched','—')))+' / '+str(fit.get('robust_inliers','—')),
               number(before.get('initial_rms_pix'),3)+' → '+number(fit.get('current_rms_pix',rms),3),
               number(final*1000,1) if final is not None else '—',reason+'. '+frame]
        rows.append([esc(v) for v in cells])
    text='<section id="cross-filter"><h2>Cross-filter alignment</h2><p>Exposure alignment within a filter and absolute Gaia/Legacy calibration are separate from alignment between filters. Per-filter science readiness does not establish precise multi-filter registration.</p>'
    if len(anchors)>1:text+='<p><strong>Separate connected filter groups:</strong> '+esc(', '.join(sorted(anchors)))+'. The saved records do not establish one common frame across disconnected groups.</p>'
    text+=table(('Filter','Reference','Internal result','Candidates / inliers','Before → after (pix)','Final RMS (mas)','Meaning'),rows)
    text+='<p class="muted">Counts are matched candidates and robust fit inliers, not independently confirmed source identities. Initial offsets can be large and still correctable; acceptance uses post-fit scatter and uncertainty. A passed provisional fit alone is not a measurement of final cross-filter accuracy.</p></section>'
    return text
