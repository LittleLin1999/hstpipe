"""Measured relative WCS selection/refinement; failed candidates never commit."""
import os,json,tempfile
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from .wcs import detector_wcs
from astropy.table import Table
from .relative import inspect_wcs
from .sources import sources,project,plot_sources
from .relative_fit import match,catalog_match,quality,MINIMUM
from ...common.target_runner import save
from ..drizzle import drizzle
from ..qa import image_plot
from ..chips import match_provenance


def evaluate(inputs,catalogs,options,paths,label,histogram=True,geometry='shift'):
    with fits.open(inputs[0]) as hd:tangent=WCS(hd['SCI'].header,hd)
    options=dict(options,pixel_scale_arcsec=float(np.sqrt(abs(np.linalg.det(tangent.pixel_scale_matrix)))*3600))
    xy=[project(c,p,tangent) for p,c in zip(inputs,catalogs)];pairs=[]
    for i in range(1,len(inputs)):
        try:
            ri,ii=catalog_match(xy[0],xy[i],catalogs[0],catalogs[i],options,histogram,matcher=match)
            error=np.hypot(catalogs[0]['centroid_error_pix'][ri],catalogs[i]['centroid_error_pix'][ii])
            q=quality(xy[0][ri],xy[i][ii],error,options,geometry,identities_verified='common_source_id' in catalogs[0].colnames and 'common_source_id' in catalogs[i].colnames)
            delta=xy[0][ri]-xy[i][ii]
            match_provenance(catalogs[0],catalogs[i],ri,ii,delta,q,options)
            table=Table(dict(reference_index=ri,source_index=ii,reference_extver=catalogs[0]['extver'][ri],source_extver=catalogs[i]['extver'][ii],dx_pix=delta[:,0],dy_pix=delta[:,1],centroid_error_pix=error))
            table.write(paths.file('intermediate','alignment',label+'_'+inputs[i].stem+'_matches','ecsv'),format='ascii.ecsv',overwrite=True)
            if len(delta):image_plot([delta],[label+' tangent-plane residuals (pixels)'],paths.file('diagnostics','alignment',label+'_'+inputs[i].stem+'_residuals','png'))
        except Exception as exc:q=dict(status='matching_failed',error=f'{type(exc).__name__}: {exc}',matched=0)
        q.update(exposure=inputs[i].name,reference=inputs[0].name);pairs.append(q)
    ok=bool(pairs) and all(p['status']=='validated' for p in pairs)
    network=None
    if not ok and len(inputs)>2:
        network=overlap_graph(inputs,catalogs,options,paths,label+'_network')
        ok=network['status']=='connected' and network['closure']['status']!='inconsistent'
    return dict(network=network,status='validated' if ok else 'unvalidated',pairs=pairs,may_combine=ok,reference=inputs[0].name,
        metric='full distortion-aware WCS to common undistorted TAN pixels',geometry=geometry)


def overlap_graph(inputs,catalogs,options,paths,label):
    """Report all internally supported pairs, including disconnected components.

    Refinement-supported edges are proposals, not established registration.
    Direct-reference validation remains conservative for chain-only fields.
    """
    with fits.open(inputs[0]) as hd:tangent=WCS(hd['SCI'].header,hd)
    options=dict(options,pixel_scale_arcsec=float(np.sqrt(abs(np.linalg.det(tangent.pixel_scale_matrix)))*3600))
    xy=[project(c,p,tangent) for p,c in zip(inputs,catalogs)];edges=[]
    for i in range(len(inputs)):
        for j in range(i+1,len(inputs)):
            try:
                ri,ii=catalog_match(xy[i],xy[j],catalogs[i],catalogs[j],options,True,matcher=match)
                error=np.hypot(catalogs[i]['centroid_error_pix'][ri],catalogs[j]['centroid_error_pix'][ii])
                measured=quality(xy[i][ri],xy[j][ii],error,options,identities_verified='common_source_id' in catalogs[i].colnames and 'common_source_id' in catalogs[j].colnames)
                match_provenance(catalogs[i],catalogs[j],ri,ii,xy[i][ri]-xy[j][ii],measured,options)
            except Exception as exc:measured=dict(status='matching_failed',matched=0,error=str(exc))
            edges.append(dict(reference=inputs[i].name,exposure=inputs[j].name,**measured))
    def components(statuses):
        pending={p.name for p in inputs};groups=[]
        while pending:
            group={min(pending)};changed=True
            while changed:
                changed=False
                for edge in edges:
                    pair={edge['reference'],edge['exposure']}
                    if edge['status'] in statuses and pair&group and not pair<=group:group|=pair;changed=True
            pending-=group;groups.append(sorted(group))
        return groups
    cycles=[]
    for i in range(len(inputs)):
        for j in range(i+1,len(inputs)):
            for k in range(j+1,len(inputs)):
                names=[(inputs[a].name,inputs[b].name) for a,b in ((i,j),(j,k),(i,k))]
                edgefits=[next((e for e in edges if (e['reference'],e['exposure'])==pair and 'shift' in e),None) for pair in names]
                if any(f is None for f in edgefits):continue
                delta=np.asarray(edgefits[0]['shift'])+np.asarray(edgefits[1]['shift'])-np.asarray(edgefits[2]['shift'])
                limit=max(options.get('closure_floor_pix',.3),3*np.sqrt(sum(f['transform_uncertainty_pix']**2 for f in edgefits)))
                cycles.append(dict(exposures=[inputs[a].name for a in (i,j,k)],closure_pix=delta.tolist(),limit_pix=float(limit),validated=bool(np.linalg.norm(delta)<=limit)))
    closure=dict(status=('consistent' if all(c['validated'] for c in cycles) else 'inconsistent') if cycles else 'not_constrained',cycles=cycles)
    current=components({'validated'});supported=components({'validated','refinement_supported'})
    result=dict(closure=closure,edges=edges,validated_components=current,fit_supported_components=supported,
        status='connected' if len(current)==1 else 'disconnected',external_catalog_used=False,
        caveat='Candidate associations require source-quality validation; chain-only fields are not silently accepted')
    save(paths.file('intermediate','alignment',label+'_overlap_graph'),result)
    return result


def internal_mosaic(inputs,prefix,reference=None):
    opts=dict(static=False,driz_separate=False,median=False,blot=False,driz_cr=False,driz_combine=True,
        final_pixfrac=1.,final_kernel='square',final_wht_type='EXP',runfile=str(prefix)+'.log')
    if reference:opts.update(final_refimage=f'{reference}[SCI,1]',final_wcs=True)
    else:opts['final_scale']=.04
    drizzle(inputs,prefix,opts)
    return prefix.with_name(prefix.name+'_drc.fits')


def image_quality(before,after,reference_catalog,options,paths,label,after_catalog=None,minimum_sources=None):
    """Independent no-CR coadds: source shape, flux and coverage comparisons."""
    measures=[];phase_cuts=[]
    detections=[];source_records=[]
    for phase,(image,cat) in enumerate(zip((before,after),(reference_catalog,reference_catalog if after_catalog is None else after_catalog))):
        widths=[];elong=[];flux=[];indices=[]
        with fits.open(image) as hd:
            h=hd['SCI'];w=WCS(h.header,hd);a=h.data;weight=hd['WHT'].data
            pos=w.all_world2pix(np.column_stack([cat['ra'],cat['dec']]),0)
            scale=float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600)
            phase_records=[];cutouts={}
            for source_index,(x,y) in enumerate(pos):
                extent=float(cat['measurement_extent_arcsec'][source_index]) if 'measurement_extent_arcsec' in cat.colnames else 0.
                x,y=int(round(x)),int(round(y));r=7+int(np.ceil(extent/scale))
                detail=dict(source_id=str(cat['source_id'][source_index]) if 'source_id' in cat.colnames else str(source_index),
                            radius_pix=r,status='insufficient_coverage')
                phase_records.append(detail)
                if x<r or y<r or x+r>=a.shape[1] or y+r>=a.shape[0]:continue
                cut=a[y-r:y+r+1,x-r:x+r+1].astype(float);g=weight[y-r:y+r+1,x-r:x+r+1]>0
                if not g.all():continue
                border=np.concatenate([cut[0],cut[-1],cut[:,0],cut[:,-1]]);sky=np.median(border)
                noise=max(float(1.4826*np.median(abs(border-sky))),1e-12)
                netflux=float(cut.sum()-sky*cut.size)
                detail['status']='insufficient_snr'
                if netflux/(noise*np.sqrt(cut.size))<options.get('qa_min_snr',5.):continue
                yy,xx=np.indices(cut.shape);rr=(xx-r)**2+(yy-r)**2
                local_weight=weight[y-r:y+r+1,x-r:x+r+1]
                outer=float(np.median(local_weight[rr>=5**2]))
                detail['status']='insufficient_core_coverage'
                if outer<=0 or np.mean(local_weight[rr<=2.5**2])/outer<options.get('qa_min_core_weight_ratio',.8):continue
                # Noise clipped to positive values has a broad apparent FWHM.
                # Only significant emission may support this image-quality test.
                support=np.ones(cut.shape,bool)
                if 'measurement_world' in cat.colnames:
                    positions=w.all_world2pix(np.asarray(json.loads(str(cat['measurement_world'][source_index]))),0)
                    support=np.zeros(cut.shape,bool)
                    for px,py in positions:
                        support|=(xx-(px-x+r))**2+(yy-(py-y+r))**2<=7**2
                signal=np.where((cut-sky>2*noise)&support,cut-sky,0);total=signal.sum()
                if total<=0:continue
                xc=(xx*signal).sum()/total;yc=(yy*signal).sum()/total
                width=float(cat['width_pix'][source_index]) if 'width_pix' in cat.colnames else 1.
                detail['status']='emission_not_centered'
                if np.hypot(xc-r,yc-r)>(max(1.,width*.5) if phase else r-1):continue
                cov=np.array([[(signal*(xx-xc)**2).sum(),(signal*(xx-xc)*(yy-yc)).sum()],[(signal*(xx-xc)*(yy-yc)).sum(),(signal*(yy-yc)**2).sum()]])/total
                ev=np.linalg.eigvalsh(cov);indices.append(source_index)
                widths.append(float(2.355*np.sqrt(ev.mean())));elong.append(float(np.sqrt(ev[-1]/max(ev[0],1e-9))));flux.append(float(cut.sum()-sky*cut.size))
                detail.update(status='measured',fwhm_pix=widths[-1],elongation=elong[-1],flux=flux[-1])
                cutouts[source_index]=cut
            phase_cuts.append(cutouts)
            source_records.append(phase_records)
            detections.append({i:(w,e,f) for i,w,e,f in zip(indices,widths,elong,flux)})
            measures.append(dict(sources=len(widths),fwhm_pix=float(np.median(widths)) if widths else None,elongation=float(np.median(elong)) if elong else None,flux=float(np.median(flux)) if flux else None,covered_pixels=int((weight>0).sum())))
    common=set(detections[0])&set(detections[1])
    for metrics,detected in zip(measures,detections):
        metrics['sources']=len(common)
        for j,key in enumerate(('fwhm_pix','elongation','flux')):metrics[key]=float(np.median([detected[i][j] for i in common])) if common else None
    a,b=measures
    passed=min(a['sources'],b['sources'])>=(options['qa_min_sources'] if minimum_sources is None else minimum_sources)
    if passed:passed=(b['fwhm_pix']<=a['fwhm_pix']*options['max_width_ratio'] and b['elongation']<=a['elongation']*options['max_width_ratio'])
    if passed and a['flux'] and b['flux']:passed=b['flux']>=a['flux']*.95
    cuts=[];labels=[]
    for i in sorted(common)[:3]:
        sid=str(reference_catalog['source_id'][i]) if 'source_id' in reference_catalog.colnames else str(i)
        for phase,name in enumerate(('Before','After')):
            cuts.append(phase_cuts[phase][i]);labels.append(name+' '+sid)
    if cuts:image_plot(cuts,labels,paths.file('diagnostics','alignment',label+'_coadd_cutouts','png'))
    return dict(status='validated' if passed else 'image_quality_unvalidated',before=a,after=b,
        sources_by_phase=dict(before=source_records[0],after=source_records[1]),
        aperture_policy='Same angular flux aperture enclosing all source positions plus seven output pixels; shape moments use seven-pixel neighborhoods of measured HST positions to exclude unrelated outer pixels',
        no_cr_rejection=True,source_selection='Requires aperture SNR, centered significant emission and adequate EXP core coverage in both coadds',caveat='Cutout moments are not a PSF model; SNR is an approximate pixel-noise diagnostic and does not model drizzle covariance')


def apply_native_fit(image,reference,fit):
    """Adapt tweakwcs's centered zero-origin fit to updatehdr's FITS convention."""
    from drizzlepac.updatehdr import updatewcs_with_shift
    matrix=np.asarray(fit['matrix']);center=np.asarray(fit['center']);shift=np.asarray(fit['shift'])
    inverse=np.linalg.inv(matrix);coefficients=inverse.T
    offset=center-inverse@(center+shift)+np.ones(2)-inverse@np.ones(2)
    native_shift=offset+reference.wcs.crpix@coefficients-reference.wcs.crpix
    updatewcs_with_shift(str(image),reference,wcsname='HP-RELATIVE',fitgeom=fit['geometry'],fit=coefficients,sciext='SCI',
        xsh=float(native_shift[0]),ysh=float(native_shift[1]),force=False)


def native_catalog_fit(inputs,catalogs,options,geometry='shift',reference_image=None,reference_catalog=None,guesses=None):
    """HAP's maintained matcher/fitter; update every SCI, including source-empty chips."""
    from stwcs.wcsutil import HSTWCS
    from .centroids import ensure_catalog
    catalogs=[ensure_catalog(p,c,options['detection']) for p,c in zip(inputs,catalogs)]
    refimage=reference_image or inputs[0];refcat=reference_catalog if reference_catalog is not None else catalogs[0]
    refcat=ensure_catalog(refimage,refcat,options['detection'])
    with fits.open(refimage) as hd:version=hd['SCI'].header.get('EXTVER',1)
    reference=HSTWCS(str(refimage),ext=('SCI',version))
    refxy=project(refcat,refimage,reference);planned=[]
    for i,(p,c) in enumerate(zip(inputs,catalogs)):
        if reference_image is None and i==0:continue
        xy=project(c,p,reference)
        ri,ii=catalog_match(refxy,xy,refcat,c,options,guesses is None,None if guesses is None else guesses[i],matcher=match)
        error=np.hypot(refcat['centroid_error_pix'][ri],c['centroid_error_pix'][ii])
        fit=quality(refxy[ri],xy[ii],error,options,geometry,identities_verified='common_source_id' in refcat.colnames and 'common_source_id' in c.colnames)
        if fit['status'] not in ('validated','refinement_supported'):raise ValueError('Native catalog fit lacks validated source support: '+str(fit))
        planned.append((p,fit))
    for p,fit in planned:apply_native_fit(p,reference,fit)
    return [dict(fit,exposure=p.name) for p,fit in planned]


def tweak(inputs,catalogs,work,options,geometry,histogram,reference_catalog=None,guess=None):
    from drizzlepac import tweakreg
    from .centroids import ensure_catalog,QUALITY_VERSION
    catalogs=[ensure_catalog(p,c,options['detection']) for p,c in zip(inputs,catalogs)]
    if reference_catalog is not None and reference_catalog.meta.get('centroid_quality',{}).get('version')!=QUALITY_VERSION:raise ValueError('Custom reference catalog requires pixel-level centroid validation on its reference image')
    catlines=[]
    for p,cat in zip(inputs,catalogs):
        names=[]
        with fits.open(p) as hd:
            for h in hd:
                if h.name!='SCI':continue
                v=h.header.get('EXTVER',1);c=cat[cat['extver']==v];f=work/(p.stem+f'_sci{v}.cat')
                if len(c)<3:raise ValueError('TweakReg catalog footprint unsupported on source-empty/sparse chip; use native HAP catalog fitting')
                np.savetxt(f,np.column_stack([c['x']+1,c['y']+1,c['flux']]))
                names.append(str(f))
        catlines.append(str(p)+' '+' '.join(names))
    catfile=work/'catalogs.txt';catfile.write_text('\n'.join(catlines))
    with fits.open(inputs[0]) as hd:pixel_scale=float(np.sqrt(abs(np.linalg.det(WCS(hd['SCI'].header,hd).pixel_scale_matrix)))*3600)
    kw=dict(catfile=str(catfile),refimage=str(inputs[0]),updatewcs=False,updatehdr=True,wcsname='HP-RELATIVE',
        interactive=False,see2dplot=False,residplot='No plot',clean=False,writecat=True,fitgeometry=geometry,
        minobj=MINIMUM[geometry],searchrad=options['search_radius_pix']*pixel_scale,searchunits='arcseconds',
        tolerance=options['match_tolerance_pix'],separation=options['separation_pix'],use2dhist=histogram,
        shiftfile=True,outshifts='relative_shifts.txt',outwcs='relative_reference.fits',runfile='relative_tweakreg.log')
    if reference_catalog is not None:
        ref=work/'mosaic_reference.cat';np.savetxt(ref,np.column_stack([reference_catalog['ra'],reference_catalog['dec']]))
        kw.update(refcat=str(ref),refxyunits='degrees');kw.pop('refimage')
    else:
        # TweakReg's separate refimage does not inherit the input catfile.
        # Supply its independently measured HST centroids explicitly as well.
        world=np.empty((len(catalogs[0]),2))
        with fits.open(inputs[0]) as hd:
            for v in np.unique(catalogs[0]['extver']):
                use=catalogs[0]['extver']==v
                world[use]=detector_wcs(hd,int(v)).all_pix2world(np.column_stack([catalogs[0]['x'][use],catalogs[0]['y'][use]]),0)
        ref=work/'internal_hst_reference.cat';np.savetxt(ref,world)
        kw.update(refcat=str(ref),refxyunits='degrees')
    if guess is not None:kw.update(xoffset=float(guess[0]),yoffset=float(guess[1]))
    cwd=Path.cwd()
    try:os.chdir(work);tweakreg.TweakReg([str(p) for p in inputs],**kw)
    finally:os.chdir(cwd)
    shifts=work/'relative_shifts.txt'
    if not shifts.exists():raise ValueError('No native transformation table')
    for line in shifts.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            v=[float(x) for x in line.split()[1:7]]
            if len(v)!=6 or not np.isfinite(v).all():raise ValueError('Nonfinite native transformation: no correction accepted')
            dx,dy,rot,scale,_,_=v
            if np.hypot(dx,dy)>options['max_shift_pix'] or abs((rot+180)%360-180)>options['max_rotation_deg'] or abs(scale-1)>options['max_scale_change']:raise ValueError('Native transformation exceeds configured physical bounds')


def refine_relative(inputs,paths,options,identification=None):
    report=dict(status='unvalidated',may_combine=False,method=None,attempts=[],wcs_solutions=inspect_wcs(inputs),warnings=[])
    full_files=(identification or {}).get('full_field',{}).get('catalogs',[])
    unified=[Table.read(p,format='ascii.ecsv') for p in full_files] if full_files else []
    point_catalogs=unified or [sources(p,options['detection'],paths=paths,label=p.stem+'_points') for p in inputs]
    identified=identification.get('catalogs',[]) if identification and options.get('use_verified_common_sources',True) else []
    catalogs=identified if identified and max(map(len,identified),default=0)>=1 else point_catalogs
    for p,c in zip(inputs,catalogs):
        c.write(paths.file('intermediate','alignment',p.stem+'_relative_points','ecsv'),format='ascii.ecsv',overwrite=True)
        plot_sources(p,c,paths.file('diagnostics','alignment',p.stem+'_relative_sources','png'))
    def completed(cats):
        report['result_state']=('validated_existing' if report.get('method') in ('existing_wcs','joint_existing_wcs') else 'validated_refined') if report.get('may_combine') else ('local_shift_only' if report.get('local_registration',{}).get('status')=='local_only' else 'unvalidated')
        report['scope']=report.get('scope','full_field') if report.get('may_combine') else 'diagnostic_only'
        report['detector_catalogs']=[]
        for p,c in zip(inputs,cats):
            output=paths.exposure_file('intermediate','alignment',fits.getheader(p).get('ROOTNAME',p.stem),'relative_sources','ecsv')
            c.write(output,format='ascii.ecsv',overwrite=True)
            report['detector_catalogs'].append(dict(exposure=p.name,path=str(output)))
        if len(inputs)>1:report['overlap_graph']=overlap_graph(inputs,cats,options,paths,'accepted' if report.get('may_combine') else 'unvalidated')
        save(paths.file('intermediate','alignment','relative_decisions'),report)
        return report
    if len(inputs)==1:
        report.update(status='single_exposure',may_combine=True,warnings=['No exposure pair available: relative registration cannot be measured']);return completed(catalogs)
    # Select the reference using measured usable sources, including compact galaxies.
    anchor_catalogs=unified or [sources(p,options['detection'],True,options['kernel_widths'][-1],paths=paths,label=p.stem+'_compact_w'+str(options['kernel_widths'][-1])) for p in inputs]
    requested=options.get('reference_exposure','auto')
    selection=identified if identified and max(map(len,identified),default=0)>=1 else anchor_catalogs
    anchor=max(range(len(inputs)),key=lambda i:len(selection[i])) if requested=='auto' else next((i for i,p in enumerate(inputs) if p.name==requested or p.stem.removesuffix('_flc')==requested),None)
    if anchor is None:raise ValueError('Requested internal reference exposure is absent')
    report['reference_selection']=dict(method='most usable measured compact sources',counts=[len(c) for c in anchor_catalogs],selected=inputs[anchor].name,identified_counts=[len(c) for c in identified])
    order=[anchor]+[i for i in range(len(inputs)) if i!=anchor]
    inputs=[inputs[i] for i in order];catalogs=[catalogs[i] for i in order]
    identified=[identified[i] for i in order] if identified else []
    compact_cache={float(options['kernel_widths'][-1]):[anchor_catalogs[i] for i in order]}
    if identified and all('common_source_id' in c.colnames for c in identified) and all(len(c) for c in identified):
        try:
            joint=joint_registration(inputs,identified,options,paths)
            report['joint_models']=joint
            report['attempts'].append(dict(method='joint_HST_common_sources',assessment=joint))
            if joint.get('may_combine'):
                report.update(status='validated',may_combine=True,method=joint['method'],scope=joint['scope'],image_quality=joint.get('image_quality'))
                report['warnings'].extend(joint.get('limitations',[]))
                return completed(identified)
        except Exception as exc:
            report['joint_models']=dict(status='unvalidated',error=f'{type(exc).__name__}: {exc}',may_combine=False)
    original=evaluate(inputs,catalogs,options,paths,'original')
    report['original_overlap_graph']=overlap_graph(inputs,[anchor_catalogs[i] for i in order],options,paths,'original_compact')
    report['original']=original;report['attempts'].append(dict(method='existing_wcs',assessment=original))
    if original['may_combine']:report.update(original,method='existing_wcs');return completed(catalogs)
    if options['mode']=='original' or not options['refine_if_needed']:return completed(catalogs)
    base=paths.area('intermediate','alignment');best=None
    with paths.scratch('relative-') as directory:
        root=Path(directory);baseline=None
        # Examine each common named alternate with the SAME detector source positions.
        common=None
        for r in report['wcs_solutions']:
            names=set.intersection(*(set(c['solutions'].values()) for c in r['chips']))
            common=names if common is None else common&names
        methods=[('alternate:'+n,None) for n in sorted(n for n in (common or []) if isinstance(n,str)) if n and n not in [r['chips'][0]['primary'] for r in report['wcs_solutions']]]
        from stwcs.wcsutil import headerlet
        hnames=[]
        for p in inputs:
            try:hnames.append(set(headerlet.get_headerlet_kw_names(str(p),kw='WCSNAME')))
            except Exception:hnames.append(set())
        for name in sorted(set.intersection(*hnames)):
            if name:methods.append(('headerlet:'+name,None))
        if identified:methods.insert(0,('identified_sources',None))
        methods += [(m,w) for m in options['methods'] if m not in ('existing_wcs',) for w in (options['kernel_widths'] if m=='improved_catalogs' else [None])]
        for index,(method,width) in enumerate(methods):
            print(f'Relative strategy {method}, kernel {width}',flush=True)
            entry=dict(method=method,kernel_width=width,status='rejected');report['attempts'].append(entry)
            work=root/str(index);work.mkdir();candidate=[]
            for p in inputs:
                dest=work/p.name;paths.copy(p,dest);candidate.append(dest)
            try:
                cats=catalogs
                if method.startswith('headerlet:'):
                    for p in candidate:
                        names=headerlet.get_headerlet_kw_names(str(p),kw='WCSNAME')
                        with fits.open(p) as hd:indices=[i for i,h in enumerate(hd) if h.name=='HDRLET']
                        headerlet.restore_from_headerlet(str(p),hdrext=indices[names.index(method.split(':',1)[1])],archive=True,force=False)
                elif method.startswith('alternate:'):
                    from stwcs.wcsutil import altwcs
                    for p in candidate:
                        with fits.open(p) as hd:ext=[('SCI',h.header.get('EXTVER',1)) for h in hd if h.name=='SCI']
                        altwcs.restoreWCS(str(p),ext=ext,wcsname=method.split(':',1)[1])
                else:
                    if method=='identified_sources':cats=identified
                    if method in ('improved_catalogs','offset_histogram','iterative_mosaic','image_cross_correlation'):
                        kernel=float(width or options['kernel_widths'][-1])
                        if kernel not in compact_cache:compact_cache[kernel]=unified or [sources(p,options['detection'],True,kernel,paths=paths,label=p.stem+'_compact_w'+str(kernel)) for p in inputs]
                        cats=compact_cache[kernel]
                    for p,c in zip(inputs,cats):c.write(paths.file('intermediate','alignment',f'rel{index}_'+p.stem+'_sources','ecsv'),format='ascii.ecsv',overwrite=True)
                    measured_original=evaluate(inputs,cats,options,paths,f'rel{index}_original',True)
                    entry['before']=measured_original
                    if measured_original['may_combine']:
                        entry['status']='original_validated_with_improved_catalogs'
                        catalogs=cats;report.update(measured_original,method='existing_wcs_improved_catalogs');best=(inputs,measured_original,method);break
                    geometries=['shift','rshift','rscale'] if options['fit_geometry']=='auto' else [options['fit_geometry']]
                    # Each geometry starts at the same original detector WCS.
                    result=dict(may_combine=False)
                    if method in ('offset_histogram','identified_sources'):
                        entry['native_histogram_fits']=native_graph_fit(candidate,cats,options) if method=='identified_sources' else native_catalog_fit(candidate,cats,options)
                        entry['geometry']='shift'
                    for geom in ([] if method in ('image_cross_correlation','offset_histogram','identified_sources') else geometries):
                        if geom!='shift' and options['fit_geometry']=='auto':
                            support=evaluate(inputs,cats,options,paths,f'rel{index}_model_support',True)
                            if not any(p.get('spatial_systematic',{}).get('significant',False) and (p.get('spatial_systematic',{}).get('amplitude_pix') or 0)>=options['min_improvement_pix'] and p.get('matched',0)>=MINIMUM[geom] for p in support['pairs']):
                                entry.setdefault('skipped_models',[]).append(dict(geometry=geom,reason='No source-supported spatial residual trend'));continue
                        for p,d in zip(inputs,candidate):paths.copy(p,d)
                        try:tweak(candidate,cats,work,options,geom,method in ('offset_histogram','iterative_mosaic'))
                        except Exception as exc:
                            entry.setdefault('fit_failures',[]).append(dict(geometry=geom,error=str(exc)))
                            if method!='improved_catalogs':continue
                            for p,d in zip(inputs,candidate):paths.copy(p,d)
                            try:entry['native_catalog_fits']=native_catalog_fit(candidate,cats,options,geom)
                            except Exception as failure:
                                entry.setdefault('fit_failures',[]).append(dict(geometry=geom,error=str(failure)));continue
                        result=evaluate(candidate,cats,options,paths,f'rel{index}_{geom}',False,geom)
                        entry.update(assessment=result,geometry=geom)
                        if result['may_combine']:break
                    if method=='iterative_mosaic' and not result['may_combine']:
                        for p,d in zip(inputs,candidate):paths.copy(p,d)
                        previous_score=float('inf')
                        for iteration in range(options['max_iterations']):
                            mosaic=internal_mosaic(candidate,work/f'iteration{iteration}')
                            ref=sources(mosaic,options['detection'],True,options['kernel_widths'][-1],paths=paths,label=f'internal_mosaic_{index}_{iteration}')
                            for p,d in zip(inputs,candidate):paths.copy(p,d)
                            try:tweak(candidate,cats,work,options,'shift',True,ref)
                            except Exception as failure:
                                entry.setdefault('mosaic_native_fallback',[]).append(str(failure))
                                for p,d in zip(inputs,candidate):paths.copy(p,d)
                                entry['mosaic_fits']=native_catalog_fit(candidate,cats,options,reference_image=mosaic,reference_catalog=ref)
                            measured=evaluate(candidate,cats,options,paths,f'rel{index}_iter{iteration}',False)
                            entry.setdefault('iterations',[]).append(measured)
                            if measured['may_combine']:result=measured;break
                            current_score=np.median([p.get('current_rms_pix',float('inf')) for p in measured['pairs']])
                            if current_score>=previous_score:break
                            previous_score=current_score
                    if method=='image_cross_correlation':
                        # Independent catalog validation remains mandatory; image shifts only seed matching.
                        from skimage.registration import phase_cross_correlation
                        reference=internal_mosaic([inputs[0]],work/'correlation_reference')
                        shifts=[]
                        for p in inputs[1:]:
                            single=internal_mosaic([p],work/(p.stem+'_single'),reference)
                            with fits.open(reference) as a,fits.open(single) as b:
                                from photutils.aperture import CircularAperture
                                masks=[]
                                for hd,cat in ((a,cats[0]),(b,cats[len(shifts)+1])):
                                    good=np.zeros(hd['SCI'].shape,bool)
                                    positions=WCS(hd['SCI'].header,hd).all_world2pix(np.column_stack([cat['ra'],cat['dec']]),0)
                                    for aperture in CircularAperture(positions,r=8).to_mask(method='center'):
                                        region=aperture.to_image(good.shape)
                                        if region is not None:good|=region>0
                                    masks.append(good&(hd['WHT'].data>0))
                                shift=phase_cross_correlation(a['SCI'].data,b['SCI'].data,reference_mask=masks[0],moving_mask=masks[1])
                                if isinstance(shift,tuple):shift=shift[0]
                                common=WCS(a['SCI'].header,a);center=common.wcs.crpix-1
                                world=common.wcs_pix2world([center,center-np.asarray(shift)[::-1]],0)
                                with fits.open(p) as detector:
                                    tangent=WCS(detector['SCI'].header,detector)
                                    xy=tangent.wcs_world2pix(world,0)
                                shifts.append(xy[1]-xy[0])
                            entry['correlation_candidates_pix']=[s.tolist() for s in shifts]
                        # Candidate offsets are never applied without source-supported fitting.
                        for i,shift in enumerate(shifts,1):
                            tweak([candidate[i]],[cats[i]],work,options,'shift',False,reference_catalog=cats[0],guess=shift)
                            (work/'relative_shifts.txt').rename(work/f'correlation_{i}_shifts.txt')
                result=evaluate(candidate,cats,options,paths,f'rel{index}_final',True,entry.get('geometry','shift'))
                entry['assessment']=result
                if not result['may_combine']:continue
                graph=overlap_graph(candidate,cats,options,paths,f'rel{index}_candidate')
                entry['overlap_graph']=graph
                if graph['closure']['status']=='inconsistent':entry['status']='inconsistent_pairwise_closure';continue
                score=lambda r:np.median([p.get('current_rms_pix',np.inf) for p in r['pairs']])
                comparable=evaluate(inputs,cats,options,paths,f'rel{index}_before',True,entry.get('geometry','shift'))
                entry['before']=comparable
                if comparable['may_combine']:
                    entry['status']='original_validated_with_improved_catalogs'
                    catalogs=cats;report.update(comparable,method='existing_wcs_improved_catalogs');best=(inputs,comparable,method);break
                improvement=score(comparable)-score(result)
                entry['improvement_pix']=float(improvement) if np.isfinite(improvement) else None
                if options['require_improvement'] and improvement<options['min_improvement_pix']:entry['status']='no_significant_improvement';continue
                if baseline is None:baseline=internal_mosaic(inputs,root/'before')
                after=internal_mosaic(candidate,work/'after',baseline)
                reference_catalog=cats[0]
                after_catalog=cats[0].copy()
                with fits.open(candidate[0]) as hd:
                    for version in np.unique(after_catalog['extver']):
                        use=after_catalog['extver']==version;h=hd['SCI',int(version)]
                        world=WCS(h.header,hd).all_pix2world(np.column_stack([after_catalog['x'][use],after_catalog['y'][use]]),0)
                        after_catalog['ra'][use]=world[:,0];after_catalog['dec'][use]=world[:,1]
                if all('common_source_id' in c.colnames for c in cats):reference_catalog,after_catalog=coadd_reference_catalog(inputs,candidate,cats)
                qa=image_quality(baseline,after,reference_catalog,options,paths,f'rel{index}',after_catalog=after_catalog)
                entry['image_quality']=qa
                if qa['status']!='validated':continue
                entry['status']='accepted';best=(candidate,result,method)
                for p,c in zip(inputs,candidate):paths.copy(c,p)
                catalogs=cats;report.update(result,method=method,image_quality=qa);break
            except Exception as exc:entry.update(error=f'{type(exc).__name__}: {exc}')
            finally:
                for f in work.glob('*'):
                    if f.is_file() and (f.suffix in ('.txt','.log','.match','.cat') or f.name=='relative_reference.fits'):
                        paths.copy(f,paths.file('intermediate','alignment',f'rel{index}_'+f.stem,f.suffix.lstrip('.')))
                import shutil
                shutil.rmtree(work)
        if best is None:
            from .target_fallback import target_alignment
            fallback,local_catalogs=target_alignment(inputs,paths,options,identification=identification,catalogs=identified,previous_joint=report.get('joint_models'))
            report['local_registration']=fallback
            if fallback.get('may_combine'):
                report.update(fallback['validation'],method='target_specific',image_quality=fallback.get('image_quality'))
                return completed(local_catalogs)
            inspection=internal_mosaic(inputs,root/'unvalidated_inspection')
            with fits.open(inspection,mode='update') as hd:
                hd[0].header['HPPASS']='RELATIVE-UNVALIDATED';hd[0].header['HPQASTAT']='REVIEW'
            inspection_output=paths.file('intermediate','drizzle','relative_inspection','fits')
            paths.copy(inspection,inspection_output);report['inspection_mosaic']=str(inspection_output)
            report['warnings'].append('Relative registration not established. Original WCS retained; absolute alignment and TweakBack blocked.')
    return completed(catalogs)


def single_source_fallback(inputs,catalogs,options,paths,reference_source=None):
    """Constrained local shift diagnostic; never authorizes full-field astrometry."""
    settings=options.get('local_fallback',{})
    if not settings.get('enabled',True) or not catalogs or any('common_source_id' not in c.colnames for c in catalogs):
        return dict(status='not_available',reason='No validated common source identities')
    common=set(str(s) for s in catalogs[0]['common_source_id'])
    for cat in catalogs[1:]:common&=set(str(s) for s in cat['common_source_id'])
    if len(common)!=1:return dict(status='not_available',reason='A unique common compact source is required across all exposures')
    sid=common.pop();rows=[c[np.asarray(c['common_source_id']).astype(str)==sid] for c in catalogs]
    with fits.open(inputs[0]) as hd:
        ref=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
        target=[hd[0].header['RA_TARG'],hd[0].header['DEC_TARG']]
        manual=settings.get('reference_positions',[])
        if len(manual)==1:target=[float(manual[0]['ra']),float(manual[0]['dec'])]
        w=detector_wcs(hd,int(rows[0]['extver'][0]));world=w.all_pix2world([[rows[0]['x'][0],rows[0]['y'][0]]],0)[0]
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    separation=float(SkyCoord(*world,unit=u.deg).separation(SkyCoord(*target,unit=u.deg)).arcsec)
    if reference_source is None and separation>settings.get('max_target_separation_arcsec',3.):
        return dict(status='not_available',reason='The lone source does not constrain target-centered registration',separation_arcsec=separation)
    from tweakwcs.linearfit import iter_linear_fit
    xy=[project(c,p,ref) for c,p in zip(rows,inputs)];fits_report=[]
    for i in range(1,len(inputs)):
        fit=iter_linear_fit(xy[0],xy[i],fitgeom='shift',nclip=0)
        if not np.isfinite(fit['shift']).all() or np.linalg.norm(fit['shift'])>options['max_shift_pix']:
            return dict(status='unvalidated',reason='Nonfinite or implausible local shift')
        fits_report.append(dict(exposure=inputs[i].name,geometry='shift',center=np.asarray(fit['center']).tolist(),matrix=np.asarray(fit['matrix']).tolist(),shift=np.asarray(fit['shift']).tolist(),
            uncertainty_pix=float(np.hypot(rows[0]['centroid_error_pix'][0],rows[i]['centroid_error_pix'][0]))))
    record=dict(status='local_only',source_id=sid,reference_source=reference_source,reference_exposure=inputs[0].name,transformations=fits_report,matched_per_pair=1,independent_scatter_pix=None,full_field_validated=False,
        caveat='One source constrains translation locally; no independent rotation, scale, residual scatter or full-field accuracy')
    norms=[float(np.linalg.norm(f['shift'])) for f in fits_report]
    record.update(initial_pair_offsets_pix=norms,existing_wcs_locally_consistent=all(n<=min(options['max_offset_pix'],3*f['uncertainty_pix']) for n,f in zip(norms,fits_report)),
        corrections_applied=False,local_state='existing_wcs_consistent' if all(n<=min(options['max_offset_pix'],3*f['uncertainty_pix']) for n,f in zip(norms,fits_report)) else 'candidate_local_translation')
    if settings.get('diagnostic_coadd',False):
        with paths.scratch('local-only-') as tmp:
            copies=[]
            for p in inputs:
                dest=Path(tmp)/p.name;paths.copy(p,dest);copies.append(dest)
            for p,fit in zip(copies[1:],fits_report):apply_native_fit(p,ref,fit)
            mosaic=internal_mosaic(copies,Path(tmp)/'local')
            with fits.open(mosaic,mode='update') as hd:hd[0].header['HPPASS']='LOCAL-ONLY';hd[0].header['HPQASTAT']='REVIEW'
            output=paths.file('intermediate','drizzle','local_inspection','fits');paths.copy(mosaic,output);record['diagnostic_coadd']=str(output)
        record['diagnostic_coadd_wcs']='Candidate local translations on disposable copies only; no working FLC correction applied'
    save(paths.file('intermediate','alignment','local_registration'),record)
    return record


def native_graph_fit(inputs,catalogs,options,geometry='shift'):
    """Traverse supported overlaps using HAP's native fitter and STWCS updates.

    Each exposure is updated once on caller-owned candidate copies. Disconnected
    or imprecise chains fail, and all cycles/image quality remain acceptance gates.
    """
    from .centroids import ensure_catalog
    from .wcs import detector_wcs
    catalogs=[ensure_catalog(p,c,options['detection']) for p,c in zip(inputs,catalogs)]
    solved={0};remaining=set(range(1,len(inputs)));variance={0:0.};records=[]
    while remaining:
        choices=[]
        for ref_index in sorted(solved):
            with fits.open(inputs[ref_index]) as hd:reference=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
            a=project(catalogs[ref_index],inputs[ref_index],reference)
            for i in sorted(remaining):
                b=project(catalogs[i],inputs[i],reference)
                try:
                    ri,ii=catalog_match(a,b,catalogs[ref_index],catalogs[i],options,True,matcher=match)
                    errors=np.hypot(catalogs[ref_index]['centroid_error_pix'][ri],catalogs[i]['centroid_error_pix'][ii])
                    fit=quality(a[ri],b[ii],errors,options,geometry,identities_verified='common_source_id' in catalogs[ref_index].colnames and 'common_source_id' in catalogs[i].colnames)
                    if fit['status'] not in ('validated','refinement_supported'):continue
                    n=fit['effective_matches']
                    edge_error=max(fit['transform_uncertainty_pix'],fit['centroid_precision_pix']/np.sqrt(n),fit['robust_scatter_pix']/np.sqrt(2*n))
                    total=variance[ref_index]+edge_error**2
                    if np.sqrt(total)>options['max_transform_error_pix']:continue
                    choices.append((fit['robust_inliers'],n,ref_index,i,fit,reference,total))
                except Exception:continue
        if not choices:raise ValueError('Overlap graph disconnected, underconstrained, or propagated uncertainty exceeds limit')
        _,_,parent,i,fit,reference,total=max(choices,key=lambda r:(r[0],r[1],-r[6],-r[2],-r[3]))
        apply_native_fit(inputs[i],reference,fit)
        variance[i]=total;solved.add(i);remaining.remove(i)
        records.append(dict(fit,exposure=inputs[i].name,parent_reference=inputs[parent].name,propagated_uncertainty_pix=float(np.sqrt(total))))
    return records


def coadd_reference_catalog(before_inputs,after_inputs,catalogs):
    """One measured object per common ID, sampled across the entire overlap graph."""
    chosen={}
    for i,c in enumerate(catalogs):
        for j,row in enumerate(c):
            sid=str(row['common_source_id'])
            if sid not in chosen or float(row['centroid_error_pix'])<chosen[sid][0]:
                chosen[sid]=(float(row['centroid_error_pix']),i,j)
    before=[];after=[]
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    for sid,(_,i,j) in sorted(chosen.items()):
        row=catalogs[i][j];xy=[[float(row['x']),float(row['y'])]];v=int(row['extver'])
        values=[]
        extents=[];all_positions=[]
        for images in (before_inputs,after_inputs):
            with fits.open(images[i]) as hd:center=detector_wcs(hd,v).all_pix2world(xy,0)[0]
            values.append(center);positions=[]
            for image,cat in zip(images,catalogs):
                with fits.open(image) as hd:
                    for other in cat[np.asarray(cat['common_source_id']).astype(str)==sid]:
                        positions.append(detector_wcs(hd,int(other['extver'])).all_pix2world([[float(other['x']),float(other['y'])]],0)[0])
            all_positions.append(positions)
            sky=SkyCoord(np.asarray(positions)[:,0]*u.deg,np.asarray(positions)[:,1]*u.deg)
            extents.append(float(sky.separation(SkyCoord(*center,unit='deg')).arcsec.max()))
        for output,world,positions in zip((before,after),values,all_positions):
            output.append(dict(ra=float(world[0]),dec=float(world[1]),source_id=sid,width_pix=float(row['width_pix']),
                               measurement_extent_arcsec=max(extents),measurement_world=json.dumps(np.asarray(positions).tolist())))
    return Table(rows=before),Table(rows=after)


def paired_comparison(before,after,catalogs,options,*,fits_report=None,inlier_mask=None,names=None):
    """Compare identical native-HST associations/inliers on every overlap.

    Corrected coordinates are never rematched: a dropped source, changed
    partner, or changed exposure-pair set cannot manufacture improvement.
    """
    allowed=None
    if fits_report is not None and inlier_mask is not None:
        allowed={};offset=0
        for i,f in enumerate(fits_report):
            for j,sid in enumerate(f.get('source_ids',[])):
                allowed[(i,str(sid))]=bool(inlier_mask[offset+j])
            offset+=len(f.get('source_ids',[]))
        if offset!=len(inlier_mask):raise ValueError('Joint inlier provenance does not match source measurements')
    residuals_before=[];residuals_after=[];pairs=[]
    for i in range(len(catalogs)):
        for j in range(i+1,len(catalogs)):
            a,b=catalog_match(before[i],before[j],catalogs[i],catalogs[j],options)
            ids=np.asarray(catalogs[i]['common_source_id']).astype(str)[a]
            if allowed is not None:
                keep=np.array([allowed.get((i,sid),False) and allowed.get((j,sid),False) for sid in ids],bool)
                a,b,ids=a[keep],b[keep],ids[keep]
            if not len(a):continue
            initial=np.asarray(before[i])[a]-np.asarray(before[j])[b]
            corrected=np.asarray(after[i])[a]-np.asarray(after[j])[b]
            if not np.isfinite(np.r_[initial.ravel(),corrected.ravel()]).all():raise ValueError('Nonfinite matched-position comparison')
            rms=lambda x:float(np.sqrt(np.mean(np.sum(x*x,axis=1))))
            residuals_before.extend(initial.tolist());residuals_after.extend(corrected.tolist())
            pairs.append(dict(reference=names[i] if names else i,exposure=names[j] if names else j,source_ids=ids.tolist(),matched=len(a),before_rms_pix=rms(initial),after_rms_pix=rms(corrected)))
    if not residuals_before:return dict(status='unavailable',matched_pair_measurements=0,pairs=pairs)
    rms=lambda x:float(np.sqrt(np.mean(np.sum(np.asarray(x)**2,axis=1))))
    a,b=rms(residuals_before),rms(residuals_after)
    return dict(status='measured',before_rms_pix=a,after_rms_pix=b,improvement_pix=a-b,exposure_pairs=len(pairs),matched_pair_measurements=len(residuals_before),pairs=pairs,
        basis='Identical HST source identities, row indices, inlier measurements and exposure pairs before/after; every available overlap included',
        caveat='Pairs sharing exposures/sources are correlated; this comparison is not independent astrometric precision')


def joint_registration(inputs,catalogs,options,paths):
    """Assess common-source models before blind strategies; commit only after QA."""
    from .relative_fit import joint_models
    with fits.open(inputs[0]) as hd:reference=detector_wcs(hd,hd['SCI'].header.get('EXTVER',1))
    options=dict(options,pixel_scale_arcsec=float(np.sqrt(abs(np.linalg.det(reference.pixel_scale_matrix)))*3600))
    xy=[project(c,p,reference) for c,p in zip(catalogs,inputs)]
    result=joint_models(xy,catalogs,options)
    for i,f in enumerate(result.get('fits',[])):f['exposure']=inputs[i].name
    result.update(method='joint_HST_common_sources',may_combine=False,corrections_applied=False)
    if result['status']!='fit_supported':return result
    if result['degrees_of_freedom']<=0:
        result.update(status='limited_solution',reason='Mathematical solution exists but has no residual redundancy; independent validation unavailable')
        save(paths.file('intermediate','alignment','joint_models'),result)
        return result
    fit_ok=True;current_ok=True
    for i,f in enumerate(result['fits']):
        if not i:continue
        limit=min(options['max_scatter_pix'],max(options['scatter_floor_pix'],options['centroid_sigma']*f['centroid_precision_pix']))
        f['scatter_limit_pix']=limit
        rot=f['rotation_deg'] or 0.
        f['fit_acceptable']=bool(not f.get('spatial_systematic',{}).get('exceeds_accuracy',False) and f['postfit_scatter_pix']<=limit and f['transform_uncertainty_pix']<=options['max_transform_error_pix'] and f['outlier_fraction']<=options['max_outlier_fraction'] and np.linalg.norm(f['shift'])<=options['max_shift_pix'] and abs(rot)<=options['max_rotation_deg'] and abs(f['scale']-1)<=options['max_scale_change'])
        # Propagate correction covariance over the actual multi-chip footprint,
        # rather than judging a rotation solely at its fitted source centroid.
        field=[]
        with fits.open(inputs[i]) as hd:
            for h in hd:
                if h.name!='SCI':continue
                ny,nx=h.shape;corners=np.array([[0.,0.],[nx-1.,0.],[0.,ny-1.],[nx-1.,ny-1.]])
                world=detector_wcs(hd,int(h.header.get('EXTVER',1))).all_pix2world(corners,0)
                field.extend(reference.wcs_world2pix(world,0).tolist())
        cov=np.asarray(f['parameter_covariance']);prediction=[]
        for position in field:
            jac=np.zeros((2,len(cov)));jac[:,:2]=np.eye(2)
            q=(np.asarray(position)-f['center'])@np.asarray(f['matrix']).T
            if f['rotation_measured']:jac[:,2]=[-q[1],q[0]]
            if f['geometry']=='rscale':jac[:,3]=q
            prediction.append(float(np.sqrt(np.linalg.eigvalsh(jac@cov@jac.T).max())))
        f['max_correction_uncertainty_over_footprint_pix']=max(prediction) if prediction else None
        f['fit_acceptable'] &= bool(prediction and max(prediction)<=options['max_transform_error_pix'])
        fit_ok &= f['fit_acceptable']
        # Existing WCS is tested using actual pair differences, not latent
        # fitted source positions that would halve residuals artificially.
        checks=[]
        for j in range(i):
            a,b=catalog_match(xy[j],xy[i],catalogs[j],catalogs[i],options)
            if not len(a):continue
            error=np.hypot(catalogs[j]['centroid_error_pix'][a],catalogs[i]['centroid_error_pix'][b]);delta=xy[j][a]-xy[i][b]
            radial=np.linalg.norm(delta,axis=1);allowed=np.minimum(options['max_current_rms_pix'],np.maximum(options['scatter_floor_pix'],options['centroid_sigma']*error))
            checks.append(dict(reference=inputs[j].name,matched=len(a),rms_pix=float(np.sqrt(np.mean(radial**2))),max_offset_pix=float(radial.max()),within_measurement_errors=bool(np.all(radial<=allowed)),source_ids=np.asarray(catalogs[i]['common_source_id'])[b].tolist()))
        f['existing_pair_checks']=checks
        current_ok &= bool(checks) and all(c['within_measurement_errors'] and c['max_offset_pix']<=options['max_offset_pix'] for c in checks)
    if not fit_ok:
        result.update(status='unvalidated',reason='Joint model fails unchanged precision, residual or physical bounds');return result
    # Record all native HST measurements and their chip identity; no catalog
    # projected position becomes an observation in the adjustment.
    measurements=[];usage={};offset=0
    for i,f in enumerate(result['fits']):
        for sid in f['source_ids']:
            usage[(i,sid)]=(bool(result['measurement_inlier_mask'][offset]),result['measurement_residuals_pix'][offset]);offset+=1
    for i,(image,c) in enumerate(zip(inputs,catalogs)):
        for row in c:
            sid=str(row['common_source_id']);used,residual=usage.get((i,sid),(False,[None,None]))
            measurements.append(dict(exposure=image.name,source_id=sid,extver=int(row['extver']),
                x=float(row['x']),y=float(row['y']),centroid_error_pix=float(row['centroid_error_pix']),
                fit_inlier=used,residual_x_pix=residual[0],residual_y_pix=residual[1],
                fit_status='inlier' if used else ('robust_outlier' if (i,sid) in usage else 'not_common')))
    Table(rows=measurements).write(paths.file('intermediate','alignment','joint_measurements','ecsv'),format='ascii.ecsv',overwrite=True)
    with paths.scratch('joint-relative-') as tmp:
        directory=Path(tmp);candidates=[]
        for p in inputs:
            dest=directory/p.name;paths.copy(p,dest);candidates.append(dest)
        for p,f in zip(candidates[1:],result['fits'][1:]):apply_native_fit(p,reference,f)
        # Verify native WCS propagation on each measured source and all SCI
        # extensions through the established updater, not isolated keywords.
        corrected=[project(c,p,reference) for c,p in zip(catalogs,candidates)]
        residuals=[]
        for i in range(1,len(inputs)):
            a,b=catalog_match(corrected[0],corrected[i],catalogs[0],catalogs[i],options)
            if len(a):residuals.extend((corrected[0][a]-corrected[i][b]).tolist())
            with fits.open(candidates[i]) as hd:result['fits'][i]['updated_sci_extensions']=[int(h.header.get('EXTVER',1)) for h in hd if h.name=='SCI']
        result['propagated_pair_rms_pix']=float(np.sqrt(np.mean(np.sum(np.asarray(residuals)**2,axis=1)))) if residuals else None
        comparison=paired_comparison(xy,corrected,catalogs,options,fits_report=result['fits'],inlier_mask=result['measurement_inlier_mask'],names=[p.name for p in inputs])
        result['paired_comparison']=comparison
        result['required_improvement_pix']=options['min_improvement_pix']
        before=internal_mosaic(inputs,directory/'before');after=internal_mosaic(candidates,directory/'after',before)
        original_cat,after_cat=coadd_reference_catalog(inputs,candidates,catalogs)
        # Every measured object is inspected. Image metrics supplement rather
        # than replace centroid, model-rank and uncertainty evidence.
        qa=image_quality(before,after,original_cat,options,paths,'joint_relative',after_catalog=after_cat,minimum_sources=1)
        result['image_quality']=qa
        baseline_qa=image_quality(before,before,original_cat,options,paths,'joint_existing',minimum_sources=1)
        result['existing_image_quality']=baseline_qa
        for phase,mosaic in (('existing',before),('candidate',after)):
            out=paths.file('intermediate','drizzle','joint_'+phase+'_inspection','fits');paths.copy(mosaic,out)
            result[phase+'_inspection']=str(out)
        if current_ok and baseline_qa['status']=='validated':
            result.update(status='validated',method='joint_existing_wcs',may_combine=True,reason='Measured pair offsets satisfy existing accuracy and centroid-error limits; original WCS retained')
        elif qa['status']=='validated':
            after_rms=comparison.get('after_rms_pix')
            improvement=comparison.get('improvement_pix',-np.inf)
            result['improvement_pix']=improvement if np.isfinite(improvement) else None
            if (not options['require_improvement'] or improvement>=options['min_improvement_pix']) and after_rms is not None and after_rms<=options['max_scatter_pix']:
                for p,c in zip(inputs,candidates):paths.copy(c,p)
                result.update(status='validated',may_combine=True,corrections_applied=True,reason='Joint source consistency and no-CR image QA support correction')
        if not result['may_combine']:result.update(status='unvalidated',reason='Source-fit proposal lacks required image quality or measured improvement')
    result['scope']='sparse_source_constrained' if any(f['validation_scope'] in ('limited_redundancy','local_shift_only') for f in result['fits']) else 'common_source_constrained'
    result['limitations']=['A fixed rotation is retained calibration, not a measured zero rotation.','Joint training residuals are not independent validation; no-CR morphology is an additional consistency check.']
    save(paths.file('intermediate','alignment','joint_models'),result)
    return result

