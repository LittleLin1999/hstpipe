"""External identity priors; all relative-fitting coordinates are measured in HST.

Catalog sky positions bound identification windows only. No external sky position
is passed to the relative fitter or used as a relative-registration residual.
"""
import json,tempfile,warnings
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord,CartesianRepresentation
from astropy.io import fits
from astropy.table import Table,vstack
from astropy.wcs import WCS
from .wcs import detector_wcs
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import gaussian_filter
from .catalogs import query_gaia_dr3,query_legacy_dr10,gaia_at_epoch,legacy_at_epoch,cached_footprint,CatalogQueryError,CatalogSchemaError
from ..fits import image_data
from .centroids import pixel_masks,validate_centroids
from ...common.cache import StageCache
from ...common.target_runner import save
from ...common.provenance import identity,versions


def footprint(inputs,margin_arcsec):
    edges=[];geometry=[]
    for image in inputs:
        with fits.open(image) as hd:
            for h in hd:
                if h.name!='SCI':continue
                ny,nx=h.shape;t=np.linspace(0,1,17)
                pixels=np.vstack([np.column_stack([t*(nx-1),np.full(17,y)]) for y in (0,ny-1)]+[np.column_stack([np.full(17,x),t*(ny-1)]) for x in (0,nx-1)])
                world=detector_wcs(hd,h.header.get('EXTVER',1)).all_pix2world(pixels,0)
                if not np.isfinite(world).all():raise ValueError('Invalid FLC footprint WCS')
                edges.extend(world);geometry.append(dict(exposure=image.name,extver=h.header.get('EXTVER',1),nx=nx,ny=ny))
    sky=SkyCoord(np.array(edges)[:,0]*u.deg,np.array(edges)[:,1]*u.deg)
    center=SkyCoord(CartesianRepresentation(sky.cartesian.xyz.mean(axis=1))).spherical
    center=SkyCoord(center.lon,center.lat)
    return dict(ra=float(center.ra.deg),dec=float(center.dec.deg),radius=float(center.separation(sky).max().deg+margin_arcsec/3600),margin_arcsec=margin_arcsec,geometry=geometry)


def epoch(image):
    h=fits.getheader(image)
    if 'EXPSTART' in h:return float(Time(h['EXPSTART'],format='mjd').tcb.jyear)
    return float(Time(h['DATE-OBS']+'T'+h.get('TIME-OBS','00:00:00')).tcb.jyear)


def merge_candidates(gaia,legacy,options):
    """Mutual, unambiguous cross-identification; unmatched candidates survive."""
    rows=[];used=set();pairs={}
    if len(gaia) and len(legacy):
        a=SkyCoord(np.asarray(gaia['ra'])*u.deg,np.asarray(gaia['dec'])*u.deg);b=SkyCoord(np.asarray(legacy['ra'])*u.deg,np.asarray(legacy['dec'])*u.deg)
        j,d,_=a.match_to_catalog_sky(b);reverse,_,_=b.match_to_catalog_sky(a)
        for i in range(len(a)):
            if d[i].arcsec>options['merge_radius_arcsec'] or reverse[j[i]]!=i:continue
            if len(b)>1 and np.sort(a[i].separation(b).arcsec)[1]<=options['merge_radius_arcsec']:continue
            if len(a)>1 and np.sort(b[j[i]].separation(a).arcsec)[1]<=options['merge_radius_arcsec']:continue
            pairs[i]=int(j[i]);used.add(int(j[i]))
    for i,r in enumerate(gaia):
        j=pairs.get(i);other=legacy[j] if j is not None else None
        rows.append(dict(source_id='G'+str(r['source_id']),gaia_id=str(r['source_id']),legacy_id=str(other['source_id']) if other is not None else '',
            ra=float(r['ra']),dec=float(r['dec']),classification='gaia_stellar_candidate',legacy_type=str(other['type']) if other is not None else '',
            release='Gaia DR3'+('/'+str(other['release']) if other is not None else ''),gaia_index=i,legacy_index=j if j is not None else -1))
    for i,r in enumerate(legacy):
        if i in used:continue
        kind=str(r['type'])
        if not options['allow_compact_galaxies'] and kind!='PSF':continue
        rows.append(dict(source_id='L'+str(r['source_id']),gaia_id='',legacy_id=str(r['source_id']),ra=float(r['ra']),dec=float(r['dec']),
            classification='legacy_point_candidate' if kind=='PSF' else 'compact_galaxy_candidate',legacy_type=kind,release=str(r['release']),gaia_index=-1,legacy_index=i,identification_prior_only=bool(r['identification_prior_only']) if 'identification_prior_only' in r.colnames else False))
    return Table(rows=rows) if rows else Table(names=['source_id','gaia_id','legacy_id','ra','dec','classification','legacy_type','release','gaia_index','legacy_index'],dtype=['U64']*3+['f8']*2+['U48']*3+['i8']*2)


def projection_models(image):
    """Read immutable detector WCS models once per identification stage."""
    models=[]
    with fits.open(image) as hd:
        for h in hd:
            if h.name!='SCI':continue
            w=detector_wcs(hd,h.header.get('EXTVER',1))
            models.append((h.header.get('EXTVER',1),h.shape[1],h.shape[0],w))
    return models


def projections(image,ra,dec,options,models=None):
    out=[]
    for version,nx,ny,w in projection_models(image) if models is None else models:
        try:x,y=w.all_world2pix([[ra,dec]],0)[0]
        except Exception:continue
        if not np.isfinite([x,y]).all():continue
        scale=np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600
        radius=min(options['search_radius_arcsec']/scale,options['max_search_radius_pix'])
        if -radius<=x<nx+radius and -radius<=y<ny+radius:
            out.append(dict(extver=version,predicted_x=float(x),predicted_y=float(y),radius_pix=float(radius),covered=bool(0<=x<nx and 0<=y<ny)))
    return out


def _local_raw(image,prediction,classification,options,detection,image_contexts=None):
    """Native Photutils local detection; predicted coordinates are never centroids."""
    import inspect
    from photutils.detection import DAOStarFinder
    from photutils.segmentation import detect_sources,deblend_sources,SourceCatalog
    from photutils.centroids import centroid_2dg
    v=prediction['extver'];radius=prediction['radius_pix'];xp=prediction['predicted_x'];yp=prediction['predicted_y']
    with fits.open(image) as hd:
        h=hd['SCI',v]
        image_contexts={} if image_contexts is None else image_contexts
        if v not in image_contexts:
            a=np.asarray(h.data,float);dq,base,mask=pixel_masks(hd,v,detection)
            err=np.asarray(image_data(hd['ERR',v]),float) if ('ERR',v) in hd else None
            image_contexts[v]=dict(data=a,dq=dq,base=base,mask=mask,error=err)
        cached=image_contexts[v];a=cached['data'];dq=cached['dq'];base=cached['base'];mask=cached['mask']
        x0=max(0,int(np.floor(xp-radius)));x1=min(a.shape[1],int(np.ceil(xp+radius+1)));y0=max(0,int(np.floor(yp-radius)));y1=min(a.shape[0],int(np.ceil(yp+radius+1)))
        if x1<=x0 or y1<=y0:return Table(),dict(reason='not_covered'),{}
        cut=a[y0:y1,x0:x1];bad=mask[y0:y1,x0:x1];_,sky,rms=sigma_clipped_stats(cut,mask=bad)
        if not np.isfinite(rms) or rms<=0:return Table(),dict(reason='insufficient_valid_pixels'),{}
        err=cached['error'] if cached['error'] is not None else np.full(a.shape,rms)
        context=dict(data=a,dq=dq,base=base,mask=mask,sky=float(sky),rms=float(rms),error=err)
        raw=[];segmaps={};count=0
        dao=DAOStarFinder(detection.get('snr',5)*rms,fwhm=detection.get('fwhm_pixels',2.5))(np.where(bad,sky,cut)-sky,mask=bad)
        if dao is not None:
            for r in dao:
                x=float(r['x_centroid'] if 'x_centroid' in dao.colnames else r['xcentroid']);y=float(r['y_centroid'] if 'y_centroid' in dao.colnames else r['ycentroid'])
                raw.append(dict(x=x+x0,y=y+y0,extver=v,flux=float(r['flux']),detection_method='local_DAO',segment_label=0,segmentation_key=0,local_sky=float(sky)))
        # Compact galaxies and sources not found by the point detector use segmentation.
        kernel=options['kernel_width_pix']/2.355;smooth=gaussian_filter(np.where(bad,0,cut-sky),kernel)
        area='n_pixels' if 'n_pixels' in inspect.signature(detect_sources).parameters else 'npixels'
        levels='n_levels' if 'n_levels' in inspect.signature(deblend_sources).parameters else 'nlevels'
        seg=detect_sources(smooth,detection.get('snr',5)*rms/max(1,2*np.sqrt(np.pi)*kernel),mask=bad,**{area:9})
        if seg is not None:
            seg=deblend_sources(smooth,seg,progress_bar=False,contrast=.005,**{area:9,levels:32})
            cat=SourceCatalog(np.where(bad,0,cut-sky),seg,error=err[y0:y1,x0:x1],mask=bad)
            segmaps[1]=dict(data=seg.data,x0=x0,y0=y0)
            for k in range(len(cat)):
                x=float((cat.x_centroid if hasattr(cat,'x_centroid') else cat.xcentroid)[k]);y=float((cat.y_centroid if hasattr(cat,'y_centroid') else cat.ycentroid)[k])
                major=float((cat.semimajor_axis if hasattr(cat,'semimajor_axis') else cat.semimajor_sigma)[k].value);minor=float((cat.semiminor_axis if hasattr(cat,'semiminor_axis') else cat.semiminor_sigma)[k].value)
                raw.append(dict(x=x+x0,y=y+y0,extver=v,flux=float(cat.segment_flux[k]),snr=float(cat.segment_flux[k]/max(float((cat.segment_flux_err if hasattr(cat,'segment_flux_err') else cat.segment_fluxerr)[k]),1e-12)),
                    width_pix=major,minor_width_pix=minor,detection_method='local_segmentation',segment_label=int(cat.labels[k]),segmentation_key=1,local_sky=float(sky)))
        if not raw:return Table(),dict(reason='covered_undetected',context=context),{}
        table=Table(rows=raw);table['source_id']=[f'local-{i}' for i in range(len(table))]
        segments={(v,1):segmaps[1]} if 1 in segmaps else {}
        accepted,audit=validate_centroids(image,table,detection,method='local',contexts={v:context},segmentations=segments,return_all=True)
        accepted['x_original']=np.asarray(accepted['x'],float);accepted['y_original']=np.asarray(accepted['y'],float);accepted['centroid_refined']=np.zeros(len(accepted),bool)
        accepted['detection_method']=np.asarray(accepted['detection_method']).astype('U64')
        # Gaussian refinement is allowed only for an already valid, point-like core.
        for i,r in enumerate(accepted):
            if not classification.startswith('gaia') and r['detection_method']!='local_DAO':continue
            if r['width_pix']>2:continue
            x=int(round(r['x']));y=int(round(r['y']));rr=5
            if min(x,y)<rr or x+rr>=a.shape[1] or y+rr>=a.shape[0]:continue
            sub=a[y-rr:y+rr+1,x-rr:x+rr+1];mb=mask[y-rr:y+rr+1,x-rr:x+rr+1]
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('error');cx,cy=centroid_2dg(sub-sky,error=err[y-rr:y+rr+1,x-rr:x+rr+1],mask=mb)
                xx=float(cx+x-rr);yy=float(cy+y-rr)
                if not np.isfinite([xx,yy]).all() or np.hypot(xx-r['x'],yy-r['y'])>max(1.,3*r['centroid_error_pix']):continue
                refined=accepted[i:i+1].copy();refined['x']=xx;refined['y']=yy
                valid=validate_centroids(image,refined,detection,method='local_Gaussian',contexts={v:context})
                if len(valid):
                    accepted['x'][i]=xx;accepted['y'][i]=yy;accepted['centroid_refined'][i]=True;accepted['detection_method'][i]='local_DAO+centroid_2dg'
            except Exception:pass
        return accepted,dict(audit=audit,context=context,segments=segments,local_sky=float(sky),local_rms=float(rms)),segmaps


def associate(image,candidate,options,detection,paths=None,label='identification',image_contexts=None):
    predicted=projections(image,float(candidate['ra']),float(candidate['dec']),options)
    base=dict(source_id=str(candidate['source_id']),exposure=image.name,classification=str(candidate['classification']),gaia_id=str(candidate['gaia_id']),legacy_id=str(candidate['legacy_id']),
        catalog_ra=float(candidate['ra']),catalog_dec=float(candidate['dec']),release=str(candidate['release']),status='not_covered',reason='',extver=-1,x=np.nan,y=np.nan,centroid_error_pix=np.nan,snr=np.nan,width_pix=np.nan,ellipticity=np.nan,minor_width_pix=np.nan,local_sky=np.nan,flux_rate=np.nan,predicted_x=np.nan,predicted_y=np.nan)
    if not predicted:return base
    base.update(nominal_detector_coverage=any(p['covered'] for p in predicted))
    base.update(predicted_x=predicted[0]['predicted_x'],predicted_y=predicted[0]['predicted_y'],extver=predicted[0]['extver'])
    valid=[];rejected=0;detected=0
    for p in predicted:
        cat,info,seg=_local_raw(image,p,str(candidate['classification']),options,detection,image_contexts)
        if 'audit' in info:
            audit=info['audit'];rejected+=int((~np.asarray(audit['centroid_valid'])).sum());detected+=len(audit)
            if paths:
                # Save full detector/QC evidence, including all rejected candidates.
                _,audit=validate_centroids(image,audit,dict(detection,diagnostic_overlay=False,max_diagnostic_cutouts=0),method='local',contexts={p['extver']:info['context']},segmentations=info['segments'],paths=paths,label=label+f'_sci{p["extver"]}',return_all=True)
        if image_contexts and p['extver'] in image_contexts and p['covered']:
            ctx=image_contexts[p['extver']];ix=int(round(p['predicted_x']));iy=int(round(p['predicted_y']))
            ix=int(np.clip(ix,0,ctx['data'].shape[1]-1));iy=int(np.clip(iy,0,ctx['data'].shape[0]-1))
            base['predicted_core_dq']=int(ctx['dq'][iy,ix]);base['predicted_core_masked']=bool(ctx['mask'][iy,ix])
        if paths:_window_plot(image,p,cat,seg,paths,label,plot=bool(candidate.get('_plot',True)))
        for row in cat:
            if np.hypot(row['x']-p['predicted_x'],row['y']-p['predicted_y'])>p['radius_pix']:continue
            # Merge DAO and segmented measurements of the same local object.
            existing=next((j for j,(r,pr) in enumerate(valid) if r['extver']==row['extver'] and np.hypot(r['x']-row['x'],r['y']-row['y'])<options['detector_merge_radius_pix']),None)
            if existing is None:valid.append((dict(zip(row.colnames,row)),p))
            elif row['centroid_error_pix']<valid[existing][0]['centroid_error_pix']:valid[existing]=(dict(zip(row.colnames,row)),p)
    base.update(detected_candidates=detected,rejected_candidates=rejected)
    if not valid:
        base.update(status='detected_rejected' if detected else 'covered_undetected',reason='No reliable HST centroid; expected catalog position is not a measurement');return base
    if len(valid)>1:
        from .stellar import stellar_options,fit_candidates,choose_profile,apply_selection,diagnostic
        if str(candidate['classification']).endswith(('point_candidate','stellar_candidate')) and stellar_options(options)['enabled']:
            audits=[]
            for v in sorted({int(row['extver']) for row,p in valid}):
                cat=Table(rows=[r for r,p in valid if int(r['extver'])==v]);cat['hst_source_id']=[f'LOCALSCI{v}_{i}' for i in range(len(cat))]
                good,audit=fit_candidates(image,cat,v,detection,options,image_contexts);audits.append(audit)
            allfits=[r for a in audits for r in a['fits']];good=[r for r in allfits if r['status']=='credible_stellar_fit']
            base.update(status='ambiguous_stellar_profiles' if len(good)>1 else 'no_reliable_stellar_profile',reason=f'{len(good)} credible stellar profiles from {len(valid)} local candidates',
                psf_audit=dict(fits=allfits,raw_candidates=len(valid),duplicates=sum(a.get('duplicates',0) for a in audits),rejected=sum(a.get('rejected',0) for a in audits)))
            apply_selection(base,choose_profile(good,None,options));diagnostic(image,base,paths,label);return base
        base.update(status='ambiguous_association',reason=f'{len(valid)} distinct valid HST detections in bounded search region');return base
    row,p=valid[0];h=fits.getheader(image)
    base.update(status='valid_detection',reason='',extver=int(row['extver']),x=float(row['x']),y=float(row['y']),centroid_error_pix=float(row['centroid_error_pix']),snr=float(row['snr']),width_pix=float(row['width_pix']),ellipticity=float(row['ellipticity']),
        minor_width_pix=float(row['minor_width_pix']),local_sky=float(row.get('local_sky',0)),flux_rate=float(row['flux'])/float(h.get('EXPTIME',1)),predicted_x=p['predicted_x'],predicted_y=p['predicted_y'],centroid_method=str(row['detection_method']),centroid_refined=bool(row.get('centroid_refined',False)),x_original=float(row.get('x_original',row['x'])),y_original=float(row.get('y_original',row['y'])),centroid_valid=True)
    return base


def _window_plot(image,p,cat,segments,paths,label,plot=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval
    v=p['extver'];radius=min(int(p['radius_pix']),30)
    with fits.open(image) as hd:
        h=hd['SCI',v];x=int(round(p['predicted_x']));y=int(round(p['predicted_y']))
        sx=slice(max(0,x-radius),min(h.shape[1],x+radius+1));sy=slice(max(0,y-radius),min(h.shape[0],y+radius+1))
        if sx.stop<=sx.start or sy.stop<=sy.start:return
        sci=h.data[sy,sx];dq=hd['DQ',v].data[sy,sx] if ('DQ',v) in hd else np.zeros(sci.shape)
        from .centroids import _segment_patch
        seg=_segment_patch(segments[1],sy,sx) if 1 in segments else np.zeros(sci.shape)
        primary=fits.PrimaryHDU();primary.header['SEGAVAIL']=bool(1 in segments);primary.header['SRCFILE']=Path(image).name
        primary.header['SCIEXT']=v;primary.header['XORIGIN']=sx.start;primary.header['YORIGIN']=sy.start
        stamps=[primary]
        planes=[('SCI',sci),('DQ',dq)]+([('SEG',seg)] if 1 in segments else [])
        for name,data in planes:
            h=fits.ImageHDU(data,name=name);h.header['XPRED']=p['predicted_x']-sx.start;h.header['YPRED']=p['predicted_y']-sy.start
            if name=='SCI':h.header['BUNIT']=hd['SCI',v].header.get('BUNIT','UNKNOWN')
            stamps.append(h)
        fits.HDUList(stamps).writeto(paths.file('intermediate','alignment',label+f'_sci{v}_prior','fits'),overwrite=True,checksum=True)
        if not plot:return
        fig,axes=plt.subplots(1,3,figsize=(10,3))
        finite=sci[np.isfinite(sci)];lo,hi=ZScaleInterval().get_limits(finite) if finite.size else (0,1)
        for ax,array,title in zip(axes,(sci,dq!=0,seg),('SCI / prior cross and HST detections','Input DQ','Segmentation')):
            ax.imshow(array,origin='lower',cmap='gray',vmin=lo if title.startswith('SCI') else None,vmax=hi if title.startswith('SCI') else None)
            ax.axvline(p['predicted_x']-sx.start,color='red',alpha=.5);ax.axhline(p['predicted_y']-sy.start,color='red',alpha=.5);ax.set_title(title,fontsize=8)
            if len(cat):ax.scatter(cat['x']-sx.start,cat['y']-sy.start,facecolors='none',edgecolors='lime',s=30)
        if 1 not in segments:
            axes[2].clear();axes[2].axis('off');axes[2].text(.5,.5,'Full-field detection: see detector QC',ha='center',va='center',transform=axes[2].transAxes,fontsize=9)
        fig.tight_layout();fig.savefig(paths.file('diagnostics','alignment',label+f'_sci{v}_prior','png'),dpi=100);plt.close(fig)


def identify(inputs,paths,astrometry):
    options=astrometry['source_identification'];cache=paths.cache()
    reportfile=paths.file('intermediate','alignment','source_identification');commonfile=paths.file('intermediate','alignment','common_sources','ecsv')
    paths.area('intermediate','alignment',create=True);paths.area('diagnostics','alignment',create=True)
    basekey=identification_key(cache,inputs,astrometry)
    from .full_field import detect_full_field
    full_catalogs,detection_report=detect_full_field(inputs,paths,astrometry['relative'])
    if not astrometry['catalogs']['refresh'] and astrometry['catalogs']['cache'] and cache.valid('uvis_identification',basekey):
        report=json.loads(reportfile.read_text());report['full_field']=detection_report
        save(reportfile,report)
        receipt=json.loads(cache.receipt('uvis_identification').read_text())
        cache.record('uvis_identification',basekey,[Path(o['path']) for o in receipt['outputs']])
        return dict(report,catalogs=[Table.read(p,format='ascii.ecsv') for p in report['hst_catalogs']])
    if not options['enabled']:return dict(status='disabled',catalogs=[],full_field=detection_report,external_positions_used_for_relative_fit=False)
    cone=footprint(inputs,options['margin_arcsec'])
    if cone['radius']>astrometry['catalogs']['max_radius_deg']:raise ValueError('Source-identification footprint exceeds configured maximum query region')
    tables={};queries=[];mean_epoch=float(np.mean([epoch(p) for p in inputs]))
    for name in options['catalogs']:
        q=dict(catalog=name,status='pending')
        try:
            provider=astrometry['catalogs'][name];quality=dict(provider['quality'])
            kwargs=dict(quality=quality,cache_dir=paths.catalog_directory(),refresh=(astrometry['catalogs']['refresh'] or not astrometry['catalogs']['cache']),offline=astrometry['catalogs']['offline'],timeout=options['timeout'],retries=options['retries'])
            if name=='gaia_dr3':quality['apply_proper_motion']=provider['apply_proper_motion']
            else:quality.update(types=options['legacy_types'],max_shape_arcsec=options['max_legacy_shape_arcsec'],allow_masked_stellar_priors=options.get('allow_masked_stellar_priors',True))
            reused=None
            if astrometry['catalogs']['cache'] and not astrometry['catalogs']['refresh']:
                try:
                    reuse=cached_footprint(name,cone['ra'],cone['dec'],cone['radius'],cache_dir=paths.catalog_directory(),quality=quality,observation_epoch=mean_epoch)
                    prior=reuse.meta['original_cone']
                    delta=SkyCoord(cone['ra']*u.deg,cone['dec']*u.deg).separation(SkyCoord(prior[0]*u.deg,prior[1]*u.deg)).deg
                    required=cone['radius']-(options['margin_arcsec']-options['search_radius_arcsec'])/3600
                    if delta+required<=prior[2]:reused=reuse;reused.meta['complete_flc_coverage']=True
                except Exception:pass
            if reused is not None:tables[name]=reused;q['status']='cached_covering_flc'
            elif name=='gaia_dr3':tables[name]=query_gaia_dr3(cone['ra'],cone['dec'],cone['radius'],mean_epoch,**kwargs)
            else:tables[name]=query_legacy_dr10(cone['ra'],cone['dec'],cone['radius'],tap_url=provider['tap_url'],**kwargs)
            q.update(status=q['status'] if q['status']=='cached_covering_flc' else ('available' if len(tables[name]) else 'no_selected_sources'),count=len(tables[name]),metadata=dict(tables[name].meta))
        except Exception as exc:
            q.update(status='query_failed' if isinstance(exc,CatalogQueryError) else 'catalog_processing_failed',failure_stage='service_request' if isinstance(exc,CatalogQueryError) else 'response_normalization',error=f'{type(exc).__name__}: {exc}')
            try:
                tables[name]=cached_footprint(name,cone['ra'],cone['dec'],cone['radius'],cache_dir=paths.catalog_directory(),quality=quality,observation_epoch=mean_epoch)
                prior=tables[name].meta['original_cone']
                delta=SkyCoord(cone['ra']*u.deg,cone['dec']*u.deg).separation(SkyCoord(prior[0]*u.deg,prior[1]*u.deg)).deg
                science_radius=cone['radius']-(options['margin_arcsec']-options['search_radius_arcsec'])/3600
                tables[name].meta['complete_flc_coverage']=bool(delta+science_radius<=prior[2])
                q.update(status='cached_after_query_failure',count=len(tables[name]),metadata=dict(tables[name].meta))
            except Exception as cache_error:q['cache_reuse_error']=str(cache_error)

        queries.append(q)
        print('Source identification:',name,q['status'],q.get('count',q.get('error','')),flush=True)
    empty=Table();gaia=tables.get('gaia_dr3',empty);legacy=tables.get('legacy_dr10',empty)
    candidates=merge_candidates(gaia,legacy,options);candidates.meta['purpose']='Identification priors only; never relative-fitting coordinates'
    candidates.write(paths.file('intermediate','alignment','identification_candidates','ecsv'),format='ascii.ecsv',overwrite=True)
    models={str(p):projection_models(p) for p in inputs}
    records=[]
    for image_index,image in enumerate(inputs):
        image_contexts={}
        ge=gaia_at_epoch(gaia,epoch(image)) if len(gaia) else gaia;le=legacy_at_epoch(legacy,epoch(image)) if len(legacy) else legacy
        for i,c in enumerate(candidates):
            c=dict(zip(c.colnames,c));position=ge[c['gaia_index']] if c['gaia_index']>=0 else le[c['legacy_index']]
            c.update(ra=float(position['ra']),dec=float(position['dec']),_plot=(i<options.get('max_window_plots',8) or bool(c.get('identification_prior_only',False)) or bool(c['gaia_id'])))
            result=associate_full_field(image,c,full_catalogs[image_index],options,paths,label=image.stem+'_identity_'+str(i+1),models=models[str(image)],detection=astrometry['relative']['detection'],contexts=image_contexts)
            records.append(result)
    from .stellar import install_measurement,profile_figure,resolve_profile_group
    stellar_summary=[]
    for candidate in candidates:
        group=[r for r in records if r['source_id']==str(candidate['source_id'])]
        resolve_profile_group(group,options)
        for r in group:
            if r.get('psf_selected'):
                i=next(i for i,p in enumerate(inputs) if p.name==r['exposure']);install_measurement(full_catalogs[i],r)
            if r.get('psf_audit'):
                stellar_summary.append(dict(exposure=r['exposure'],source_id=r['source_id'],status=r['status'],extver=r['extver'],predicted_x=r['predicted_x'],predicted_y=r['predicted_y'],
                    raw_candidates=r['psf_audit']['raw_candidates'],duplicates=r['psf_audit'].get('duplicates',0),rejected=r['psf_audit'].get('rejected',0),
                    selected=r.get('psf_selected'),reason=r['reason'],profile_consensus=r.get('profile_consensus'),audit_file=r.get('psf_audit_file')))
    for r in records:
        selected=r.pop('psf_selected',None);r.pop('psf_audit',None)
        r['centroid_method']=selected.get('detection_method','PSFPhotometry') if selected else 'full_field_detector'
        r['psf_fwhm_pix']=selected['fwhm_pix'] if selected else np.nan
        r.pop('profile_consensus',None)
        r['psf_audit_file']=r.get('psf_audit_file','')
    associations=Table(rows=records) if records else Table(names=['source_id','exposure','status'],dtype=['U64','U64','U64'])
    associations['status']=np.asarray(associations['status']).astype('U64')
    verified=[]
    for sid in candidates['source_id']:
        idx=np.flatnonzero((associations['source_id']==sid)&(associations['status']=='valid_detection'))
        if len(idx)<2:continue
        widths=np.asarray(associations['width_pix'][idx]);flux=np.asarray(associations['flux_rate'][idx]);ellip=np.asarray(associations['ellipticity'][idx])
        stable=widths.max()/max(widths.min(),1e-6)<=options['max_width_ratio'] and np.ptp(ellip)<=options['max_ellipticity_difference'] and flux.min()>0 and flux.max()/flux.min()<=options['max_flux_ratio']
        if stable:associations['status'][idx]='confirmed_common';verified.append(str(sid))
        else:associations['status'][idx]='inconsistent_morphology'
    # One HST centroid cannot identify two distinct external objects.
    for image in inputs:
        idx=np.flatnonzero((associations['exposure']==image.name)&np.isin(associations['status'],['valid_detection','confirmed_common']))
        for k,i in enumerate(idx):
            for j in idx[k+1:]:
                if associations['extver'][i]==associations['extver'][j] and np.hypot(associations['x'][i]-associations['x'][j],associations['y'][i]-associations['y'][j])<options['detector_merge_radius_pix']:
                    associations['status'][[i,j]]='ambiguous_duplicate_association'
    verified=sorted({str(sid) for sid in associations['source_id'] if sum((associations['source_id']==sid)&(associations['status']=='confirmed_common'))>=2})
    for i,r in enumerate(associations):
        if r['status']=='confirmed_common' and r['source_id'] not in verified:associations['status'][i]='valid_not_common'
    from .common_sources import build_common,resolve_repeated_counterparts
    catalogs,hst_common,common_report=build_common(inputs,full_catalogs,associations,astrometry['relative'],association_options=options)
    from .source_recovery import recover_missing
    centroid_recovery=recover_missing(inputs,full_catalogs,hst_common,astrometry['relative']['detection'],astrometry['relative'])
    if any(r['status']=='measured_candidate' for r in centroid_recovery):
        catalogs,hst_common,common_report=build_common(inputs,full_catalogs,associations,astrometry['relative'],association_options=options)
    common_report['centroid_recovery']=centroid_recovery
    repeated_recovery=resolve_repeated_counterparts(inputs,full_catalogs,hst_common,associations,options)
    hst_common.write(paths.file('intermediate','alignment','hst_common_sources','ecsv'),format='ascii.ecsv',overwrite=True)
    files=[]
    for image,cat in zip(inputs,catalogs):
        out=paths.file('intermediate','alignment',image.stem+'_identified_sources','ecsv');cat.write(out,format='ascii.ecsv',overwrite=True);files.append(str(out))
    associations.write(commonfile,format='ascii.ecsv',overwrite=True)
    for s in stellar_summary:
        row=next(r for r in associations if str(r['source_id'])==s['source_id'] and str(r['exposure'])==s['exposure'])
        s['status']=str(row['status'])
        if str(row['centroid_method'])=='HST_repeat_association':
            s['reason']=str(row['reason'])
            s['HST_repeat_recovery']=dict(source_id=str(row['hst_source_id']),x=float(row['x']),y=float(row['y']),extver=int(row['extver']),centroid_error_pix=float(row['centroid_error_pix']))
        if s['status'] not in ('valid_detection','confirmed_common'):s['reason']=s['status']+'; '+s['reason']
    for sid in sorted({s['source_id'] for s in stellar_summary}):
        group=[]
        for s in stellar_summary:
            if s['source_id']!=sid:continue
            r=dict(s,psf_audit=json.loads(Path(s['audit_file']).read_text()))
            if s['status'] in ('valid_detection','confirmed_common'):r['psf_selected']=s.get('selected')
            if s.get('HST_repeat_recovery'):r['selected_hst']=s['HST_repeat_recovery']
            group.append(r)
        profile_figure(inputs,group,paths,sid)
    verified=sorted({str(r['source_id']) for r in associations if r['status']=='confirmed_common' and sum((associations['source_id']==r['source_id'])&(associations['status']=='confirmed_common'))>=2})
    by_exposure={p.name:{str(status):int(sum((associations['exposure']==p.name)&(associations['status']==status))) for status in set(associations['status'])} for p in inputs}
    kinds={str(r['source_id']):str(r['classification']) for r in candidates}
    common_counts=dict(gaia_stellar_candidates=sum(kinds[s]=='gaia_stellar_candidate' for s in verified),legacy_point_candidates=sum(kinds[s]=='legacy_point_candidate' for s in verified),compact_galaxy_candidates=sum(kinds[s]=='compact_galaxy_candidate' for s in verified))
    actual=[c for c in candidates if any(any(p['covered'] for p in projections(image,float(c['ra']),float(c['dec']),options,models=models[str(image)])) for image in inputs)]
    inventory=dict(gaia_in_footprints=sum(bool(c['gaia_id']) for c in actual),legacy_psf_in_footprints=sum(bool(c['legacy_id']) and c['legacy_type']=='PSF' for c in actual),legacy_compact_in_footprints=sum(bool(c['legacy_id']) and c['legacy_type']!='PSF' for c in actual))
    inventory['gaia_in_footprints']=inventory['gaia_in_footprints'] if 'gaia_dr3' in tables else None
    inventory['legacy_psf_in_footprints']=inventory['legacy_psf_in_footprints'] if 'legacy_dr10' in tables else None
    inventory['legacy_compact_in_footprints']=inventory['legacy_compact_in_footprints'] if 'legacy_dr10' in tables else None
    summary=dict(HST_repeat_recovery=repeated_recovery,stellar_resolution=stellar_summary,full_field=detection_report,hst_association=common_report,per_exposure=by_exposure,common_source_counts=common_counts,confirmation_definition='Unique image-supported counterpart, independent valid HST centroids in >=2 exposures, stable morphology/flux; relative fit still requires validation',measurement_frame='Read-only calibrated FLC detector pixels with local backgrounds',inventory=inventory,status='identified' if verified else 'no_verified_common_sources',queries=queries,footprint=cone,verified_common_ids=verified,hst_catalogs=files,
        external_positions_used_for_relative_fit=False,relative_coordinate_definition='Independent measured HST x/y transformed through HST distortion-aware WCS; IDs only constrain associations',warnings=[] if tables else ['No external catalog available: retaining source-only relative-registration fallback'])
    save(reportfile,summary)
    outputs=[reportfile,commonfile,paths.file('intermediate','alignment','hst_common_sources','ecsv')]+[Path(p) for p in files]+[Path(p) for p in detection_report['catalogs']]+[p for p in paths._files if p.name.endswith('_quality.ecsv') and p.exists()]
    outputs+=[Path(t.meta['cache_file']) for t in tables.values() if t.meta.get('cache_file')]
    cache.record('uvis_identification',basekey,list(dict.fromkeys(outputs)))
    return dict(summary,catalogs=catalogs)


def associate_full_field(image,candidate,catalog,options,paths=None,label='identification',models=None,detection=None,contexts=None):
    """Identify an independently measured full-field detection; never force nearest."""
    predicted=projections(image,float(candidate['ra']),float(candidate['dec']),options,models=models)
    result=dict(source_id=str(candidate['source_id']),exposure=image.name,gaia_id=str(candidate['gaia_id']),legacy_id=str(candidate['legacy_id']),classification=str(candidate['classification']),
        catalog_ra=float(candidate['ra']),catalog_dec=float(candidate['dec']),release=str(candidate['release']),status='not_covered',reason='',extver=-1,
        predicted_x=np.nan,predicted_y=np.nan,x=np.nan,y=np.nan,snr=np.nan,width_pix=np.nan,ellipticity=np.nan,centroid_error_pix=np.nan,flux_rate=np.nan,hst_source_id='')
    if not predicted:return result
    possible=[]
    psf_records=[]
    with fits.open(image) as hd:
        for p in predicted:
            v=p['extver'];x=p['predicted_x'];y=p['predicted_y'];result.update(extver=v,predicted_x=x,predicted_y=y)
            h=hd['SCI',v];dq=image_data(hd['DQ',v]) if ('DQ',v) in hd else np.zeros(h.shape,'u2');r=p['radius_pix']
            sx=slice(max(0,int(x-r)),min(h.shape[1],int(x+r+1)));sy=slice(max(0,int(y-r)),min(h.shape[0],int(y+r+1)))
            use=(catalog['extver']==v)&(np.hypot(catalog['x']-x,catalog['y']-y)<=r)
            possible.extend(list(catalog[use]))
            if str(candidate['classification']).endswith(('point_candidate','stellar_candidate')):
                from .stellar import fit_candidates,stellar_options
                if stellar_options(options)['enabled']:
                    from .centroids import QUALITY_DEFAULTS
                    measured,audit=fit_candidates(image,catalog[use],v,detection or dict(QUALITY_DEFAULTS,snr=5,fwhm_pixels=2.5),options,contexts)
                    psf_records.append(audit)
            if paths:_window_plot(image,p,catalog[use],{},paths,label,plot=bool(candidate.get('_plot',False)))
        if psf_records:
            from .stellar import apply_selection,choose_profile,diagnostic
            fits_records=[q for a in psf_records for q in a['fits']]
            audit=dict(status='candidate_psf_assessment',raw_candidates=len(possible),duplicates=sum(a.get('duplicates',0) for a in psf_records),rejected=sum(a.get('rejected',0) for a in psf_records),fits=fits_records,
                model='; '.join(dict.fromkeys(a.get('model','Native HST profile assessment') for a in psf_records)))
            good=[q for q in fits_records if q['status']=='credible_stellar_fit']
            result.update(status='ambiguous_stellar_profiles' if len(good)>1 else 'no_reliable_stellar_profile',reason=f'{len(good)} credible stellar profiles from {len(possible)} HST candidates',psf_audit=audit)
            apply_selection(result,choose_profile(good,None,options));diagnostic(image,result,paths,label)
            return result
        if len(possible)==1:
            row=possible[0];result.update(status='valid_detection',extver=int(row['extver']),x=float(row['x']),y=float(row['y']),snr=float(row['snr']),width_pix=float(row['width_pix']),ellipticity=float(row['ellipticity']),
                centroid_error_pix=float(row['centroid_error_pix']),flux_rate=float(row['flux'])/float(hd[0].header.get('EXPTIME',1)),hst_source_id=str(row['hst_source_id']))
        elif len(possible)>1:result.update(status='ambiguous_association',reason=f'{len(possible)} vetted full-field detections; no nearest-object assumption')
        else:result.update(status='covered_undetected',reason='No quality-approved full-field measurement; inspect detector QC tables for rejected candidates')
    return result


def identification_key(cache,inputs,astrometry):
    options=astrometry['source_identification']
    code=[Path(__file__),Path(__file__).with_name('stellar.py'),Path(__file__).with_name('source_recovery.py'),Path(__file__).with_name('centroids.py'),Path(__file__).with_name('catalogs.py'),Path(__file__).with_name('full_field.py'),Path(__file__).with_name('common_sources.py'),Path(__file__).with_name('wcs.py'),Path(__file__).with_name('relative_fit.py'),Path(__file__).with_name('spatial_residuals.py'),Path(__file__).parents[1]/'dq.py']
    return cache.key(inputs+code,dict(options=options,detection=astrometry['relative']['detection'],catalog_options=astrometry['catalogs'],versions=versions()))
