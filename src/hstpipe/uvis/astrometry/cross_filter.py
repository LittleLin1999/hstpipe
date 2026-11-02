"""Internal provisional-image registration. No external catalog dependencies."""
import json,tempfile
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from stwcs.wcsutil import altwcs
from .sources import sources,project
from .relative_fit import match,quality
from .relative_engine import native_catalog_fit
from .wcs import signature
from ...common.paths import ReductionPaths
from ...common.cache import StageCache
from ...common.provenance import identity
from ...common.target_runner import save,target_lock


def internal_cross_filter(directory,selected,options):
    """Measure cross-filter transforms on copies; retain original provisional WCS.

    Accepted internal relationships are saved, not prematurely propagated to FLCs.
    The absolute stage validates the combined transfer against its solved anchor.
    """
    directory=Path(directory);statefile=directory/'manifest.json'
    with target_lock(directory):
        state=json.loads(statefile.read_text());available={};outputs={}
        for key,r in state['filters'].items():
            a=r.get('scientific_assessment',{});rel=a.get('relative_astrometry',a.get('alignment',{}))
            p=ReductionPaths(directory.parent.parent,directory.parent.name,'uvis',r['filter'],r.get('setup')).file('intermediate','drizzle','provisional','fits')
            if rel.get('status') in ('validated','single_exposure') and p.exists():available[key]=(p,r)
        # Normalize archived WCS only on temporary copies; earlier absolute zero points
        # must not make internal matching depend on an external solution.
        with tempfile.TemporaryDirectory(prefix='hstpipe-cross-filter-') as tmp:
            normalized={};catalogs={}
            for key,(p,r) in available.items():
                copy=Path(tmp)/p.name
                paths=ReductionPaths(directory.parent.parent,directory.parent.name,'uvis',r['filter'],r.get('setup'))
                paths.copy(p,copy)
                with fits.open(copy) as hd:v=hd['SCI'].header.get('EXTVER',1);original=hd['SCI'].header.get('HPORIGK')
                if original:altwcs.restoreWCS(str(copy),ext=[('SCI',v)],wcskey=original)
                normalized[key]=copy;catalogs[key]=sources(copy,r.get('configuration',{}).get('options',{}).get('astrometry',options)['relative']['detection'],True,6,paths=paths,label='cross_filter_sources')
            requested=options['cross_filter']['reference_filter']
            pool=[k for k,(_,r) in available.items() if requested=='auto' or r['filter']==requested]
            # Source quality, not wavelength, controls the initial internal anchor.
            score=lambda k:sum(1/max(float(e),.08)**2 for e in catalogs[k]['centroid_error_pix'])
            anchor=max(pool,key=score) if pool else None
            for key in selected:
                r=state['filters'][key];paths=ReductionPaths(directory.parent.parent,directory.parent.name,'uvis',r['filter'],r.get('setup'))
                effective=r.get('configuration',{}).get('options',{}).get('astrometry',options)
                record=dict(status='skipped',anchor_group=anchor,external_catalog_used=False)
                if key not in available:record.update(reason='Within-filter relative registration is unvalidated')
                elif not options['cross_filter']['enabled']:record.update(reason='Cross-filter registration disabled')
                elif len(available)<2:record.update(status='not_applicable',reason='Only one validated provisional image')
                elif anchor is None:record.update(reason='Requested internal reference filter unavailable')
                elif key==anchor:record.update(status='anchor',source_count=len(catalogs[key]))
                else:
                    cache=paths.cache()
                    deps=[Path(__file__),Path(__file__).with_name('sources.py'),Path(__file__).with_name('relative_fit.py'),Path(__file__).with_name('relative_engine.py'),Path(__file__).with_name('wcs.py'),Path(__file__).with_name('centroids.py')]
                    params=dict(target=signature(available[key][0]),anchor=signature(available[anchor][0]),options=options['cross_filter'],relative=effective['relative'],detection=effective['relative']['detection'])
                    cachekey=cache.key(deps,params);output=paths.file('intermediate','alignment','cross_filter')
                    if cache.valid('uvis_cross_filter',cachekey):record=json.loads(output.read_text())
                    else:
                        try:
                            ref=normalized[anchor];image=normalized[key];a=catalogs[anchor];b=catalogs[key]
                            with fits.open(ref) as hd:tangent=WCS(hd['SCI'].header,hd)
                            x=project(a,ref,tangent);y=project(b,image,tangent)
                            ri,ii=match(x,y,effective['relative'],True)
                            with fits.open(image) as hd:target_scale=np.sqrt(abs(np.linalg.det(WCS(hd['SCI'].header,hd).pixel_scale_matrix)))
                            ref_scale=np.sqrt(abs(np.linalg.det(tangent.pixel_scale_matrix)))
                            error=np.hypot(a['centroid_error_pix'][ri],b['centroid_error_pix'][ii]*target_scale/ref_scale)
                            before=quality(x[ri],y[ii],error,effective['relative'])
                            record.update(before=before,anchor_signature=params['anchor'],target_signature=params['target'],anchor_filter=available[anchor][1]['filter'])
                            pairs=Table(dict(reference_index=ri,source_index=ii,dx_pix=(x[ri]-y[ii])[:,0],dy_pix=(x[ri]-y[ii])[:,1],centroid_error_pix=error))
                            pairs.write(paths.file('intermediate','alignment','cross_filter_matches','ecsv'),format='ascii.ecsv',overwrite=True)
                            if before['status'] in ('validated','refinement_supported'):
                                fits_report=native_catalog_fit([image],[b],effective['relative'],reference_image=ref,reference_catalog=a)
                                after=project(b,image,tangent)
                                validation=quality(x[ri],after[ii],error,effective['relative'])
                                record.update(status='validated' if validation['status']=='validated' else 'unvalidated',validation=validation,transformations=fits_report,
                                    caveat='Cross-filter morphology enters the centroid error floor; fit-star scatter is not independent absolute accuracy')
                            else:record.update(status='unvalidated',reason='No reliable common-source transformation')
                        except Exception as exc:record.update(status='unvalidated',error=f'{type(exc).__name__}: {exc}')
                        save(output,record);cache.record('uvis_cross_filter',cachekey,[p for p in paths._files if p.is_file() and p.parent!=Path(tmp)])
                output=paths.file('intermediate','alignment','cross_filter')
                # Do not rewrite a valid receipt's identical output.
                if not output.exists() or json.loads(output.read_text())!=record:save(output,record)
                outputs[key]=record
                qa=paths.file('diagnostics','drizzle','qa_summary');assessment=r.get('scientific_assessment',{})
                assessment['cross_filter_alignment']=record
                assessment.setdefault('stages',{})['cross_filter_alignment']=dict(status=record['status'],inputs=[str(available[k][0]) for k in (key,anchor) if k in available],outputs=[str(output)],external_catalog_used=False)
                r['scientific_assessment']=assessment
                if qa.exists():save(qa,assessment)
                artifacts={a['path']:a for a in r.get('artifacts',[]) if Path(a['path']).exists()}
                for p in paths._files|{qa}:
                    if p.is_file() and not p.is_relative_to(Path(tmp)):artifacts[str(p.resolve())]=identity(p)
                r['artifacts']=list(artifacts.values())
                r['intermediate_files']=[a for a in artifacts.values() if 'intermediate' in Path(a['path']).relative_to(directory).parts]
                state['filters'][key]=r
        save(statefile,state)
    return outputs
