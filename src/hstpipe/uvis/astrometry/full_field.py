"""Cached full-field detector catalogs, independent of external reference data."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table,vstack
from .sources import sources
from .centroids import QUALITY_VERSION
from ...common.cache import StageCache
from ...common.target_runner import save
from ...common.provenance import versions

def merge_detectors(catalogs,radius=2.):
    """Deduplicate measured positions, preferring the better centroid precision."""
    nonempty=[c for c in catalogs if len(c)]
    if not nonempty:
        c=catalogs[0].copy();c['hst_source_id']=np.array([],dtype='U128');return c
    # Metadata are kept explicitly per detector in the stage report.
    parts=[]
    for c in nonempty:
        c=c.copy();c.meta={};parts.append(c)
    allrows=vstack(parts,join_type='outer')
    order=np.argsort(np.asarray(allrows['centroid_error_pix'],float));chosen=[]
    for k in order:
        row=allrows[k]
        if any(int(row['extver'])==int(allrows[j]['extver']) and np.hypot(row['x']-allrows[j]['x'],row['y']-allrows[j]['y'])<radius for j in chosen):continue
        chosen.append(int(k))
    chosen.sort(key=lambda k:(int(allrows[k]['extver']),float(allrows[k]['x']),float(allrows[k]['y'])))
    out=allrows[chosen];out['detector_source_id']=np.asarray(out['source_id']).astype('U128')
    out['hst_source_id']=[f"{r['exposure_rootname']}_sci{int(r['extver'])}_{k:05d}" for k,r in enumerate(out)]
    out.meta['centroid_quality']=dict(version=QUALITY_VERSION,accepted=len(out),purpose='Deduplicated full-field HST measurements')
    return out


def detect_full_field(inputs,paths,options):
    paths.area('intermediate','alignment',create=True);paths.area('diagnostics','alignment',create=True)
    cache=paths.cache();base=Path(__file__).parent
    key=detection_key(cache,inputs,options)
    reportfile=paths.file('intermediate','alignment','full_field_sources')
    if cache.valid('uvis_sources',key):
        import json
        r=json.loads(reportfile.read_text());catalogs=[Table.read(p,format='ascii.ecsv') for p in r['catalogs']]
        from ..chips import source_products
        for image,cat,record in zip(inputs,catalogs,r['exposures']):record['chips']=source_products(image,cat,record['methods'],paths)
        save(reportfile,r)
        receipt=__import__('json').loads(cache.receipt('uvis_sources').read_text())
        cache.record('uvis_sources',key,[Path(o['path']) for o in receipt['outputs']])
        return catalogs,r
    catalogs=[];records=[];outputs=[]
    for image in inputs:
        method_cats=[];methods=[]
        for extended,width in [(False,0)]+[(True,float(w)) for w in options['kernel_widths']]:
            label=image.stem+('_full_DAO' if not extended else f'_full_seg_w{width}')
            cat=sources(image,options['detection'],extended,width or 3.5,paths=paths,label=label)
            method_cats.append(cat);methods.append(dict(method='DAOStarFinder' if not extended else 'Photutils segmentation',
                kernel_fwhm_pix=width if extended else options['detection']['fwhm_pixels'],threshold_snr=options['detection']['snr'],
                **cat.meta.get('centroid_quality',{}),quality_file=cat.meta.get('quality_file'),chips=cat.meta.get('detection_chips',[])))
        cat=merge_detectors(method_cats,options.get('dedup_radius_pix',2.))
        output=paths.file('intermediate','alignment',image.stem+'_full_sources','ecsv');cat.write(output,format='ascii.ecsv',overwrite=True)
        catalogs.append(cat);outputs.append(output)
        with fits.open(image) as hd:units={str(h.header.get('EXTVER',1)):h.header.get('BUNIT') for h in hd if h.name=='SCI'}
        from ..chips import source_products
        chips=source_products(image,cat,methods,paths)
        records.append(dict(exposure=image.name,units=units,methods=methods,unique_accepted=len(cat),chips=chips))
    report=dict(status='complete',catalogs=[str(p) for p in outputs],exposures=records,
        scope='Entire SCI extensions on native detector pixels; no rectification/resampling',source_identification_required=False)
    save(reportfile,report)
    outputs+=[Path(m['quality_file']) for r in records for m in r['methods'] if m.get('quality_file')]
    cache.record('uvis_sources',key,outputs+[reportfile])
    return catalogs,report


def detection_key(cache,inputs,options):
    base=Path(__file__).parent
    deps=[base/n for n in ('full_field.py','sources.py','matching.py','centroids.py','wcs.py')]+[base.parent/'dq.py']
    return cache.key(inputs+deps,dict(detection=options['detection'],widths=options['kernel_widths'],versions=versions()))
