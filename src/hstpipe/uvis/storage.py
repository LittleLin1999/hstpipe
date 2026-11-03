"""UVIS checkpoint identities and serialization; no scientific computations."""
from pathlib import Path
from astropy.table import Table
from ..common.provenance import versions
from .amplifiers import reference_path

def code_inputs():
    # Presentation changes regenerate the report on a normal cache hit. They
    # must not invalidate scientific reduction receipts by themselves.
    return [Path(__file__).parents[1]/'common/calibration.py']+sorted(p for p in Path(__file__).parent.rglob('*.py')
                  if 'reporting' not in p.relative_to(Path(__file__).parent).parts
                  and p.name!='reporting.py' and not p.name.startswith('report_') and not p.name.startswith('._'))

def references(options,rows):
    refs=[]
    if options['alignment'].get('refcat'):refs.append(Path(options['alignment']['refcat']))
    from ..common.calibration import reference_requests,reference_location
    import os
    directory=Path(options.get('calibration',{}).get('directory') or os.environ.get('iref') or '.')
    for item in reference_requests([e.path for e in rows],bool(options['background']['amplifier_correction'] or options['background'].get('amplifier_assess'))):
        try:path=reference_location(item['header_value'],directory)
        except ValueError:continue  # Preparation records the actionable error.
        if path.is_file():refs.append(path)
    return list(dict.fromkeys(refs))

def stage_keys(cache,rows,options,copies,identification=None):
    package=Path(__file__).parent
    background_code=[package/n for n in ('background.py','amplifiers.py','fits.py','inventory.py','dq.py')]
    alignment_code=[package/n for n in ('alignment.py','fits.py','inventory.py','dq.py')]+list((package/'astrometry').glob('relative*.py'))+[package/'astrometry/sources.py',package/'astrometry/matching.py',package/'astrometry/spatial_residuals.py',package/'astrometry/centroids.py']
    cr_code=[package/n for n in ('drizzle.py','qa.py')]
    raw=[e.path for e in rows]
    refs=references(options,rows)
    background_refs=[reference_path(e.metadata['PFLTFILE']) for e in rows] if options['background']['amplifier_correction'] or options['background'].get('amplifier_assess') else []
    background=cache.key(raw+background_code+background_refs,
                         dict(background=options['background'],versions=versions()))
    alignment_parameters=dict(options['alignment'],relative=options['astrometry']['relative'],source_identification=options['astrometry']['source_identification'])
    idfiles=[Path(p) for p in (identification or {}).get('hst_catalogs',[])]
    alignment_code+=[package/'astrometry/identification.py',package/'astrometry/stellar.py',package/'astrometry/source_recovery.py',package/'astrometry/target_fallback.py',package/'astrometry/full_field.py',package/'astrometry/common_sources.py',package/'astrometry/wcs.py']
    if options.get('astrometry',{}).get('method')=='two_pass':
        alignment_parameters['mode']=options['astrometry']['relative']['mode'] if options['astrometry']['relative']['refine_if_needed'] else 'original'
    alignment=cache.key(raw+alignment_code+refs+idfiles,dict(background=background,alignment=alignment_parameters,versions=versions()))
    cr_options={k:v for k,v in options['drizzle'].items() if not k.startswith('final_')}
    cr=cache.key(copies+cr_code,dict(options=cr_options,versions=versions())) if all(p.exists() for p in copies) else None
    return background,alignment,cr

def reduction_key(cache,rows,options):
    dependencies=[]
    receipt=cache.receipt('uvis_absolute')
    if receipt.exists():dependencies=[Path(p) for p in __import__('json').loads(receipt.read_text()).get('dependencies',[]) if Path(p).exists()]
    return cache.key([e.path for e in rows]+code_inputs()+references(options,rows)+dependencies,dict(options=options,versions=versions()))

def registration_table(path,report):
    names=('exposure','nmatch','dx','dy','rotation','scale','xrms','yrms','accepted')
    rows=[[r.get(k) for k in names] for r in report.get('fits',[])]
    table=Table(rows=rows,names=names) if rows else Table(names=names,dtype=['U128','i8','f8','f8','f8','f8','f8','f8','bool'])
    table.meta={k:v for k,v in report.items() if k!='fits'}
    table.write(path,format='ascii.ecsv',overwrite=True)
