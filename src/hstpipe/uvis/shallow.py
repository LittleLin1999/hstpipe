"""Recoverable path-only migration of target-level UVIS artifacts and receipts."""
import json,os,re,shutil,time
from pathlib import Path
from astropy.table import Table
from ..common.paths import ReductionPaths
from ..common.cache import StageCache
from ..common.provenance import identity,versions
from ..common.target_runner import target_lock,save

def _replace(value,mapping):
    if isinstance(value,dict):
        out={_replace(k,mapping):_replace(v,mapping) for k,v in value.items()}
        if {'path','size','sha256'}<=set(out) and Path(out['path']).is_file():out.update(identity(out['path']))
        return out
    if isinstance(value,list):return [_replace(v,mapping) for v in value]
    if isinstance(value,str):
        for old,new in mapping:value=value.replace(old,new)
    return value

def _destination(paths,rel):
    p=rel.parts
    if '.history' in p:return None
    if p[0]=='.cache' and len(p)==3:return paths.area('intermediate','alignment')/(p[1]+'_cache_'+p[2])
    if p[0] in ('diagnostics','logs'):
        stage=p[1] if len(p)>2 else None
        return paths.area(p[0])/paths.qualified_filename(p[0],stage,rel.name)
    if p[0]=='intermediate':
        if 'catalogs' in p:
            name='catalog_cache_'+rel.name if '.cache' in p else rel.name
            return paths.catalog_directory()/name
        stage='working' if rel.name.endswith('_flc.fits') else p[1]
        if stage not in ('catalogs','background','alignment','drizzle','working'):stage='alignment'
        return paths.area('intermediate',stage)/rel.name
    if len(p)==1 and p[0] not in ('config.yaml','manifest.json'):
        return paths.area('intermediate','alignment')/rel.name
    return paths.directory/rel

def _rekey(paths,record,verified):
    """Rebase only verified path-only checkpoints; no science acceptance upgrade."""
    from .config import settings
    from .storage import stage_keys,reduction_key
    from .inventory import discover
    from .astrometry.identification import identification_key
    from .astrometry.full_field import detection_key
    cfg=record.get('configuration',{});opts=settings(cfg.get('options',{}));names={Path(p['path']).name for p in record.get('inputs',[])}
    rows=[e for e in discover(cfg['input']) if e.path.name in names]
    active=[e for e in rows if e.path.name not in opts['exclude'] and e.metadata['ROOTNAME'] not in opts['exclude']]
    raw=[e.path for e in active];cache=paths.cache()
    copies=[paths.exposure_file('intermediate','working',e.metadata['ROOTNAME'],'flc') for e in active]
    idfile=paths.file('intermediate','alignment','source_identification')
    identification=json.loads(idfile.read_text()) if idfile.exists() else {}
    bg,al,cr=stage_keys(cache,active,opts,copies,identification)
    # A layout migration cannot certify a different scientific DQ/sky policy.
    from .astrometry.centroids import QUALITY_VERSION
    sourcefile=paths.file('intermediate','alignment','full_field_sources')
    sources=json.loads(sourcefile.read_text()) if sourcefile.exists() else {}
    source_current=bool(sources.get('exposures')) and all(m.get('version')==QUALITY_VERSION for e in sources['exposures'] for m in e.get('methods',[]))
    bg_current=False
    for e in active:
        sky=paths.exposure_file('intermediate','background',e.metadata['ROOTNAME'],'sky')
        if sky.exists():
            from astropy.io import fits
            with fits.open(sky) as hd:bg_current=all('HPSKYBG' in h.header for h in hd if h.name=='SCI')
        if not bg_current:break
    keys={}
    if bg_current:keys['uvis_background']=bg
    if source_current:
        keys.update(uvis_sources=detection_key(cache,raw,opts['astrometry']['relative']),uvis_identification=identification_key(cache,raw,opts['astrometry']))
    # Failed registration is deliberately retried: the target-local fallback may help.
    assessment=record.get('scientific_assessment',{});relative=assessment.get('relative_astrometry',assessment.get('alignment',{}))
    if source_current and bg_current and relative.get('may_combine') and relative.get('status') in ('validated','single_exposure'):
        keys['uvis_alignment']=al
        if cr:keys['uvis_cr']=cr
        if copies and all(p.exists() for p in copies):
            options={k:v for k,v in opts['drizzle'].items() if not k.startswith('final_')}
            keys['uvis_provisional']=cache.key(copies+[Path(__file__).with_name('drizzle.py')],dict(cr=options,provisional=opts['astrometry']['provisional'],versions=versions()))
        if assessment.get('astrometry',{}).get('status')=='validated':keys['reduction']=reduction_key(cache,rows,opts)
    rebased=[]
    for stage,key in keys.items():
        receipt=cache.receipt(stage)
        if stage in verified and receipt.exists():
            data=json.loads(receipt.read_text())
            if not data.get('outputs') or not all(Path(o['path']).is_file() and identity(o['path'])==o for o in data['outputs']):continue
            data['layout_previous_key']=data['key'];data['key']=key;data['layout_rebased']=True
            save(receipt,data);rebased.append(stage)
    return rebased

def migrate_target(directory):
    directory=Path(directory)
    if not directory.exists() or not (directory/'manifest.json').exists():return {'status':'not_needed'}
    nested=any(p.is_dir() for area in ('diagnostics','logs') for p in (directory/area).glob('*'))
    nested|=any(p.is_dir() for p in (directory/'intermediate').glob('*/*'))
    nested|=(directory/'.cache').exists() or (directory/'.history').exists()
    nested|=any((directory/'intermediate/alignment').glob('*_flc.fits'))
    if not nested:return {'status':'already_shallow'}
    with target_lock(directory):
        state=json.loads((directory/'manifest.json').read_text());paths=ReductionPaths(directory.parent.parent,state['target'],'uvis','layout')
        files=[p for p in directory.rglob('*') if p.is_file()];plan={};owners={};verified={}
        for key,r in state.get('filters',{}).items():
            group=ReductionPaths(paths.root,paths.target,'uvis',r['filter'],r.get('setup'));old=StageCache(directory/'.cache'/group.group_id)
            verified[key]=[]
            for receipt in old.directory.glob('*.json'):
                try:
                    data=json.loads(receipt.read_text())
                    if old.valid(receipt.stem,data['key']):verified[key].append(receipt.stem)
                except (OSError,ValueError,KeyError):pass
        for source in files:
            dest=_destination(paths,source.relative_to(directory))
            if dest is None:continue
            if dest in owners and source.read_bytes()!=owners[dest].read_bytes():
                raise ValueError(f'Ambiguous shallow migration mapping: {source} and {owners[dest]} -> {dest}')
            plan[source]=dest;owners[dest]=source
        backup=paths.recovery_directory('layout-'+str(time.time_ns()))
        mapping=sorted([(str(s),str(d)) for s,d in plan.items()],key=lambda kv:len(kv[0]),reverse=True)
        # Historical products stay recoverable outside the active shallow tree.
        mapping.append((str(directory/'.history'),str(backup/'.history')))
        directory.rename(backup);directory.mkdir()
        try:
            for source,dest in plan.items():
                original=backup/source.relative_to(directory);dest.parent.mkdir(parents=True,exist_ok=True)
                if dest.exists():continue
                if original.suffix.lower() in ('.json','.ecsv','.cat','.txt','.yaml','.log'):
                    text=original.read_text(errors='strict')
                    for old,new in mapping:text=text.replace(old,new)
                    dest.write_text(text)
                else:os.link(original,dest)
            # Refresh identities after every textual reference has been translated.
            for dest in set(plan.values()):
                if dest.suffix=='.json':
                    try:save(dest,_replace(json.loads(dest.read_text()),mapping))
                    except (ValueError,UnicodeError):pass
            state=json.loads((directory/'manifest.json').read_text());rebased={}
            for key,r in state.get('filters',{}).items():
                group=ReductionPaths(paths.root,paths.target,'uvis',r['filter'],r.get('setup'))
                rebased[key]=_rekey(group,r,verified.get(key,[]))
            summary={'status':'verified','backup':str(backup),'files':len(plan),'rebased_stages':rebased,'fits_bytes_unchanged':True}
            save(paths.file('intermediate','alignment','migration'),summary)
            state['layout_migration']=summary;save(directory/'manifest.json',state)
            return summary
        except BaseException:
            failed=paths.recovery_directory('migration-incomplete-'+str(time.time_ns()));directory.rename(failed);backup.rename(directory)
            raise
