"""Target-scoped UVIS transactions with independent per-filter receipts."""
import copy,json,os,time
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
import yaml
from .paths import ReductionPaths
from .inventory import discover
from .grouping import groups
from .cache import StageCache
from .provenance import identity,versions
from .logging import reduction_logger
from .streams import stage_log

def merge(base,override):
    result=copy.deepcopy(base)
    for k,v in override.items():
        result[k]=merge(result[k],v) if isinstance(v,dict) and isinstance(result.get(k),dict) else copy.deepcopy(v)
    return result

def save(path,data,yaml_format=False):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(yaml.safe_dump(data,sort_keys=False) if yaml_format else json.dumps(data,indent=2,default=str,allow_nan=False))
    os.replace(temp,path)

@contextmanager
def target_lock(directory,group=None):
    import fcntl
    directory.mkdir(parents=True,exist_ok=True)
    name='.uvis.lock' if group is None else '.uvis-'+group+'.lock'
    with (directory.parent/name).open('a') as stream:
        try:fcntl.flock(stream,fcntl.LOCK_EX | (fcntl.LOCK_NB if group is not None else 0))
        except BlockingIOError:raise RuntimeError(f'Another reduction/migration is writing {directory}')
        try:yield
        finally:fcntl.flock(stream,fcntl.LOCK_UN)

def run_uvis(config,engine,dry_run=False,_retry_anchors=True,_retry_group=None,_inventory=None):
    from ..uvis.config import settings
    from ..uvis.storage import reduction_key
    from ..uvis.inventory import scan
    exposures,inventory_errors=_inventory if _inventory is not None else scan(config.input)
    allgroups=groups(exposures,'uvis')
    for error in inventory_errors:
        key=(error['target'],error['filter'],None)
        if not any(k[:2]==key[:2] for k in allgroups):allgroups[key]=[]
    selected=[(k,v) for k,v in allgroups.items() if (not config.target or k[0]==config.target) and (not config.element or k[1]==config.element) and (_retry_group is None or k==_retry_group)]
    if not selected:raise ValueError('No compatible selected exposures')
    astro=config.options.get('astrometry',{})
    if astro.get('method')=='two_pass' and astro.get('cross_filter',{}).get('enabled',False):
        import re
        requested=astro['cross_filter']['reference_filter']
        def priority(item):
            (t,f,s),rows=item;number=re.search(r'F(\d+)',f);wave=int(number[1]) if number else 0
            return (t,0 if requested==f else 1,0,-sum(e.metadata['EXPTIME'] for e in rows),f,s or '')
        selected.sort(key=priority)
    if config.setup and len({k[:2] for k,v in selected})!=len(selected):raise ValueError('Explicit setup would merge incompatible fields')
    roots={}
    for (t,f,s),es in allgroups.items():
        for e in es:
            key=(t,e.metadata['ROOTNAME'])
            if key in roots and roots[key]!=f:raise ValueError(f'Exposure ROOTNAME {key} is ambiguous across filters')
            roots[key]=f
    results=[]
    for (target,filt,setup),rows in selected:
        existing=Path(config.output)/target/'uvis'/'config.yaml'
        inherited=yaml.safe_load(existing.read_text()).get('filter_overrides',{}) if existing.exists() else {}
        override=merge(inherited.get(filt,{}),config.filter_overrides.get(filt,{}))
        options=settings(merge(config.options,override))
        # Resolve the shared reference directory before checking cached stages.
        options['calibration']['directory']=str(Path(options['calibration']['directory'] or os.environ.get('iref') or Path(config.output).parent/'calibration/uvis').expanduser().resolve())
        effective=replace(config,options=options,target=target,element=filt,setup=config.setup or setup)
        paths=ReductionPaths(Path(config.output),target,'uvis',filt,effective.setup)
        if dry_run:
            results.append(dict(inputs=[str(e.path) for e in rows],output=str(paths.directory),product=str(paths.product('drc')),group=[target,filt,setup],configuration=effective.effective()));continue
        # Legacy layouts need an explicit audited migration, never accidental overwrite.
        if any(paths.directory.glob('*/manifest.json')):raise RuntimeError(f'Legacy UVIS layout at {paths.directory}; run hstpipe uvis migrate first')
        from ..uvis.shallow import migrate_target
        migrate_target(paths.directory)
        with target_lock(paths.directory,paths.group_key):
            with target_lock(paths.directory):
                statepath=paths.directory/'manifest.json'
                state=json.loads(statepath.read_text()) if statepath.exists() else dict(schema_version=2,instrument='uvis',target=target,filters={})
                state.pop('workflow_error',None)
                if state.get('schema_version')!=2:raise ValueError('Unknown target manifest schema')
                cache=paths.cache();key=reduction_key(cache,rows,options)
                previous=state['filters'].get(paths.group_key,{})
                if paths.product('drc').exists() and not previous:raise FileExistsError(f'Untracked product would be overwritten: {paths.product("drc")}')
                if previous.get('status')=='complete' and config.resume and not config.force and not options.get('astrometry',{}).get('catalogs',{}).get('refresh',False) and cache.valid('reduction',key):
                    if paths.directory not in results:results.append(paths.directory)
                    from ..uvis.reporting import generate_report
                    generate_report(paths.directory)
                    continue
                if previous and not (config.resume or config.force):raise FileExistsError('Existing filter reduction; choose --resume or --force')
                if config.force and previous:
                    backup=paths.recovery_directory('history-'+paths.group_id+'-'+str(time.time_ns()));backup.mkdir(parents=True)
                    # Archive only this filter/setup. Other filters and shared state remain in place.
                    for artifact in previous.get('artifacts',[]):
                        p=Path(artifact['path'])
                        if p.is_file() and p.is_relative_to(paths.directory) and p.parent!=paths.catalog_directory():
                            dest=backup/p.relative_to(paths.directory);dest.parent.mkdir(parents=True,exist_ok=True);p.rename(dest)
                    save(backup/'manifest.json',previous)
                    for receipt in cache.directory.glob(cache.prefix+'*.json'):
                        receipt.rename(backup/receipt.name)
                cfgpath=paths.directory/'config.yaml'
                stored=yaml.safe_load(cfgpath.read_text()) if cfgpath.exists() else {}
                shared=config.effective();shared.update(target=target,element=None,setup=None,resume=True,force=False)
                shared['filter_overrides']=merge(stored.get('filter_overrides',{}),config.filter_overrides)
                save(cfgpath,shared,True)
                previous_products=previous.get('products') or previous.get('previous_products',[])
                # Existing products remain available even when this attempt fails.
                inputs=[identity(e.path) for e in rows]
                record=dict(filter=filt,setup=effective.setup,status='running',started_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),inputs=inputs,configuration=effective.effective(),
                            previous_products=previous_products,observing_configurations=[e.metadata for e in rows],grouping_decisions=[e.metadata.get('grouping_reason') for e in rows])
                state['filters'][paths.group_key]=record;save(statepath,state)
                log=reduction_logger(paths);log.info('Starting UVIS %s',paths.group_key)
            try:
                bad=[e for e in inventory_errors if e['target']==target and e['filter']==filt]
                if bad:
                    record['input_validation_errors']=bad
                    raise ValueError('Invalid FLC inputs block this filter: '+str(bad))
                with stage_log(paths.file('logs',None,'stages','log')):products=engine(effective,paths,rows)
                if [identity(e.path) for e in rows]!=inputs:raise RuntimeError('Original observations changed')
                qa=paths.file('diagnostics','drizzle','qa_summary')
                assessment=json.loads(qa.read_text()) if qa.exists() else {}
                files=set(paths._files)|set(products)|{Path(a['path']) for a in previous.get('artifacts',[])}
                artifacts=[identity(p) for p in sorted(files) if p.is_file() and p.is_relative_to(paths.directory) and 'logs' not in p.relative_to(paths.directory).parts]
                record.update(status=('complete' if products else assessment.get('status','incomplete')),versions=versions(),scientific_assessment=assessment,products=[identity(p) for p in products],artifacts=artifacts,
                              intermediate_files=[a for a in artifacts if 'intermediate' in Path(a['path']).relative_to(paths.directory).parts])
                with target_lock(paths.directory):
                    state=json.loads(statepath.read_text());state['filters'][paths.group_key]=record;save(statepath,state)
                if products:cache.record('reduction',reduction_key(cache,rows,options),[Path(a['path']) for a in artifacts])
                log.info('UVIS %s processing status: %s',paths.group_key,record['status'])
                if paths.directory not in results:results.append(paths.directory)
            except Exception as exc:
                record.update(status='failed',error=str(exc))
                with target_lock(paths.directory):
                    state=json.loads(statepath.read_text());state['filters'][paths.group_key]=record;save(statepath,state)
                log.exception('UVIS group failed')
                if paths.directory not in results:results.append(paths.directory)
            finally:
                from ..uvis.reporting import generate_report
                with target_lock(paths.directory):generate_report(paths.directory)
    if not dry_run and _retry_anchors and astro.get('method')=='two_pass' and astro.get('cross_filter',{}).get('enabled',False):
        # Optical suitability is known only after fitting. Give earlier failed
        # selected filters one chance to use a subsequently validated anchor.
        for (t,f,s),rows in selected:
            paths=ReductionPaths(Path(config.output),t,'uvis',f,config.setup or s)
            state=json.loads((paths.directory/'manifest.json').read_text())
            current=state['filters'].get(paths.group_key,{})
            assessment=current.get('scientific_assessment',{}).get('astrometry',{})
            if assessment.get('status')=='validated' or not assessment.get('attempts'):continue
            anchors=[v for k,v in state['filters'].items() if k!=paths.group_key and v.get('products') and v.get('scientific_assessment',{}).get('astrometry',{}).get('status')=='validated']
            if anchors:
                retry=replace(config,target=t,element=f,resume=True,force=False)
                run_uvis(retry,engine,False,_retry_anchors=False,_retry_group=(t,f,s),_inventory=(exposures,inventory_errors))
    return results
