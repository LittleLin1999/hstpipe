"""Target-local barrier between internal preparation and external astrometry."""
from dataclasses import replace
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from functools import partial
import json
from pathlib import Path
from ...common.runner import run
from ...common.paths import ReductionPaths
from ...common.target_runner import run_uvis
from .cross_filter import internal_cross_filter


def _prepare_engine(engine,c,paths,rows):
    return engine(c,paths,rows,prepare_only=True)

def _run_filter(config,engine,inventory,plan,prepare=False):
    c=plan['configuration']
    worker=partial(_prepare_engine,engine) if prepare else engine
    return run_uvis(replace(config,target=c['target'],element=c['element'],resume=True,force=False),worker,
        _retry_anchors=False,_retry_group=tuple(plan['group']),_inventory=inventory)

def _filters(config,engine,inventory,plans,prepare=False):
    workers=int(config.options.get('execution',{}).get('filter_workers',1))
    if workers==1:return [_run_filter(config,engine,inventory,p,prepare) for p in plans]
    # Native STScI routines change cwd and global streams: isolate each filter.
    with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
        futures=[pool.submit(_run_filter,config,engine,inventory,p,prepare) for p in plans]
        return [f.result() for f in futures]


def relative_then_absolute(config,engine):
    # Inventory/grouping and per-filter overrides stay owned by the common runner.
    from ..inventory import scan
    inventory=scan(config.input)
    plans=run_uvis(config,engine,dry_run=True,_inventory=inventory);targets=list(dict.fromkeys(p['configuration']['target'] for p in plans));outputs=[]
    target_workers=int(config.options.get('execution',{}).get('target_workers',1))
    if target_workers>1 and len(targets)>1:
        options=dict(config.options,execution=dict(config.options.get('execution',{}),target_workers=1))
        with ProcessPoolExecutor(max_workers=target_workers,mp_context=get_context('spawn')) as pool:
            jobs=[pool.submit(relative_then_absolute,replace(config,target=t,options=options),engine) for t in targets]
            return [root for job in jobs for root in job.result()]
    for target in targets:
        directory=ReductionPaths(Path(config.output),target,'uvis','layout').directory
        try:
            selected=[p for p in plans if p['configuration']['target']==target]
            cfg=replace(config,target=target)
            if config.force:
                roots=run_uvis(cfg,partial(_prepare_engine,engine),_retry_anchors=False,_inventory=inventory)
            else:
                roots=[root for group in _filters(cfg,engine,inventory,selected,True) for root in group]
            for root in roots:
                if root not in outputs:outputs.append(root)
            directory=ReductionPaths(Path(config.output),target,'uvis','layout').directory
            keys=[p['configuration']['element']+('_'+p['configuration']['setup'] if p['configuration']['setup'] else '') for p in selected]
            cross=internal_cross_filter(directory,keys,config.options['astrometry'])
            # Solve the measured internal anchor first, without a wavelength preference.
            selected.sort(key=lambda p:0 if cross.get(p['configuration']['element']+('_'+p['configuration']['setup'] if p['configuration']['setup'] else ''),{}).get('status')=='anchor' else 1)
            def finish(p):
                c=p['configuration'];key=c['element']+('_'+c['setup'] if c['setup'] else '')
                state=json.loads((directory/'manifest.json').read_text())['filters'][key]
                a=state.get('scientific_assessment',{})
                if a.get('alignment',{}).get('status') not in ('validated','single_exposure'):return
                # Force archives once in preparation; finalization resumes those checkpoints.
                run_uvis(replace(cfg,element=c['element'],resume=True,force=False),engine,_retry_anchors=False,
                         _retry_group=tuple(p.get('group',(target,c['element'],c['setup']))),_inventory=inventory)
            # Parent WCS must be absolute before a descendant transfers it.
            # Independent components and equal-depth filters may run together.
            def depth(plan):
                c=plan['configuration'];key=c['element']+('_'+c['setup'] if c['setup'] else '')
                return cross.get(key,{}).get('depth',0 if cross.get(key,{}).get('status')=='anchor' else 1)
            for level in sorted({depth(p) for p in selected}):
                state=json.loads((directory/'manifest.json').read_text())['filters']
                eligible=[]
                for plan in selected:
                    c=plan['configuration'];key=c['element']+('_'+c['setup'] if c['setup'] else '')
                    if depth(plan)==level and state[key].get('scientific_assessment',{}).get('alignment',{}).get('status') in ('validated','single_exposure'):
                        eligible.append(plan)
                _filters(cfg,engine,inventory,eligible)
            # If the initial anchor could not be tied externally but another member could,
            # give selected, internally connected filters one transfer attempt.
            state=json.loads((directory/'manifest.json').read_text())['filters']
            if any(r.get('scientific_assessment',{}).get('astrometry',{}).get('status')=='validated' for r in state.values()):
                for plan in selected:
                    c=plan['configuration'];key=c['element']+('_'+c['setup'] if c['setup'] else '')
                    a=state[key].get('scientific_assessment',{})
                    if cross.get(key,{}).get('status') in ('validated','anchor') and a.get('astrometry',{}).get('status') not in ('validated','final_validation_failed'):finish(plan)
        except Exception as exc:
            from ...common.target_runner import save
            directory.mkdir(parents=True,exist_ok=True)
            statefile=directory/'manifest.json'
            state=json.loads(statefile.read_text()) if statefile.exists() else dict(schema_version=2,instrument='uvis',target=target,filters={})
            state['workflow_error']=f'{type(exc).__name__}: {exc}';save(statefile,state)
            print(f'UVIS {target}: {state["workflow_error"]}',flush=True)
        finally:
            from ..reporting import generate_report
            if (directory/'manifest.json').exists():generate_report(directory)
            if directory not in outputs:outputs.append(directory)
    return outputs


def resume_cross_filters(directory, selected=None):
    """Refresh cross-filter/final stages from validated provisional checkpoints.

    No inventory, background or within-filter registration is recomputed. Old
    products remain published when a new candidate fails any scientific gate.
    """
    import traceback
    from contextlib import redirect_stdout, redirect_stderr
    from ...common.paths import ReductionPaths
    from ...common.inventory import Exposure
    from ...common.provenance import identity
    from ...common.target_runner import save, target_lock
    from ..config import settings
    from ..reporting import generate_report
    from .workflow import finalize
    directory=Path(directory);statefile=directory/'manifest.json';results={}
    try:
        cross=internal_cross_filter(directory,selected)
        for key in sorted(cross,key=lambda k:(cross[k].get('depth',0),k)):
            edge=cross[key]
            if edge['status']!='validated':
                results[key]=dict(status='unchanged',cross_filter=edge['status'],reason=edge.get('reason','Reference frame anchor'))
                continue
            with target_lock(directory):state=json.loads(statefile.read_text())
            record=state['filters'][key];assessment=record.get('scientific_assessment',{})
            prior=results.get(edge['anchor_group'],{})
            if prior.get('status') not in (None,'unchanged','science_ready','complete_with_limitations','validated_existing_product'):
                results[key]=dict(status='blocked',reason='Parent final validation did not pass in this refresh');continue
            parent=state['filters'].get(edge['anchor_group'],{}).get('scientific_assessment',{}).get('astrometry',{})
            if parent.get('status')!='validated':
                results[key]=dict(status='blocked',reason='Parent absolute WCS is not validated');continue
            paths=ReductionPaths.from_directory(directory,record['filter'],record.get('setup'),target=state['target'])
            paths.area('logs',create=True);log=paths.file('logs',None,'cross_filter_resume','log')
            try:
                with log.open('w') as stream,redirect_stdout(stream),redirect_stderr(stream):
                    opts=settings(record['configuration'].get('options',{}))
                    # This focused resume transfers an existing absolute frame;
                    # it must not fall back to a new external fit if transfer fails.
                    checkpoint=paths.file('intermediate','drizzle','provisional_inputs')
                    expected=json.loads(checkpoint.read_text());working=[Path(r['path']) for r in expected]
                    if [identity(p) for p in working]!=expected:raise ValueError('Relative FLC checkpoint differs from provisional drizzle inputs')
                    rows=[Exposure(Path(i['path']),metadata,'flc',i) for i,metadata in zip(record['inputs'],record['observing_configurations'])]
                    if [identity(e.path) for e in rows]!=record['inputs']:raise ValueError('Original FLC identity changed since inventory')
                    provisional=paths.file('intermediate','drizzle','provisional','fits')
                    from .alignment import _anchor,observation_epoch
                    from .cross_filter import transfer_absolute_frame
                    anchor=_anchor(paths,observation_epoch(rows)[0],opts['astrometry'])
                    if anchor is None:raise ValueError('No validated compatible parent absolute frame')
                    from .cross_filter import validate_existing_frame
                    prior_product=record.get('products',[])
                    existing=validate_existing_frame(Path(prior_product[0]['path']),paths,opts['astrometry'],anchor) if prior_product else {}
                    if existing.get('status')=='validated':
                        edge['final_validation']=existing
                        assessment['cross_filter_alignment']=edge
                        results[key]=dict(status='validated_existing_product',products=[],validation=existing,log=str(log))
                        record['scientific_assessment']=assessment
                        artifacts={item['path']:item for item in record.get('artifacts',[]) if Path(item['path']).exists()}
                        for name,suffix in [('cross_filter_existing_final','json'),('cross_filter_existing_matches','ecsv')]:
                            artifact=paths.file('intermediate','alignment',name,suffix)
                            if artifact.exists():artifacts[str(artifact)]=identity(artifact)
                        record['artifacts']=list(artifacts.values())
                        record['intermediate_files']=[item for item in artifacts.values() if '/intermediate/' in item['path']]
                        with target_lock(directory):
                            state=json.loads(statefile.read_text());state['filters'][key]=record;save(statefile,state)
                            save(paths.file('diagnostics','drizzle','qa_summary'),assessment)
                        print(paths.target,key,'validated_existing_product',flush=True)
                        continue
                    transfer_absolute_frame(provisional,rows,paths,opts['astrometry'],anchor)
                    products,result=finalize(provisional,working,rows,paths,opts['astrometry'],opts['drizzle'],
                        assessment.get('relative_astrometry',assessment['alignment']),opts['target_cr_fraction'],reuse_absolute=True)
                    if [identity(e.path) for e in rows]!=record['inputs']:raise ValueError('Original FLC identity changed during finalization')
                results[key]=dict(status=result['status'],products=[str(p) for p in products],astrometry=result.get('astrometry',{}),log=str(log))
                if products:
                    assessment.update(result)
                    assessment['provisional_image']=dict(path=str(provisional),relative_status=assessment['relative_astrometry']['status'])
                    for stale in ('unvalidated_final','source_review','independent_cr_recheck','cross_filter_reprocessing'):assessment.pop(stale,None)
                    assessment.setdefault('stages',{}).update(absolute_alignment=dict(status='validated',reference_mode='cross_filter'),tweakback=dict(status=result['tweakback']['status']),final_drizzle=dict(status='complete',outputs=[str(p) for p in products]),final_validation=dict(status=result['astrometry']['final_validation']['status']))
                    record.update(products=[identity(p) for p in products],status='complete')
                else:
                    # The previous production assessment remains true for the
                    # retained product. The new failed attempt is separate.
                    assessment['cross_filter_reprocessing']=result
                assessment['cross_filter_alignment']=edge
            except Exception as exc:
                with log.open('a') as stream:stream.write(traceback.format_exc())
                results[key]=dict(status='failed',reason=f'{type(exc).__name__}: {exc}',log=str(log))
                assessment['cross_filter_reprocessing']=results[key]
            record['scientific_assessment']=assessment
            with target_lock(directory):
                state=json.loads(statefile.read_text());state['filters'][key]=record;save(statefile,state)
                save(paths.file('diagnostics','drizzle','qa_summary'),assessment)
                artifacts={item['path']:item for item in record.get('artifacts',[]) if Path(item['path']).exists()}
                for artifact in paths._files:
                    if artifact.is_file() and artifact.is_relative_to(directory):artifacts[str(artifact)]=identity(artifact)
                record['artifacts']=list(artifacts.values())
                record['intermediate_files']=[item for item in artifacts.values() if '/intermediate/' in item['path']]
                state['filters'][key]=record;save(statefile,state)
            print(paths.target,key,results[key]['status'],flush=True)
        with target_lock(directory):
            state=json.loads(statefile.read_text());state['cross_filter_refresh']=results;save(statefile,state)
    finally:
        if statefile.exists():generate_report(directory)
    return results
