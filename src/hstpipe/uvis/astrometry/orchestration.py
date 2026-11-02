"""Target-local barrier between internal preparation and external astrometry."""
from dataclasses import replace
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from functools import partial
import json
from pathlib import Path
from ...common.runner import run
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
        directory=Path(config.output)/target/'uvis'
        try:
            selected=[p for p in plans if p['configuration']['target']==target]
            cfg=replace(config,target=target)
            if config.force:
                roots=run_uvis(cfg,partial(_prepare_engine,engine),_retry_anchors=False,_inventory=inventory)
            else:
                roots=[root for group in _filters(cfg,engine,inventory,selected,True) for root in group]
            for root in roots:
                if root not in outputs:outputs.append(root)
            directory=Path(config.output)/target/'uvis'
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
            # Complete the selected anchor before concurrent independent filters,
            # preserving cross-filter absolute-frame transfer dependencies.
            if selected:finish(selected[0])
            eligible=[]
            state=json.loads((directory/'manifest.json').read_text())['filters']
            for plan in selected[1:]:
                c=plan['configuration'];key=c['element']+('_'+c['setup'] if c['setup'] else '')
                if state[key].get('scientific_assessment',{}).get('alignment',{}).get('status') in ('validated','single_exposure'):eligible.append(plan)
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
