"""Target-local barrier between internal preparation and external astrometry."""
from dataclasses import replace
import json
from pathlib import Path
from ...common.runner import run
from ...common.target_runner import run_uvis
from .cross_filter import internal_cross_filter


def relative_then_absolute(config,engine):
    # Inventory/grouping and per-filter overrides stay owned by the common runner.
    from ..inventory import scan
    inventory=scan(config.input)
    plans=run_uvis(config,engine,dry_run=True,_inventory=inventory);targets=list(dict.fromkeys(p['configuration']['target'] for p in plans));outputs=[]
    for target in targets:
        directory=Path(config.output)/target/'uvis'
        try:
            selected=[p for p in plans if p['configuration']['target']==target]
            cfg=replace(config,target=target)
            def prepare(c,paths,rows):return engine(c,paths,rows,prepare_only=True)
            roots=run_uvis(cfg,prepare,_retry_anchors=False,_inventory=inventory)
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
            for plan in selected:finish(plan)
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
