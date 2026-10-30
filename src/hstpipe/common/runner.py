from __future__ import annotations
"""Shared run transaction; a failed run never becomes a valid cache entry."""
import json,shutil
from pathlib import Path
import yaml
from .config import ReductionConfig
from .inventory import discover
from .grouping import groups
from .paths import ReductionPaths
from .cache import StageCache
from .provenance import identity,write_manifest,versions
from .logging import reduction_logger

def run(config,engine,dry_run=False):
    if isinstance(config,(str,Path)):config=ReductionConfig.read(config)
    elif isinstance(config,dict):config=ReductionConfig(**config)
    exposures=discover(config.input,config.instrument);allgroups=groups(exposures,config.instrument)
    selected=[(k,v) for k,v in allgroups.items() if (not config.target or k[0]==config.target) and (not config.element or k[1]==config.element)]
    if not selected:raise ValueError('No compatible selected exposures')
    pairs=[k[:2] for k,v in selected]
    if len(pairs)!=len(set(pairs)):
        raise ValueError('Multiple incompatible setups share target/element; select a homogeneous input set before assigning setup')
    outputs=[]
    for (target,element,setup),rows in selected:
        paths=ReductionPaths(Path(config.output),target,config.instrument,element,config.setup)
        plan=dict(inputs=[str(e.path) for e in rows],output=str(paths.directory),configuration=config.effective())
        if dry_run:outputs.append(plan);continue
        cache=StageCache(paths.directory/'.cache')
        options=config.effective();options.pop('force');options.pop('resume')
        code=list((Path(__file__).parents[1]/config.instrument).glob('*.py'))
        refs=[]
        for key in ('reference_dir','control_dir'):
            if config.options.get(key):refs.extend(p for p in Path(config.options[key]).rglob('*') if p.is_file())
        key=cache.key([e.path for e in exposures]+code+refs,dict(config=options,versions=versions()))
        if config.resume and cache.valid('reduction',key):outputs.append(paths.directory);continue
        if paths.directory.exists():
            if not config.force:raise FileExistsError('Existing or stale reduction; choose --resume with matching cache or --force')
            # Preserve prior run for review, rather than overwrite or delete it.
            import time
            backup=paths.directory.with_name(paths.directory.name+'.previous-'+str(time.time_ns()))
            paths.directory.rename(backup)
        paths.directory.mkdir(parents=True)
        log=reduction_logger(paths);log.info('Starting %s: %s',config.instrument,plan['inputs'])
        (paths.directory/'config.yaml').write_text(yaml.safe_dump(config.effective(),sort_keys=False))
        inputs=[identity(e.path) for e in exposures]
        manifest=paths.directory/'manifest.json'
        write_manifest(manifest,status='running',inputs=inputs,configuration=config.effective())
        try:
            products=engine(config,paths,rows)
            if [identity(e.path) for e in exposures]!=inputs:raise RuntimeError('Original observations changed')
            write_manifest(manifest,status='complete',inputs=inputs,configuration=config.effective(),
                products=[identity(p) for p in products],calibration=[identity(p) for p in refs])
            artifacts=[p for p in paths.directory.rglob('*') if p.is_file() and '.cache' not in p.parts and 'logs' not in p.parts]
            cache.record('reduction',key,artifacts);log.info('Complete: %s',products)
            outputs.append(paths.directory)
        except Exception as exc:
            write_manifest(manifest,status='failed',inputs=inputs,configuration=config.effective(),error=str(exc));log.exception('Reduction failed');raise
    return outputs
