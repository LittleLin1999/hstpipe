"""Audited, recoverable migration of legacy UVIS layouts; FITS bytes never change."""
import copy,json,os,re,shutil,time
from pathlib import Path
from collections import defaultdict
import yaml
from ..common.paths import ReductionPaths
from ..common.cache import StageCache
from ..common.provenance import identity,versions
from ..common.config import ReductionConfig
from ..common.target_runner import target_lock,save
from .config import settings
from .inventory import discover
from .grouping import group_exposures
from .storage import stage_keys,reduction_key,registration_table


def _delta(base,value):
    out={}
    for k,v in value.items():
        if isinstance(v,dict) and isinstance(base.get(k),dict):
            d=_delta(base[k],v)
            if d:out[k]=d
        elif base.get(k)!=v:out[k]=v
    return out


def _destination(paths,relative,roots):
    parts=relative.parts
    area=parts[0]
    if area not in ('products','intermediate','diagnostics','logs'):return None
    stage=parts[1] if area in ('intermediate','diagnostics') and len(parts)>2 else None
    if any(p.startswith(('run-','candidate-')) for p in parts):return None
    name=relative.name;stem=Path(name).stem;suffix=Path(name).suffix.lstrip('.')
    if area=='products':
        if name.endswith('_drc.fits'):return paths.product('drc')
        raise ValueError(f'Unexpected final product, requires review: {relative}')
    if area=='intermediate':
        for root in roots:
            for kind in ('sky','mask','amp','crmask','flc'):
                if name==f'{root}_{kind}.fits':return paths.exposure_file(area,stage,root,kind)
            match=re.fullmatch(re.escape(root)+r'_flc_(sci\d+_sources)\.cat',name)
            if match:return paths.exposure_file(area,stage,root,match[1],'cat')
    if area=='diagnostics':
        for root in roots:
            match=re.fullmatch(re.escape(root)+r'_flc_(.+)\.png',name)
            if match:return paths.exposure_file(area,stage,root,match[1],'png',group=True)
        if stage=='alignment' and stem=='fits':stem='alignment'
    prefix=f'{paths.target}_{paths.element}_'
    if stem.startswith(prefix):stem=stem[len(prefix):]
    return paths.file(area,stage,stem,suffix)


def _translate(value,mapping,legacy,paths):
    if isinstance(value,dict):return {k:_translate(v,mapping,legacy,paths) for k,v in value.items()}
    if isinstance(value,list):return [_translate(v,mapping,legacy,paths) for v in value]
    if not isinstance(value,str):return value
    if value in mapping:return str(mapping[value])
    # Candidate-directory references in alignment catalogs refer to exposures now
    # retained as prepared checkpoints. Preserve the recorded numerical solution.
    if value.startswith(str(legacy)):
        for old,new in mapping.items():
            if Path(value).name==Path(old).name:return str(new)
        return value.replace(str(legacy),str(paths.directory),1)
    return value


def migrate_uvis(output,*,target=None,dry_run=False):
    output=Path(output).expanduser().resolve();reports=[];inventories={}
    candidates=sorted(output.glob('*/uvis'))
    for directory in candidates:
        if target and directory.parent.name!=target:continue
        legacy=sorted(directory.glob('*/manifest.json'))+sorted(directory.glob('*/*/manifest.json'))
        legacy=[p for p in legacy if p.parent.relative_to(directory).parts[0] not in ('intermediate','diagnostics','logs','.cache')]
        if not legacy:continue
        plans=[];destinations={};owners={};common=None;overrides={}
        if (directory/'manifest.json').exists():raise ValueError(f'Mixed old/new manifests at {directory}; resolve before migrating')
        for manifest in legacy:
            old=json.loads(manifest.read_text());cfg=old.get('configuration') or yaml.safe_load((manifest.parent/'config.yaml').read_text())
            cfg=ReductionConfig(**cfg);opts=settings(cfg.options)
            if plans and cfg.input!=plans[0]['cfg'].input:raise ValueError('Legacy input directories differ; select a common raw directory before migration')
            if cfg.input not in inventories:inventories[cfg.input]=discover(cfg.input)
            allrows=inventories[cfg.input];filt=manifest.parent.relative_to(directory).parts[0]
            rows=[e for e in allrows if e.metadata['TARGNAME']==directory.parent.name and e.metadata['FILTER']==filt]
            inputnames=set(old.get('scientific_assessment',{}).get('inventory',{}))
            if not inputnames:inputnames={Path(i['path']).name for i in old.get('inputs',[])}
            if inputnames:rows=[e for e in rows if e.path.name in inputnames]
            if not rows:raise ValueError(f'Cannot identify FLC group for {manifest}')
            compatible=group_exposures([e for e in allrows if e.metadata['TARGNAME']==directory.parent.name and e.metadata['FILTER']==filt])
            matches=[(k,v) for k,v in compatible.items() if {e.path.name for e in v}=={e.path.name for e in rows}]
            setup=cfg.setup or (matches[0][0][2] if len(matches)==1 else None)
            if len(legacy)>1 and any(p['paths'].element==filt and p['paths'].setup==setup for p in plans):raise ValueError(f'Ambiguous setup identity for {manifest}')
            paths=ReductionPaths(output,directory.parent.name,'uvis',filt,setup)
            roots=[e.metadata['ROOTNAME'] for e in rows];mapping={}
            for src in sorted(manifest.parent.rglob('*')):
                if not src.is_file() or '.cache' in src.relative_to(manifest.parent).parts:continue
                if src.parent==manifest.parent:continue
                dst=_destination(paths,src.relative_to(manifest.parent),roots)
                if dst is None:continue
                ident=identity(src)
                if dst in owners and owners[dst]!=paths.group_key:raise ValueError(f'Ambiguous cross-filter/setup artifact: {dst}')
                owners[dst]=paths.group_key
                if dst in destinations and destinations[dst]['sha256']!=ident['sha256']:raise ValueError(f'Conflicting files map to {dst}')
                if dst.exists():raise FileExistsError(f'Will not overwrite {dst}')
                destinations[dst]=ident;mapping[str(src)]=dst
            cache=StageCache(manifest.parent/'.cache');valid={}
            for stage in ('reduction','uvis_upstream','uvis_cr'):
                receipt=cache.directory/(stage+'.json')
                valid[stage]=receipt.exists() and cache.valid(stage,json.loads(receipt.read_text())['key'])
            if common is None:common=opts
            diff=_delta(common,opts)
            if filt in overrides and overrides[filt]!=diff:raise ValueError(f'Setup-specific configuration conflict for {filt}; preserve legacy configuration and resolve explicitly')
            overrides[filt]=diff
            plans.append(dict(old=old,cfg=cfg,options=opts,rows=rows,paths=paths,mapping=mapping,legacy=manifest.parent,valid=valid))
        summary=dict(target=directory.parent.name,groups=[p['paths'].group_key for p in plans],files=len(destinations),
                     fits=sum(Path(k).suffix.lower()=='.fits' for p in plans for k in p['mapping']),status='planned')
        if dry_run:reports.append(summary);continue
        with target_lock(directory):
            backup=directory.with_name('.uvis-legacy-'+str(time.time_ns()))
            directory.rename(backup);directory.mkdir()
            try:
                state=dict(schema_version=2,instrument='uvis',target=directory.parent.name,filters={},migration=dict(backup=str(backup),fits_bytes_unchanged=True))
                for plan in plans:
                    paths=plan['paths'];mapping=plan['mapping'];legacydir=plan['legacy']
                    for oldpath,dest in mapping.items():
                        src=backup/Path(oldpath).relative_to(directory)
                        dest.parent.mkdir(parents=True,exist_ok=True)
                        if dest.exists():continue # identical planned shared artifact, verified below
                        if src.suffix.lower()=='.fits':os.link(src,dest)
                        elif src.suffix=='.json':save(dest,_translate(json.loads(src.read_text()),mapping,legacydir,paths))
                        elif src.suffix in ('.txt','.cat','.ecsv'):
                            text=src.read_text()
                            for old,new in sorted(mapping.items(),key=lambda a:len(a[0]),reverse=True):text=text.replace(old,str(new))
                            for old,new in mapping.items():
                                text=re.sub(re.escape(str(legacydir))+r'/intermediate/alignment/candidate-[^/\s]+/'+re.escape(Path(old).name),str(new),text)
                            text=re.sub(r'(?<![A-Za-z0-9_])reference_wcs\.fits',paths.file('intermediate','alignment','reference_wcs','fits').name,text)
                            dest.write_text(text)
                        else:shutil.copy2(src,dest)
                    for oldpath,dest in mapping.items():
                        if Path(oldpath).suffix.lower()=='.fits' and identity(dest)['sha256']!=destinations[dest]['sha256']:raise RuntimeError(f'FITS changed during migration: {dest}')
                    old=plan['old'];opts=plan['options'];rows=plan['rows']
                    oldup=backup/legacydir.relative_to(directory)/'intermediate/alignment/upstream.json'
                    upstream=_translate(json.loads(oldup.read_text()),mapping,legacydir,paths) if oldup.exists() else {}
                    bg=paths.file('intermediate','background','background');al=paths.file('intermediate','alignment','alignment')
                    reg=paths.file('intermediate','alignment','registration','ecsv')
                    if upstream:
                        save(bg,{k:upstream[k] for k in ('exposures','background')});save(al,upstream['alignment']);registration_table(reg,upstream['alignment'])
                    active=[e for e in rows if e.path.name not in opts['exclude'] and e.metadata['ROOTNAME'] not in opts['exclude']]
                    copies=[paths.exposure_file('intermediate','alignment',e.metadata['ROOTNAME'],'flc') for e in active]
                    cache=paths.cache();bgkey,alkey,crkey=stage_keys(cache,active,opts,copies)
                    if plan['valid']['uvis_upstream'] and upstream:
                        models=[dst for src,dst in mapping.items() if dst.parent==paths.area('intermediate','background') and dst.suffix=='.fits']
                        cache.record('uvis_background',bgkey,models+[bg]);cache.record('uvis_alignment',alkey,copies+[al,reg])
                    cr=paths.file('intermediate','drizzle','cr')
                    if plan['valid']['uvis_cr'] and crkey and cr.exists():
                        masks=[paths.exposure_file('intermediate','drizzle',e.metadata['ROOTNAME'],'crmask') for e in active]
                        cache.record('uvis_cr',crkey,masks+[cr])
                    qa=paths.file('diagnostics','drizzle','qa_summary')
                    assessment=json.loads(qa.read_text()) if qa.exists() else _translate(old.get('scientific_assessment',{}),mapping,legacydir,paths)
                    artifacts=[identity(dst) for dst in sorted(set(mapping.values())|{bg,al,reg}) if dst.exists() and 'logs' not in dst.relative_to(directory).parts]
                    products=[identity(paths.product('drc'))] if paths.product('drc').exists() else []
                    effective=plan['cfg'].effective();effective.update(target=paths.target,element=paths.element,setup=paths.setup,options=opts)
                    record=dict(filter=paths.element,setup=paths.setup,status=old['status'],inputs=[identity(e.path) for e in rows],configuration=effective,
                        observing_configurations=[e.metadata for e in rows],grouping_decisions=[e.metadata.get('grouping_reason') for e in rows],
                        scientific_assessment=assessment,products=products,artifacts=artifacts,
                        intermediate_files=[a for a in artifacts if 'intermediate' in Path(a['path']).relative_to(directory).parts],
                        versions=old.get('versions',{}),migration=dict(legacy=str(backup/legacydir.relative_to(directory)),receipts_preserved=plan['valid']))
                    state['filters'][paths.group_key]=record
                    if old['status']=='complete' and plan['valid']['reduction'] and products:
                        cache.record('reduction',reduction_key(cache,rows,opts),[Path(a['path']) for a in artifacts])
                config=plans[0]['cfg'].effective();config.update(target=directory.parent.name,element=None,setup=None,options=common,filter_overrides={k:v for k,v in overrides.items() if v},resume=True,force=False)
                save(directory/'config.yaml',config,True);save(directory/'manifest.json',state)
                summary.update(status='verified',backup=str(backup));save(directory/'migration.json',summary)
            except BaseException:
                failed=directory.with_name('.uvis-migration-incomplete-'+str(time.time_ns()));directory.rename(failed);backup.rename(directory)
                raise
        reports.append(summary)
    return reports
