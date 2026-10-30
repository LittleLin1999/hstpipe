"""FLC-only imaging workflow; no detector calibration or CTE rerun."""
import json,shutil,os
from pathlib import Path
from ..common.runner import run
from .qa import inspect_flc,validate_drizzle
from .background import assess_regions
from .alignment import align
from .drizzle import drizzle

def _engine(config,paths,exposures):
    opts=config.options
    if set(opts)-{'alignment','drizzle','align','background_masks','amplifier_masks'}:raise ValueError('Unsupported UVIS option')
    directory=paths.area('intermediate','alignment',create=True)
    copies=[];qa={}
    for e in exposures:
        p=directory/e.path.name;shutil.copy2(e.path,p);copies.append(p);qa[p.name]=inspect_flc(p)
        for option in ('background_masks','amplifier_masks'):
            if e.path.name in opts.get(option,{}):
                assessment=assess_regions(p,opts[option][e.path.name])
                folder='background' if option=='background_masks' else 'alignment'
                (paths.area('diagnostics',folder,create=True)/(p.stem+'_'+option+'.json')).write_text(json.dumps(assessment,indent=2))
    (paths.area('diagnostics','alignment',create=True)/'input_quality.json').write_text(json.dumps(qa,indent=2))
    # DrizzlePac writes auxiliaries in cwd; confine its global cwd convention.
    previous=Path.cwd();work=paths.area('intermediate','drizzle',create=True)
    try:
        os.chdir(work)
        if opts.get('align',True):align(copies,opts.get('alignment'))
        prefix=work/f'{paths.target}_{paths.element}'
        drizzle(copies,prefix,opts.get('drizzle'))
        product=prefix.with_name(prefix.name+'_drc.fits')
        if not product.exists():raise FileNotFoundError(f'AstroDrizzle did not create {product}')
        validate_drizzle(product)
        final=paths.product('drc');final.parent.mkdir(parents=True,exist_ok=True);shutil.move(product,final)
        (paths.area('diagnostics','drizzle',create=True)/'validation.json').write_text(json.dumps({'SCI_WHT_CTX_valid':True,'alignment':opts.get('align',True),'ninput':len(copies)}))
        return [final]
    finally:os.chdir(previous)

def reduce_uvis(config,*,dry_run=False):return run(config,_engine,dry_run)
