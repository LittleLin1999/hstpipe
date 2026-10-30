"""Infrastructure adapter around the preserved six-stage science workflow."""
from pathlib import Path
from ..common.runner import run

def _engine(config,paths,exposures):
    from .config import Config
    from .workflow import run as legacy_run
    # A private flat input view selects exactly this group without renaming data.
    view=paths.area('intermediate','flt/inputs',create=True)
    for e in exposures:
        for suffix in ('flt','x2d'):
            source=e.path.with_name(e.path.name.replace('_flt.fits','_'+suffix+'.fits'))
            if not source.exists():raise FileNotFoundError(f'Required STIS counterpart: {source}')
            (view/source.name).symlink_to(source)
    options=dict(config.options)
    options.update(project_dir=str(paths.directory),input_glob='intermediate/flt/inputs/*_flt.fits',
                   mast_dir=str(view),output_dir=str(paths.directory),name=f'{paths.target}_{paths.element}')
    legacy_run(Config(**options))
    # Keep native manifest intact alongside the shared provenance manifest.
    (paths.directory/'manifest.json').rename(paths.area('intermediate')/'stis_manifest.json')
    return [paths.product(kind) for kind in ('s2d_drz','x2d_coadd','x1d') if paths.product(kind).exists()]

def reduce_stis(config,*,dry_run=False):return run(config,_engine,dry_run)

def stages():
    """Explicit legacy scientific stages accepting Config and their native inputs."""
    from . import workflow
    return {name:getattr(workflow,'stage'+str(i)+'_'+name) for i,name in enumerate(
        ('prepare','background','joint_drizzle','mast_coadd','validate','extract'),1)}
