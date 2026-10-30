"""Restartable six-stage orchestration around unchanged STIS science functions."""
from pathlib import Path
import json, shutil, time
from ..common.runner import run
from ..common.cache import StageCache
from ..common.provenance import identity, versions

def stages():
    from . import workflow
    return {name:getattr(workflow,"stage"+str(i)+"_"+name) for i,name in enumerate(
        ("prepare","background","joint_drizzle","mast_coadd","validate","extract"),1)}

def _engine(config,paths,exposures):
    from .config import Config
    from .io import read_inputs
    from .geometry import Grid
    from astropy.table import Table
    root=paths.directory
    paths.area("products",create=True)
    for area in ("background","coadd","extraction","fluxcal"):
        paths.area("diagnostics",area,create=True)
    view=paths.area("intermediate","flt/inputs",create=True)
    originals=[]
    for e in exposures:
        for suffix in ("flt","x2d"):
            source=e.path.with_name(e.path.name.replace("_flt.fits","_"+suffix+".fits"))
            if not source.exists():raise FileNotFoundError(source)
            q=view/source.name
            if q.exists() and q.resolve()!=source.resolve():raise ValueError("Stale input view")
            if not q.exists():q.symlink_to(source)
            originals.append(source)
    options=dict(config.options)
    # HSTCAL is resolved from the active standalone environment PATH.
    options.update(project_dir=str(root),input_glob="intermediate/flt/inputs/*_flt.fits",
        mast_dir=str(view),output_dir=str(root),name=f"{paths.target}_{paths.element}")
    options.setdefault("control_dir",str(root/"intermediate/flt/registration/official"))
    native=Config(**options)
    (root/"config.json").write_text(json.dumps(native.serializable(),indent=2))
    inputs=read_inputs(sorted(view.glob("*_flt.fits")))
    functions=stages();cache=StageCache(root/".cache/stages")
    code=list(Path(__file__).parent.glob("*.py"))
    refs=list(native.path(native.reference_dir).glob("*"))
    refs=[p for p in refs if p.is_file()]
    executables=[Path(shutil.which(name)) for name in ("cs6.e","cs7.e") if shutil.which(name)]
    previous=cache.key(originals+refs+code+executables,dict(config=native.serializable(),versions=versions()))
    def stage(name,targets,callback):
        nonlocal previous
        key=cache.key([],dict(upstream=previous,stage=name))
        if cache.valid(name,key):
            print("CACHE",name,flush=True);previous=key;return
        for relative in targets:
            p=root/relative
            if p.exists():
                backup=root/".previous"/(str(time.time_ns())+"-"+name)/relative
                backup.parent.mkdir(parents=True,exist_ok=True);shutil.move(p,backup)
        print("STAGE",name,flush=True);callback()
        outputs=[]
        for relative in targets:
            p=root/relative
            outputs.extend([p] if p.is_file() else [a for a in p.rglob("*") if a.is_file()])
        cache.record(name,key,outputs);previous=key
    stage("prepare",["intermediate/flt/registration"],lambda:functions["prepare"](native,inputs))
    meta=json.loads((root/"intermediate/flt/registration/inputs.json").read_text())
    grid=Grid(**meta["grid"])
    stage("background",["intermediate/background"],lambda:functions["background"](native,inputs))
    reports=json.loads((root/"intermediate/background/models.json").read_text())
    clean=paths.product("s2d_drz");mast=paths.product("x2d_coadd")
    stage("joint_drizzle",["intermediate/flt/coadd",str(clean.relative_to(root))],lambda:functions["joint_drizzle"](native,inputs,grid))
    stage("mast_coadd",["intermediate/x2d" ,str(mast.relative_to(root)),"diagnostics/coadd"],lambda:functions["mast_coadd"](native,inputs,grid))
    hashes={str(p):identity(p)["sha256"] for p in originals}
    saved=root/"intermediate/stis_manifest.json"
    def validate():
        paths.area("diagnostics","fluxcal",create=True)
        functions["validate"](native,[clean,mast],grid,hashes,reports)
        shutil.copy2(root/"manifest.json",saved)
        for name in ('REPORT.md','report.json','comparison.png'):
            q=root/'diagnostics/coadd'/name
            if q.exists():shutil.move(q,root/'diagnostics/fluxcal'/name)
    stage("validate",["intermediate/stis_manifest.json","diagnostics/fluxcal"],validate)
    shutil.copy2(saved,root/"manifest.json")
    if native.extraction.enabled:
        stage("extract",["intermediate/x1d","diagnostics/extraction",str(paths.product("x1d").relative_to(root))],lambda:functions["extract"](native))
    rows=[dict(input=e.path.name,postarg2_arcsec=float(e.p["POSTARG2"]),
        configured_offset_arcsec=float(native.source_offsets.get(e.path.name,0)),
        spatial_scale_arcsec=float(meta["native_spatial_scale"])) for e in inputs]
    Table(rows=rows).write(root/"intermediate/registration.ecsv",format="ascii.ecsv",overwrite=True)
    # Display artifacts belong in diagnostics; numerical exposure models remain intermediate.
    for directory in (root/"intermediate/background",root/"intermediate/x2d"):
        for p in directory.glob("*.png"):
            shutil.copy2(p,paths.area("diagnostics","background")/p.name)
    return [p for p in (clean,mast,paths.product("x1d")) if p.exists()]

def reduce_stis(config,*,dry_run=False):return run(config,_engine,dry_run)
