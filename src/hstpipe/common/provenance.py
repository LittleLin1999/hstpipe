from __future__ import annotations
import hashlib,json,platform
from pathlib import Path
from importlib.metadata import version,PackageNotFoundError

def identity(path):
    path=Path(path);h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return dict(path=str(path.resolve()),size=path.stat().st_size,sha256=h.hexdigest())

def versions():
    result={'python':platform.python_version()}
    for name in ('hstpipe','numpy','astropy','scipy','matplotlib','photutils','drizzle','stistools','drizzlepac'):
        try:result[name]=version(name)
        except PackageNotFoundError:pass
    return result

def write_manifest(path,**records):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(dict(versions=versions(),**records),indent=2,default=str))
