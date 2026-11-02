from __future__ import annotations
"""Content-addressed stage receipts. No automatic deletion of observations."""
import hashlib,json
from pathlib import Path
from .provenance import identity

class StageCache:
    def __init__(self,directory,prefix=''):
        self.directory=Path(directory);self.prefix=prefix
    def receipt(self,stage):return self.directory/(self.prefix+stage+'.json')
    def key(self,inputs,parameters):
        payload={'inputs':[identity(p) for p in inputs],'parameters':parameters}
        return hashlib.sha256(json.dumps(payload,sort_keys=True,default=str).encode()).hexdigest()
    def valid(self,stage,key):
        p=self.receipt(stage)
        if not p.exists():return False
        data=json.loads(p.read_text())
        return data['key']==key and bool(data['outputs']) and all(Path(a['path']).is_file() and identity(a['path'])==a for a in data['outputs'])
    def record(self,stage,key,outputs):
        self.directory.mkdir(parents=True,exist_ok=True)
        self.receipt(stage).write_text(json.dumps({'key':key,'outputs':[identity(p) for p in outputs]},indent=2))
