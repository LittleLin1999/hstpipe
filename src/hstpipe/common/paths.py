from __future__ import annotations
from dataclasses import dataclass,field
from pathlib import Path
import re

def label(value):
    if not value or not re.fullmatch(r'[A-Za-z0-9_.+-]+',value) or value in ('.','..'):
        raise ValueError(f'Unsafe directory label: {value!r}')
    return value

@dataclass(frozen=True)
class ReductionPaths:
    root: Path
    target: str
    instrument: str
    element: str
    setup: str | None = None
    _files: set = field(default_factory=set,compare=False,repr=False)
    @property
    def directory(self):
        p=Path(self.root)/label(self.target)/label(self.instrument)
        if self.instrument=='uvis':return p
        p=p/label(self.element)
        return p/label(self.setup) if self.setup else p
    @property
    def group_id(self):
        return f'{label(self.target)}_{label(self.element)}'+(f'_{label(self.setup)}' if self.setup else '')
    @property
    def group_key(self):return label(self.element)+(f'_{label(self.setup)}' if self.setup else '')
    @property
    def cache_directory(self):
        if self.instrument=='uvis':return self.area('intermediate','alignment')
        p=self.directory/'.cache'
        return p
    def cache(self):
        from .cache import StageCache
        return StageCache(self.cache_directory,prefix=self.group_id+'_cache_' if self.instrument=='uvis' else '')
    def report(self):
        return self._file(self.area('products')/f'{label(self.target)}_uvis_report.html')
    def scratch(self,prefix='run-'):
        """Ephemeral native-tool workspaces live outside the reduction tree."""
        import tempfile
        return tempfile.TemporaryDirectory(prefix='hstpipe-'+prefix)
    def recovery_directory(self,tag):
        return self.directory.with_name('.uvis-'+label(tag))
    def area(self,name,stage=None,create=False):
        if name not in ('products','intermediate','diagnostics','logs'):raise ValueError(name)
        p=self.directory/name
        if self.instrument=='uvis':
            if name in ('diagnostics','logs'):stage=None
            elif name=='intermediate':
                if stage not in ('catalogs','background','alignment','drizzle','working'):
                    raise ValueError(f'UVIS intermediate stage must be shallow: {stage}')
        if stage:
            stage_path=Path(stage)
            if stage_path.is_absolute():raise ValueError('Stage must be relative')
            p=p.joinpath(*(label(part) for part in stage_path.parts))
        if create:p.mkdir(parents=True,exist_ok=True)
        return p
    def _file(self,path):
        self._files.add(path);return path
    def qualified_filename(self,area,stage,name):
        if self.instrument=='uvis' and area in ('diagnostics','logs') and stage:
            return label(stage)+'_'+name
        return name
    def file(self,area,stage,kind,suffix='json'):
        name=self.qualified_filename(area,stage,f'{self.group_id}_{label(kind)}.{label(suffix)}')
        return self._file(self.area(area,stage)/name)
    def exposure_file(self,area,stage,root,kind,suffix='fits',group=False):
        if self.instrument=='uvis' and area=='intermediate' and kind=='flc' and suffix=='fits':stage='working'
        stem=f'{label(root)}'+(f'_{label(self.setup)}' if self.setup else '')
        if group:stem=f'{self.group_id}_{label(root)}'
        name=self.qualified_filename(area,stage,f'{stem}_{label(kind)}.{label(suffix)}')
        return self._file(self.area(area,stage)/name)
    def product(self,suffix):
        if self.instrument=='uvis':return self.file('products',None,suffix,'fits')
        return self.area('products')/f'{label(self.target)}_{label(self.element)}_{suffix}.fits'
    def detach(self,path):
        """Discard a checkpoint before recreating it, leaving migration backups intact."""
        path=Path(path)
        if path.exists():path.unlink()

    def copy(self,source,destination):
        """Atomically replace an output inode; recovery hard links stay immutable."""
        import os,shutil,tempfile
        destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
        fd,name=tempfile.mkstemp(prefix='.copy-',dir=destination.parent);os.close(fd)
        temp=Path(name)
        try:
            shutil.copy2(source,temp);os.replace(temp,destination)
        finally:
            if temp.exists():temp.unlink()
        self._files.add(destination)
        return destination

    def catalog_directory(self,create=False):
        path=self.area('intermediate','catalogs') if self.instrument=='uvis' else self.area('intermediate','alignment')/'catalogs'
        if create:path.mkdir(parents=True,exist_ok=True)
        return path
