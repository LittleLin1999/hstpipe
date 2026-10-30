from __future__ import annotations
from dataclasses import dataclass
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
    @property
    def directory(self):
        p=Path(self.root)/label(self.target)/label(self.instrument)/label(self.element)
        return p/label(self.setup) if self.setup else p
    def area(self,name,stage=None,create=False):
        if name not in ('products','intermediate','diagnostics','logs'):raise ValueError(name)
        p=self.directory/name
        if stage:p=p/stage
        if create:p.mkdir(parents=True,exist_ok=True)
        return p
    def product(self,suffix):return self.area('products')/f'{label(self.target)}_{label(self.element)}_{suffix}.fits'
