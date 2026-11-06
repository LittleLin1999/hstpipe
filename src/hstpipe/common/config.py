from __future__ import annotations
from dataclasses import dataclass, field, asdict
from pathlib import Path
import yaml

@dataclass
class ReductionConfig:
    input: str
    output: str
    instrument: str
    target: str | None = None
    element: str | None = None
    setup: str | None = None
    options: dict = field(default_factory=dict)
    resume: bool = False
    force: bool = False
    filter_overrides: dict = field(default_factory=dict)
    field: str | None = None

    def __post_init__(self):
        if self.field is not None:
            from .paths import label
            label(self.field)
            if self.instrument!='uvis':raise ValueError('Explicit field grouping is supported for UVIS only')
            if self.target not in (None,self.field):raise ValueError('--target selects the output field; with --field it must match that field')
        self.input=str(Path(self.input).expanduser().resolve())
        self.output=str(Path(self.output).expanduser().resolve())
        if self.instrument not in ('stis','uvis'):raise ValueError('instrument must be stis or uvis')
        if Path(self.output).is_relative_to(Path(self.input)) or Path(self.input).is_relative_to(Path(self.output)):
            raise ValueError('Raw and reduction directories must be disjoint')

    @classmethod
    def read(cls,path):
        path=Path(path).resolve();data=yaml.safe_load(path.read_text())
        if data.get('instrument')=='uvis' and 'astrometry' in data:
            options=data.setdefault('options',{})
            if 'astrometry' in options and options['astrometry']!=data['astrometry']:raise ValueError('Conflicting top-level/options astrometry settings')
            options['astrometry']=data.pop('astrometry')
        for key in ('input','output'):
            p=Path(data[key]).expanduser();data[key]=str(p if p.is_absolute() else path.parent/p)
        # Calibration paths are relative to YAML, never the installed package.
        for key in ('reference_dir','control_dir','software_dir','source_model'):
            if key in data.get('options',{}):
                p=Path(data['options'][key]).expanduser()
                data['options'][key]=str(p if p.is_absolute() else path.parent/p)
        for key in ('background_masks','amplifier_masks'):
            for name,value in data.get('options',{}).get(key,{}).items():
                p=Path(value).expanduser()
                data['options'][key][name]=str(p if p.is_absolute() else path.parent/p)
        if data.get('instrument')=='uvis' and data.get('options',{}).get('alignment',{}).get('refcat'):
            p=Path(data['options']['alignment']['refcat']).expanduser()
            data['options']['alignment']['refcat']=str(p if p.is_absolute() else path.parent/p)
        for options in [data.get('options',{}),*data.get('filter_overrides',{}).values()]:
            value=options.get('calibration',{}).get('directory')
            if value:
                p=Path(value).expanduser();options['calibration']['directory']=str(p if p.is_absolute() else path.parent/p)
        for override in data.get('filter_overrides',{}).values():
            value=override.get('alignment',{}).get('refcat')
            if value:
                p=Path(value).expanduser();override['alignment']['refcat']=str(p if p.is_absolute() else path.parent/p)
        return cls(**data)

    def effective(self):
        data=asdict(self)
        if self.field is None:data.pop('field')
        if self.instrument=='stis':data.pop('filter_overrides',None)
        return data
