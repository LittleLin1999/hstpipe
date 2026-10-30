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

    def __post_init__(self):
        self.input=str(Path(self.input).expanduser().resolve())
        self.output=str(Path(self.output).expanduser().resolve())
        if self.instrument not in ('stis','uvis'):raise ValueError('instrument must be stis or uvis')
        if Path(self.output).is_relative_to(Path(self.input)) or Path(self.input).is_relative_to(Path(self.output)):
            raise ValueError('Raw and reduction directories must be disjoint')

    @classmethod
    def read(cls,path):
        path=Path(path).resolve();data=yaml.safe_load(path.read_text())
        for key in ('input','output'):
            p=Path(data[key]).expanduser();data[key]=str(p if p.is_absolute() else path.parent/p)
        # Calibration paths are relative to YAML, never the installed package.
        for key in ('reference_dir','control_dir','software_dir','source_model'):
            if key in data.get('options',{}):
                p=Path(data['options'][key]).expanduser()
                data['options'][key]=str(p if p.is_absolute() else path.parent/p)
        return cls(**data)

    def effective(self):return asdict(self)
