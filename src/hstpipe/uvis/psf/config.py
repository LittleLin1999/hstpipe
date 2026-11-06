"""Portable configuration for the archival UVIS2 workflow."""
import json
import math
import re
from pathlib import Path
import yaml


def read_config(filename):
    filename = Path(filename).expanduser().resolve()
    c = yaml.safe_load(filename.read_text())
    if not isinstance(c, dict):
        raise ValueError('PSF configuration must be a mapping')
    allowed = {'workdir', 'targets', 'manifest', 'raw_root', 'anderson_library',
               'filters', 'products', 'models', 'goal', 'limit', 'chip', 'review_sample'}
    if set(c) - allowed:
        raise ValueError(f'Unknown UVIS PSF options: {sorted(set(c) - allowed)}')
    if not c.get('workdir'):
        raise ValueError('workdir is required')
    if c.get('chip', 2) != 2:
        raise ValueError('This workflow currently supports UVIS2 only')
    def path(value, base=filename.parent):
        value = Path(value).expanduser()
        return (base / value).resolve() if not value.is_absolute() else value.resolve()
    c['workdir'] = path(c['workdir'])
    defaults = dict(targets=c['workdir'] / 'targets.json',
                    review_sample=c['workdir'] / 'review_sample_manifest.json',
                    manifest=c['workdir'] / 'drizzled/input_products/manifest.json',
                    raw_root=c['workdir'].parent / 'raw/uvis',
                    anderson_library=c['workdir'] / 'anderson_comparison/library')
    for key, default in defaults.items():
        c[key] = path(c.get(key, default))
    for key, choices, default in (
        ('products', {'hstpipe', 'mast'}, ['hstpipe']),
        ('models', {'empirical_v2', 'anderson'}, ['empirical_v2']),
    ):
        value = c.get(key, default)
        if not isinstance(value, list) or not value or any(v not in choices for v in value):
            raise ValueError(f'{key} must be a nonempty list drawn from {sorted(choices)}')
        c[key] = list(dict.fromkeys(value))
    filters = c.get('filters', [])
    if not isinstance(filters, list) or any(not isinstance(f, str) or not re.fullmatch(r'F\d{3}[A-Z]', f) for f in filters):
        raise ValueError('filters must be a list of exact filter names, e.g. [F275W]')
    c['filters'] = list(dict.fromkeys(filters))
    for key, default in [('goal', 100), ('limit', 800)]:
        value = c.get(key, default)
        if type(value) is not int or value < 10:
            raise ValueError(f'{key} must be an integer >= 10')
        c[key] = value
    if c['goal'] > c['limit']:
        raise ValueError('goal must not exceed limit')
    return c


def read_targets(c):
    targets = json.loads(c['targets'].read_text())
    if not isinstance(targets, list) or not targets:
        raise ValueError('targets must contain a nonempty JSON list')
    for target in targets:
        if not re.fullmatch(r'F\d{3}[A-Z]', target['filter']) or not target.get('exposures'):
            raise ValueError('Each target requires an exact filter and detector positions')
        for exposure in target['exposures']:
            if exposure.get('chip', 2) != 2:
                raise ValueError('Only UVIS2 detector positions are supported')
            if not all(math.isfinite(float(exposure[k])) for k in ('detector_x', 'detector_y')):
                raise ValueError('Nonfinite target detector position')
    filters = c['filters'] or sorted({t['filter'] for t in targets})
    missing = set(filters) - {t['filter'] for t in targets}
    if missing:
        raise ValueError(f'No detector positions for filters {sorted(missing)}')
    return targets, filters


def read_records(c):
    manifest = json.loads(c['manifest'].read_text())
    records = []
    for record in manifest['records']:
        if c['filters'] and record['filter'] not in c['filters']:
            continue
        if (not re.fullmatch(r'F\d{3}[A-Z]', record['filter'])
                or not re.fullmatch(r'[A-Za-z0-9_+.\-]+', record['target'])
                or record['target'] in ('.', '..')):
            raise ValueError('Invalid target or filter in science manifest')
        for kind in c['products']:
            if kind in record['products']:
                records.append((record, kind))
    if not records:
        raise ValueError('No science products match the requested filters/products')
    missing = set(c['filters']) - {r['filter'] for r, _ in records}
    if missing:
        raise ValueError(f'No matching science products for {sorted(missing)}')
    return records
