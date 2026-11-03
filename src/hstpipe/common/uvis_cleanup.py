"""Remove regenerable UVIS comparison images, never science or review inputs."""
from contextlib import ExitStack
from datetime import datetime, timezone
import json
import os
from pathlib import Path

from .paths import ReductionPaths
from .provenance import identity
from .target_runner import save, target_lock

# Exact group-qualified names only. Provisional, review, CR and working files
# deliberately have no matching rule. No age-based or recursive directory delete.
COMPARISONS = ('joint_existing_inspection', 'joint_candidate_inspection',
               'relative_inspection')


def _groups(root, target):
    manifests = ([root / 'manifest.json'] if root.name == 'uvis' else
                 sorted(root.rglob('uvis/manifest.json')))
    groups = []
    for manifest in manifests:
        if not manifest.is_file() or manifest.is_symlink():
            continue
        if any(p.is_symlink() for p in [manifest.parent, *manifest.parent.parents]
               if p.is_relative_to(root)):
            continue
        state = json.loads(manifest.read_text())
        if state.get('schema_version') != 2 or (target and state.get('target') != target):
            continue
        for key, record in state.get('filters', {}).items():
            paths = ReductionPaths(manifest.parent.parent.parent, state['target'],
                                   'uvis', record['filter'], record.get('setup'))
            if paths.directory != manifest.parent or paths.group_key != key:
                raise ValueError(f'Inconsistent group path in {manifest}: {key}')
            groups.append((paths, record))
    return groups


def _signature(paths, record):
    inputs = record.get('inputs', [])
    if not inputs or any(not v.get('sha256') for v in inputs):
        return None
    return (paths.target, paths.group_key,
            tuple(sorted((v['path'], v['sha256'], v['size']) for v in inputs)))


def _plan(root, groups):
    published = {}
    for paths, record in groups:
        if record.get('status') != 'complete' or record.get('scientific_assessment', {}).get('status') != 'science_ready':
            continue
        product = paths.product('drc')
        # A report-only view referencing another reduction is not a publisher.
        expected = next((v for v in record.get('products', [])
                         if v.get('path') == str(product)), None)
        if expected is None or not product.is_file() or product.is_symlink():
            continue
        if any(p.is_symlink() for p in product.parents if p.is_relative_to(root)):
            continue
        if identity(product) != expected:
            raise ValueError(f'Published product differs from its manifest: {product}')
        signature = _signature(paths, record)
        if signature:
            published[signature] = expected
    selected, preserved = [], []
    for paths, record in groups:
        proof = published.get(_signature(paths, record))
        for kind in COMPARISONS:
            path = paths.file('intermediate', 'drizzle', kind, 'fits')
            if not path.is_file() or path.is_symlink():
                continue
            if path.resolve() != path or not path.is_relative_to(root):
                raise ValueError(f'Comparison file escapes reduction directory: {path}')
            if not proof:
                preserved.append(dict(path=str(path), reason='No verified publication for this group and exact input set'))
                continue
            selected.append(dict(**identity(path), published_product=proof['path'],
                                 target=paths.target, filter=paths.group_key,
                                 directory=str(paths.directory)))
    names = {v['path'] for v in selected}
    aliases = []
    for directory, dirs, files in os.walk(root, followlinks=False):
        for name in dirs + files:
            path = Path(directory) / name
            if path.is_symlink() and str(path.resolve()) in names:
                if path.parent.name != 'drizzle' or path.parent.parent.name != 'intermediate':
                    raise ValueError(f'Unexpected alias to a comparison file: {path}')
                aliases.append(str(path))
    receipts = []
    for paths, _ in groups:
        for receipt in paths.cache_directory.glob(paths.group_id + '_cache_*.json'):
            if receipt.is_symlink():
                continue
            outputs = json.loads(receipt.read_text()).get('outputs', [])
            if any(v.get('path') in names for v in outputs):
                receipts.append(str(receipt))
    return dict(files=selected, bytes=sum(v['size'] for v in selected),
                aliases=sorted(set(aliases)), preserved=preserved,
                affected_cache_receipts=sorted(set(receipts)))


def cleanup_uvis(output, *, target=None, apply=False):
    """Preview by default; apply removes only published groups' comparison FITS.

    `output` may be a redux root, a target's uvis directory, or a run containing
    retry redux roots. Publication is verified against the exact input set and
    product checksum. Scientific cache receipts are not rewritten: missing
    diagnostic outputs may require regeneration on a later resume.
    """
    root = Path(output).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f'Reduction directory does not exist: {root}')
    groups = _groups(root, target)
    if not groups:
        raise ValueError(f'No matching UVIS target manifests under {root}')
    with ExitStack() as locks:
        if apply:
            # Same locks/order as reduction and human-review publication.
            for paths, _ in groups:
                locks.enter_context(target_lock(paths.directory, paths.group_key))
            for directory in sorted({paths.directory for paths, _ in groups}):
                locks.enter_context(target_lock(directory))
            groups = _groups(root, target)
        plan = _plan(root, groups)
        result = dict(mode='apply' if apply else 'preview', root=str(root), **plan,
                      removed=[], receipts=[],
                      cache_note='Scientific receipts are unchanged. Receipts listing deleted diagnostics will fail validation until those outputs are regenerated; other stage receipts are retained.')
        if not apply or not plan['files']:
            return result
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        result['timestamp_utc'] = stamp
        for directory in sorted({v['directory'] for v in plan['files']}):
            paths = next(p for p, _ in groups if str(p.directory) == directory)
            receipt = paths.file('logs', None, 'cleanup_' + stamp)
            result['receipts'].append(str(receipt))
        # Persist the full deletion plan before unlinking anything. A partial
        # failure remains recoverable as a recorded operation, not silent loss.
        def checkpoint():
            for receipt in result['receipts']:
                save(receipt, result)
        checkpoint()
        try:
            for item in plan['files']:
                path = Path(item['path'])
                if path.is_symlink() or identity(path) != {k: item[k] for k in ('path', 'size', 'sha256')}:
                    raise ValueError(f'Comparison changed during cleanup: {path}')
                path.unlink()
                result['removed'].append(str(path))
            for alias in plan['aliases']:
                path = Path(alias)
                if path.is_symlink() and str(path.resolve()) in result['removed']:
                    path.unlink()
        except Exception as error:
            result['error'] = str(error)
            raise
        finally:
            checkpoint()
        return result
