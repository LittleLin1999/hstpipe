"""Explicit UVIS PSF stages; independent of science reduction."""
import hashlib
import json
import platform
from contextlib import contextmanager
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from .config import read_config, read_targets, read_records

STAGES = ('all', 'pilot', 'empirical', 'drizzle', 'validate')


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, default=str, allow_nan=False) + '\n')
    temporary.replace(path)


@contextmanager
def lock(root):
    import fcntl
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.hstpipe-uvis-psf.lock').open('a') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Another PSF run is writing {root}') from exc
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def run(config_path, stage='all'):
    if stage not in STAGES:
        raise ValueError(f'Unknown PSF stage: {stage}')
    c = read_config(config_path)
    with lock(c['workdir']):
        receipt = dict(stage=stage, config=c, started=datetime.now(timezone.utc).isoformat(),
                       status='running', python=platform.python_version(), packages={},
                       modules_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                       for p in Path(__file__).parent.glob('*.py')})
        for package in ('hstpipe', 'numpy', 'astropy', 'scipy', 'photutils', 'drizzlepac', 'stwcs'):
            try:
                receipt['packages'][package] = version(package)
            except PackageNotFoundError:
                receipt['packages'][package] = None
        receipt['inputs_sha256'] = {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(config_path).resolve(), c['targets'], c['manifest'], c['review_sample']) if p.is_file()
        }
        receipt_path = c['workdir'] / 'runs' / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%f') + '.json')
        save(receipt_path, receipt)
        try:
            selected = ('empirical', 'drizzle', 'validate') if stage == 'all' else (stage,)
            if stage == 'all':
                read_records(c)  # Fail before downloading if the science selection is unusable.
            for name in selected:
                if name == 'pilot':
                    pilot(c)
                elif name == 'empirical' and 'empirical_v2' in c['models']:
                    empirical(c)
                elif name == 'empirical':
                    if stage != 'all':
                        raise ValueError('empirical stage requires empirical_v2 in models')
                elif name == 'drizzle':
                    drizzle(c)
                elif name == 'validate':
                    from .validation import validate
                    validate(c)
            receipt['status'] = 'complete'
        except Exception as exc:
            receipt.update(status='failed', error=f'{type(exc).__name__}: {exc}')
            raise
        finally:
            receipt['finished'] = datetime.now(timezone.utc).isoformat()
            save(receipt_path, receipt)
    return receipt


def empirical(c):
    from .empirical import EmpiricalBuilder
    targets, filters = read_targets(c)
    builder = EmpiricalBuilder(c['workdir'], targets)
    summary = []
    for band in filters:
        try:
            result = builder.run(band, c['limit'], c['goal'])
        except Exception as exc:
            result = dict(filter=band, status='failed', error=str(exc))
        summary.append(result)
        save(builder.out / 'summary.json', summary)
    failures = [r for r in summary if r['status'] not in ('trial_stack', 'limited_sample_trial')]
    if failures:
        raise RuntimeError(f'Empirical PSF generation failed: {failures}')
    return summary


def drizzle(c):
    from .drizzle import DrizzleBuilder
    records = read_records(c)
    builder = DrizzleBuilder(c['workdir'], c['manifest'], c['raw_root'],
                             c['anderson_library'], c['models'])
    summary = []
    for record, kind in records:
        try:
            result = builder.process(record, kind)
            result['status'] = 'complete'
        except Exception as exc:
            result = dict(target=record['target'], filter=record['filter'], product=kind,
                          status='failed', error=str(exc))
        summary.append(result)
        save(builder.out / 'summary.json', summary)
    failures = [r for r in summary if r['status'] == 'failed']
    if failures:
        raise RuntimeError(f'Drizzled PSF generation failed: {failures}')
    return summary


def pilot(c):
    """Process the saved review sample without changing any stacked PSF."""
    from concurrent.futures import ThreadPoolExecutor
    from .empirical import EmpiricalBuilder
    rows = json.loads(c['review_sample'].read_text())
    if not isinstance(rows, list) or not rows:
        raise ValueError('review_sample must contain a nonempty list of catalog records')
    if c['filters']:
        rows = [r for r in rows if r['filter'] in c['filters']]
    if not rows:
        raise ValueError('No review stars match the requested filters')
    builder = EmpiricalBuilder(c['workdir'], [])
    with ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(builder.process, rows))
    save(builder.out / 'pilot.json', results)
    errors = [r for r in results if r.get('reason') == 'processing_error']
    if errors:
        raise RuntimeError(f'Review sample processing failed for {len(errors)} stars')
    return results
