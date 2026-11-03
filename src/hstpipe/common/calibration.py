"""Retrieve exact header-named HST references, without recalibration or bestrefs."""
from datetime import datetime, timezone
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from astropy.io import fits
from .provenance import identity
from .target_runner import save

DEFAULTS = dict(auto_download=True, offline=False, directory=None, timeout=120, retries=1)
KEYS = ('IDCTAB', 'NPOLFILE', 'D2IMFILE', 'CCDTAB', 'MDRIZTAB')


def reference_requests(inputs, amplifier=False, instrument="uvis"):
    requests = []
    for image in inputs:
        header = fits.getheader(image)
        keys = KEYS + (('PFLTFILE',) if amplifier else ())
        if instrument == 'stis':
            keys = [key for key,value in header.items() if isinstance(value,str) and
                    (value.startswith('oref$') or key.endswith(('FILE','TAB')))]
        for key in keys:
            value = str(header.get(key, '')).strip()
            if not value or value.upper() in ('N/A', 'NONE', 'INDEF'):
                continue
            if instrument == 'stis' and not value.lower().endswith('.fits'):
                continue
            requests.append(dict(input=str(image), keyword=key, header_value=value,
                                 required=key != 'MDRIZTAB'))
    return requests


def reference_location(value, directory, prefix_name="iref"):
    if '$' in value:
        prefix, name = value.split('$', 1)
        if prefix != prefix_name:
            raise ValueError(f'Unsupported reference prefix: {prefix}')
        path = directory / name
    elif Path(value).is_absolute():
        path = Path(value)
    else:
        name = value
        path = directory / name
    if not re.fullmatch(r'[A-Za-z0-9_.+-]+\.fits', path.name, flags=re.I):
        raise ValueError(f'Invalid reference filename: {value}')
    if not Path(value).is_absolute() and path.parent != directory:
        raise ValueError(f'Reference must be a plain filename: {value}')
    return path


def _verify(path):
    with fits.open(path, checksum=True) as hd:
        hd.verify('exception')
        for hdu in hd:
            if 'CHECKSUM' in hdu.header and hdu.verify_checksum() != 1:
                raise ValueError(f'Invalid FITS checksum: {path}')
            if 'DATASUM' in hdu.header and hdu.verify_datasum() != 1:
                raise ValueError(f'Invalid FITS data checksum: {path}')


def _fetch(name, directory, timeout):
    env = dict(os.environ, CRDS_SERVER_URL='https://hst-crds.stsci.edu')
    env.setdefault('CRDS_PATH', str(directory.parent.parent / '.crds-cache'))
    command = [sys.executable, '-m', 'crds.sync', '--hst', '--files', name,
               '--output-dir', str(directory)]
    result = subprocess.run(command, env=env, capture_output=True, text=True,
                            timeout=timeout, check=False)
    if result.returncode:
        raise RuntimeError(f'CRDS download failed for {name}: ' +
                           (result.stderr or result.stdout)[-2000:])
    if not (directory / name).is_file():
        raise RuntimeError(f'CRDS did not supply {name}')


def ensure_references(inputs, config, paths):
    """Reuse verified files; atomically install missing exact references under locks."""
    import fcntl
    options = dict(DEFAULTS, **config.options.get('calibration', {}))
    if set(options) - set(DEFAULTS):
        raise ValueError('Unknown calibration settings')
    if any(type(options[k]) is not bool for k in ('auto_download','offline')):
        raise ValueError('Calibration switches must be boolean')
    if not isinstance(options['timeout'],(int,float)) or not 0 < options['timeout'] <= 600:
        raise ValueError('Calibration timeout must be in (0,600] seconds')
    if type(options['retries']) is not int or not 0 <= options['retries'] <= 3:
        raise ValueError('Calibration retries must be an integer from 0 to 3')
    instrument = config.instrument
    prefix = 'iref' if instrument == 'uvis' else 'oref'
    explicit = config.options.get('reference_dir') if instrument == 'stis' else None
    if explicit and options['directory'] and Path(explicit).resolve() != Path(options['directory']).resolve():
        raise ValueError('STIS reference_dir and calibration.directory must agree')
    directory = Path(explicit or options['directory'] or os.environ.get(prefix) or
                     Path(config.output).parent / ('calibration/' + instrument)).expanduser().resolve()
    raw = Path(config.input).resolve()
    if directory.is_relative_to(raw):
        raise ValueError('Calibration downloads cannot be stored in the raw input directory')
    requests = reference_requests(inputs, bool(config.options.get('background', {}).get('amplifier_correction') or
                                               config.options.get('background', {}).get('amplifier_assess')), instrument=instrument)
    report = dict(server='https://hst-crds.stsci.edu', directory=str(directory),
                  timestamp_utc=datetime.now(timezone.utc).isoformat(),
                  selection='Exact filenames from input headers; no bestrefs or header updates',
                  references=[], warnings=[], status='checking')
    receipt = paths.file('logs', None, 'calibration_references')
    os.environ[prefix] = str(directory) + '/'
    grouped = {}
    for row in requests:
        path = reference_location(row['header_value'], directory, prefix)
        if path.resolve().is_relative_to(raw):
            raise ValueError('A calibration reference resolves into raw inputs')
        grouped.setdefault(path, []).append(row)
    try:
        for path, usage in grouped.items():
            item = dict(path=str(path), uses=usage, status='pending')
            report['references'].append(item)
            try:
                if path.is_file():
                    _verify(path)
                    item.update(status='reused', identity=identity(path))
                    continue
                if options['offline'] or not options['auto_download']:
                    raise FileNotFoundError(f'Missing reference {path}; calibration downloads are disabled')
                if path.parent != directory:
                    raise FileNotFoundError(f'Missing absolute reference {path}; restore it at its recorded location')
                directory.mkdir(parents=True, exist_ok=True)
                with (directory / ('.' + path.name + '.lock')).open('a') as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX)
                    if path.is_file():
                        _verify(path)
                        item.update(status='reused', identity=identity(path))
                        continue
                    for attempt in range(options['retries'] + 1):
                        item['attempts'] = attempt + 1
                        try:
                            with tempfile.TemporaryDirectory(prefix='.crds-', dir=directory) as temp:
                                temp = Path(temp)
                                _fetch(path.name, temp, options['timeout'])
                                _verify(temp / path.name)
                                os.replace(temp / path.name, path)
                            break
                        except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired):
                            if attempt == options['retries']:
                                raise
                    item.update(status='downloaded', identity=identity(path))
            except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
                item.update(status='unavailable', error=str(error))
                if any(row['required'] for row in usage):
                    raise RuntimeError(f'Required calibration reference unavailable: {path.name}. {error}') from error
                report['warnings'].append(f'Optional MDRIZTAB unavailable: {path.name}: {error}')
        report['status'] = 'ready' if not report['warnings'] else 'ready_with_warnings'
    except Exception:
        report['status'] = 'failed'
        raise
    finally:
        save(receipt, report)
    return report


def ensure_uvis_references(inputs, config, paths):
    return ensure_references(inputs, config, paths)
