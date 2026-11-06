"""Validate selected drizzled products and their recorded input identities."""
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from .config import read_records
from .cli import save


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(c):
    root = c['workdir'] / 'drizzled'
    checks = []
    for record, kind in read_records(c):
        base = f"{record['target']}_{record['filter']}_{kind}"
        provenance = json.loads((root / 'products' / record['target'] / (base + '_provenance.json')).read_text())
        reference = c['manifest'].parent / record['products'][kind]['path']
        if sha(reference) != record['products'][kind]['cutout']['sha256']:
            raise ValueError(f'Science cutout changed: {reference}')
        replay = provenance['validation']
        if (replay['context_equal_fraction'] != 1 or replay['weight_relative_core_max'] > .005
                or replay['weight_relative_median'] > .002):
            raise ValueError(f'Failed drizzle replay: {base}')
        outputs = {r['model']: r for r in provenance['outputs']}
        for model in c['models']:
            if model not in outputs:
                # Legacy provenance predates skipped_models. Determine availability
                # from the configured exact-filter library, never substitute filters.
                library = c['anderson_library'] / f"STDPSF_WFC3UV_{record['filter']}.fits"
                if model == 'anderson' and len(c['models']) > 1 and not library.is_file():
                    checks.append(dict(product=base, model=model, status='skipped_no_library'))
                    continue
                raise ValueError(f'Missing {model} PSF for {base}')
            result = outputs[model]
            path = root / result['path']
            if sha(path) != result['sha256'] or sha(result['model_source']) != result['model_sha256']:
                raise ValueError(f'PSF or source model changed: {path}')
            with fits.open(path, checksum=True) as hdus:
                if any(h.verify_checksum() != 1 or h.verify_datasum() != 1 for h in hdus):
                    raise ValueError(f'Invalid FITS checksum: {path}')
                psf = hdus[0].data
                if psf.shape != (61, 61) or not np.isfinite(psf).all() or not np.isclose(psf.sum(), 1, atol=1e-6):
                    raise ValueError(f'Invalid normalized PSF: {path}')
                if any(hdus[name].data.shape != psf.shape for name in ('WHT', 'CTX')):
                    raise ValueError(f'Invalid WHT/CTX shape: {path}')
            checks.append(dict(product=base, model=model, status='passed', path=str(path)))
    save(root / 'validation.json', checks)
    (root / 'README.md').write_text(
        '# hstpipe UVIS2 PSFs\n\n' +
        '\n'.join(f"- {r['product']} / {r['model']}: {r['status']}" for r in checks) +
        '\n\nUnit-sum finite 61×61 stamps. Use PSFCENX/Y to retain fractional phase. '
        'WHT describes exposure weighting, not PSF uncertainty. Focus, epoch and color '
        'are unmatched. Native empirical stacks contain registration blur and are not '
        'oversampled ePSFs. See docs/uvis-psf.md for method and limitations.\n')
    return checks
