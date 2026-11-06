"""Replay retained science geometry and weights for target-position UVIS2 PSFs."""
import json, re, hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from stwcs.wcsutil import HSTWCS
from drizzlepac.adrizzle import do_driz
from drizzlepac import cdriz
from photutils.psf import GriddedPSFModel
from stsci.imagestats import ImageStats
from .render import render_empirical


WCSKEY = re.compile('^(CRPIX[12]|CRVAL[12]|CD[12]_[12]|CTYPE[12]|CUNIT[12]|WCSAXES|WCSNAME|[AB]P?_.*|IDCSCALE|LTV[12]|LTM[12]_[12]|ORIENTAT)$')

def save(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def arrays():
    return (np.zeros((61, 61), np.float32), np.zeros((61, 61), np.float32), np.zeros((61, 61), np.int32))

class FullDetectorMap:

    def __init__(self, w, outwcs):
        if not hasattr(w, '_fullmap'):
            full = w._full_wcs
            w._fullmap = cdriz.DefaultWCSMapping(full, outwcs, full.pixel_shape[0], full.pixel_shape[1], 10)
        self.mapping = w._fullmap
        self.origin = w._input_origin

    def forward(self, x, y):
        return self.mapping(x + self.origin[0], y + self.origin[1])

def drizzle(a, w, weight, outwcs, outputs, d):
    do_driz(np.asarray(a, dtype=np.float32), w, np.asarray(weight, dtype=np.float32), outwcs, *outputs, expin=d['exptime'], in_units='cps', wt_scl=d['weight_scale'], wcslin_pscale=d['input_scale'], uniqid=d['drizzle_id'], pixfrac=d['pixfrac'], kernel=d['kernel'], fillval='0', stepsize=10, wcsmap=FullDetectorMap)

class DrizzleBuilder:
    def __init__(self, root, manifest, raw_root, library, models):
        self.root = Path(root)
        self.out = self.root / 'drizzled'
        self.out.mkdir(parents=True, exist_ok=True)
        self.input = Path(manifest).parent
        self.raw_root = Path(raw_root)
        self.library = Path(library)
        self.models = list(models)

    def static_masks(self, record, kind, header, table, working):
        cache = self.out / 'static_masks' / f"{record['target']}_{record['filter']}_{kind}"
        cache.parent.mkdir(exist_ok=True)
        masks = {}
        stats = []
        sources = []
        for i in range(1, int(header['NDRIZIM']) + 1):
            name, ext = re.match('(.*?)\\[sci,(\\d+)\\]', header[f'D{i:03d}DATA'], re.I).groups()
            ext = int(ext)
            name = Path(name).name
            row = next((r for r in table if r['ROOTNAME'].strip() == name.split('_')[0] and int(r['EXTVER']) == ext))
            if int(row['CCDCHIP']) != 2:
                continue
            source = working / name if kind == 'hstpipe' else self.raw_root / name
            sources.append((source, ext))
        for source, ext in sources:
            with fits.open(source, memmap=False) as h:
                a = h['SCI', ext].data
                stat = ImageStats(a, nclip=3, fields='mode', lower=np.nanmin(a), upper=np.nanmax(a))
                mode, rms, nbins = (float(stat.mode), float(stat.stddev), len(stat.histogram))
                if a.shape not in masks:
                    masks[a.shape] = np.ones(a.shape, bool)
                if nbins >= 2:
                    masks[a.shape] &= ~(a < mode - 4 * rms)
                stats.append(dict(source=str(source), extver=ext, mode=mode, rms=rms, nbins=nbins, static_sigma=4))
        np.savez_compressed(Path(str(cache) + '.npz'), **{'x'.join(map(str, k)): v for k, v in masks.items()})
        save(Path(str(cache) + '.json'), stats)
        return masks

    def setup(self, record, kind):
        product = record['products'][kind]
        path = self.input / product['path']
        if sha(path) != product['cutout']['sha256']:
            raise ValueError(f'Science cutout checksum mismatch: {path}')
        with fits.open(path) as h:
            header = h[0].header.copy()
            if int(header['NDRIZIM']) > 32:
                raise ValueError('Multi-plane CTX (>32 inputs) is not supported')
            output = HSTWCS(h, ext=('SCI', 1))
            x, y = output.all_world2pix([record['center_ra_deg']], [record['center_dec_deg']], 0)
            x, y = (float(np.ravel(x)[0]), float(np.ravel(y)[0]))
            ox = int(round(x)) - 30
            oy = int(round(y)) - 30
            if ox < 0 or oy < 0 or ox + 61 > h['SCI'].data.shape[-1] or (oy + 61 > h['SCI'].data.shape[-2]):
                raise ValueError('Science cutout must contain the full 61x61 PSF stamp')
            output = output.slice((slice(oy, oy + 61), slice(ox, ox + 61)))
            reference = {name: h[name].data[oy:oy + 61, ox:ox + 61].copy() for name in ['SCI', 'WHT', 'CTX']}
            table = h['HDRTAB'].data.copy()
        working = Path(product.get('working_directory', Path(product['source_resolved']).parent.parent / 'intermediate' / 'working'))
        static = self.static_masks(record, kind, header, table, working)
        data = []
        for i in range(1, int(header['NDRIZIM']) + 1):
            prefix = f'D{i:03d}'
            name, ext = re.match('(.*?)\\[sci,(\\d+)\\]', header[prefix + 'DATA'], re.I).groups()
            ext = int(ext)
            name = Path(name).name
            root = name.split('_')[0]
            row = next((r for r in table if r['ROOTNAME'].strip() == root and int(r['EXTVER']) == ext))
            if int(row['CCDCHIP']) != 2:
                continue
            working = Path(product.get('working_directory', Path(product['source_resolved']).parent.parent / 'intermediate' / 'working'))
            source = working / name if kind == 'hstpipe' else self.raw_root / name
            if not source.exists():
                raise FileNotFoundError(source)
            with fits.open(source, memmap=False) as h:
                sh = h['SCI', ext].header
                for k in table.names:
                    if WCSKEY.match(k):
                        sh[k] = row[k].item() if isinstance(row[k], np.generic) else row[k]
                w = HSTWCS(h, ext=('SCI', ext))
                tx, ty = w.all_world2pix([record['center_ra_deg']], [record['center_dec_deg']], 0)
                tx, ty = (float(np.ravel(tx)[0]), float(np.ravel(ty)[0]))
                ix = int(round(tx)) - 55
                iy = int(round(ty)) - 55
                ny, nx = h['SCI', ext].data.shape
                ix = max(0, ix)
                iy = max(0, iy)
                ex = min(nx, ix + 111)
                ey = min(ny, iy + 111)
                if not (0 <= tx < nx and 0 <= ty < ny):
                    raise ValueError('target outside expected UVIS2 input')
                sci = h['SCI', ext].data[iy:ey, ix:ex].astype(float) - float(row['MDRIZSKY'])
                dq = h['DQ', ext].data[iy:ey, ix:ex].copy()
                exp = float(header[prefix + 'DEXP'])
                maskpath = working / Path(header[prefix + 'MASK']).name
                if kind == 'hstpipe':
                    maskpath = working.parent / 'drizzle' / f'{root}_final_crmask.fits'
                    badcr = fits.getdata(maskpath, ext=('SCI', ext))[iy:ey, ix:ex].astype(bool)
                    bits = int(header['HPDQSCI'])
                    weight = (dq.astype(np.int32) & (65535 ^ (bits | 4096)) == 0) & ~badcr
                    weight = weight.astype(float) * exp
                    method = 'retained_final_CR_plus_calibration_DQ_plus_rebuilt_static'
                else:
                    weight = (dq.astype(np.int32) & (65535 ^ 336) == 0).astype(float) * exp
                    method = 'raw_FLC_DQ_ignore336_plus_rebuilt_static_times_exptime'
                weight *= static[ny, nx][iy:ey, ix:ex]
                sliced = w.slice((slice(iy, ey), slice(ix, ex)))
                sliced._full_wcs = w
                sliced._input_origin = (ix, iy)
                a = np.array([0.0, 30.0, 60.0])
                b = np.array([0.0, 30.0, 60.0])
                wa = w.all_pix2world(a + ix, b + iy, 0)
                wb = sliced.all_pix2world(a, b, 0)
                if np.max(abs(np.array(wa) - np.array(wb))) >= 1e-09:
                    raise ValueError('Sliced input WCS differs from the full-detector WCS')
                detector_x = tx - float(sh.get('LTV1', 0))
                detector_y = ty - float(sh.get('LTV2', 0))
                d = dict(root=root, extver=ext, chip=2, drizzle_id=i, source=str(source), weight_file=str(maskpath) if kind == 'hstpipe' else None, weight_method=method, input_xy=[tx, ty], detector_xy=[detector_x, detector_y], input_origin_xy=[ix, iy], output_xy=[x - ox, y - oy], exptime=exp, sky_subtracted=float(row['MDRIZSKY']), pixfrac=float(header[prefix + 'PIXF']), kernel=header[prefix + 'KERN'].strip(), input_scale=float(header[prefix + 'ISCL']), output_scale=float(header[prefix + 'SCAL']), weight_scale=float(header[prefix + 'WTSC']))
                data.append((d, sliced, sci, weight))
        if not data:
            raise ValueError('No contributing UVIS2 exposures')
        return (path, header, output, reference, data, [ox, oy], [x - ox, y - oy])

    def process(self, record, kind, probe=False):
        path, header, outwcs, reference, data, origin, center = self.setup(record, kind)
        target, band = (record['target'], record['filter'])
        folder = self.out / 'products' / target
        folder.mkdir(parents=True, exist_ok=True)
        base = f'{target}_{band}_{kind}'
        replay = arrays()
        for d, w, sci, weight in data:
            drizzle(sci / d['exptime'], w, weight, outwcs, replay, d)
        good = reference['WHT'] > 0
        rel = (replay[1] - reference['WHT']) / np.maximum(reference['WHT'], 1)
        yy, xx = np.indices((61, 61))
        rad = np.hypot(xx - center[0], yy - center[1])
        core = rad <= 6
        validation = dict(weight_relative_median=float(np.median(abs(rel[good]))), weight_relative_max=float(np.max(abs(rel[good]))), weight_relative_core_max=float(np.max(abs(rel[good & core]))), context_equal_fraction=float(np.mean(replay[2] == reference['CTX'])), science_replay_median_offset=float(np.median(replay[0] - reference['SCI'])), science_replay_mad_after_offset=float(np.median(abs(replay[0] - reference['SCI'] - np.median(replay[0] - reference['SCI'])))), science_replay_median_abs=float(np.median(abs(replay[0] - reference['SCI']))), science_replay_max_abs=float(np.max(abs(replay[0] - reference['SCI']))))
        summary = dict(target=target, filter=band, product=kind, reference=str(path.relative_to(self.input)), output_origin_in_reference_xy=origin, psf_center_zero_based=center, input_exposures=len(data), pixfrac=data[0][0]['pixfrac'], kernel=data[0][0]['kernel'], pixel_scale_arcsec=data[0][0]['output_scale'], validation=validation, inputs=[d for d, *_ in data], outputs=[], focus_matched=False)
        save(folder / (base + '_provenance.json'), summary)
        if probe:
            fits.HDUList([fits.PrimaryHDU(replay[0]), fits.ImageHDU(replay[1], name='WHT'), fits.ImageHDU(reference['WHT'], name='REFWHT'), fits.ImageHDU(replay[2], name='CTX'), fits.ImageHDU(reference['CTX'], name='REFCTX')]).writeto(folder / (base + '_replay.fits'), overwrite=True)
            print('PROBE', base, json.dumps(validation), flush=True)
            return summary
        if validation['context_equal_fraction'] != 1:
            raise ValueError(f'Context replay mismatch: {validation}')
        if validation['weight_relative_core_max'] > 0.005 or validation['weight_relative_median'] > 0.002:
            raise ValueError(f'weight replay mismatch: {validation}')
        for method in self.models:
            lib = self.library / f'STDPSF_WFC3UV_{band}.fits'
            if method == 'anderson' and (not lib.exists()):
                if self.models == ['anderson']:
                    raise FileNotFoundError(lib)
                summary.setdefault('skipped_models', []).append(dict(model=method, reason='No exact-filter library file', path=str(lib)))
                continue
            if method == 'anderson':
                model = GriddedPSFModel.read(lib, format='stdpsf', detector_id=1)
                modelpath = lib
            else:
                modelpath = self.root / 'refined_v2' / 'psfs' / f'MAST_UVIS2_{band}_v2_trial.fits'
                model = fits.getdata(modelpath).astype(float)
            result = arrays()
            input_sums = []
            for d, w, sci, weight in data:
                iy, ix = np.indices(sci.shape)
                ix = ix + d['input_origin_xy'][0]
                iy = iy + d['input_origin_xy'][1]
                dx = ix - d['input_xy'][0]
                dy = iy - d['input_xy'][1]
                if method == 'anderson':
                    px, py = d['detector_xy']
                    render = model.evaluate(px + dx, py + dy, 1.0, px, py)
                else:
                    render = render_empirical(model, dx, dy)
                total = float(render.sum())
                if total <= 0 or not np.isfinite(render).all():
                    raise ValueError('Invalid rendered PSF')
                input_sums.append(total)
                render /= total
                drizzle(render, w, weight, outwcs, result, d)
            total = float(result[0].sum())
            if not np.isfinite(total) or total <= 0:
                raise ValueError('No positive drizzled PSF flux')
            psf = result[0] / total
            if not np.isfinite(psf).all() or abs(float(psf.sum()) - 1) > 1e-06:
                raise ValueError('Invalid normalized PSF')
            hdr = outwcs.to_header(relax=True)
            for k, v in dict(FILTER=band, PSFMODEL=method, PIXFRAC=data[0][0]['pixfrac'], KERNEL=data[0][0]['kernel'], PIXSCALE=data[0][0]['output_scale'], NEXP=len(data), PSFCENX=center[0], PSFCENY=center[1], SCICUTX=origin[0], SCICUTY=origin[1], DRIZMAT=True, FOCMATCH=False, BUNIT='FRACTION').items():
                hdr[k] = v
            hdr['HISTORY'] = 'Actual exposure WCS from science HDRTAB, distortion tables, dithers and final weights.'
            hdr['HISTORY'] = 'Normalized finite 61x61 stamp; PSFCENX/Y zero based, fractional phase retained.'
            if method == 'empirical_v2':
                hdr['HISTORY'] = 'Native v2 stack already includes registration blur; cubic phase interpolation is approximate.'
            outpath = folder / (base + f'_psf_{method}.fits')
            fits.HDUList([fits.PrimaryHDU(psf.astype('f4'), hdr), fits.ImageHDU(result[1], name='WHT'), fits.ImageHDU(result[2], name='CTX')]).writeto(outpath, overwrite=True, checksum=True)
            summary['outputs'].append(dict(model=method, path=str(outpath.relative_to(self.out)), model_source=str(modelpath), model_sha256=sha(modelpath), input_finite_sums=input_sums, output_sum_before_normalization=total, sha256=sha(outpath)))
        save(folder / (base + '_provenance.json'), summary)
        print('DONE', base, json.dumps(validation), flush=True)
        return summary
