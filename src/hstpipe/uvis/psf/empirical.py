"""MAST UVIS2 empirical v2 algorithm, ported from refine_mast_psfs.py.

Scientific cuts and registration are preserved. Each builder owns its paths
and locks; importing this module never creates files or reads a dataset.
"""
import json, hashlib, io, time, threading, warnings
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import requests
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
from scipy.ndimage import convolve, binary_dilation
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from photutils.segmentation import detect_sources, deblend_sources

VERSION = 'hstpipe-uvis2-v2-1'
COLS = 'id,rootname,filter,chip,psf_x_center,psf_y_center,x_cal,y_cal,qfit,iso_index,pixc,psf_flux,sky,focus,mjd,exptime,subarray,psf_ra,psf_dec'
Y, X = np.indices((51, 51))
R = np.hypot(X - 25, Y - 25)
YY, XX = np.indices((41, 41))
RR = np.hypot(XX - 20, YY - 20)
ky, kx = np.indices((7, 7))
K = np.exp(-((kx - 3) ** 2 + (ky - 3) ** 2) / (2 * 0.8 ** 2))
K /= K.sum()

def disk(radius):
    y, x = np.mgrid[-radius:radius + 1, -radius:radius + 1]
    return x * x + y * y <= radius * radius

def save(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + f'.tmp-{threading.get_ident()}')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)

def rng(k, a, b):
    return {'paramName': k, 'values': [{'min': a, 'max': b}]}

def val(k, v):
    return {'paramName': k, 'values': [str(v)]}

def sky_groups(rows):
    if not rows:
        return np.array([], dtype=int)
    c = SkyCoord([r['psf_ra'] for r in rows] * u.deg, [r['psf_dec'] for r in rows] * u.deg)
    i, j, _, _ = search_around_sky(c, c, 0.3 * u.arcsec)
    _, labels = connected_components(coo_matrix((np.ones(len(i)), (i, j)), shape=(len(c), len(c))), directed=False)
    return labels

def registration(dx, dy):

    def weights(delta):
        t = np.arange(41)[:, None] + 5 - delta - np.arange(51)[None, :]
        w = np.maximum(1 - abs(t), 0)
        w /= w.sum(axis=1, keepdims=True)
        return w
    return (weights(dx), weights(dy))

def segment_mask(b, err, bad, cx, cy, neighbors, row):
    sm = convolve(np.where(bad, 0, b), K, mode='reflect')
    smooth_err = np.sqrt(convolve(np.where(np.isfinite(err), err, 1000000.0) ** 2, K * K, mode='reflect'))
    seg = detect_sources(sm, 2.5 * smooth_err, n_pixels=5, mask=bad)
    labels = np.zeros(b.shape, dtype=int)
    central = 0
    if seg is not None:
        seg = deblend_sources(sm, seg, n_pixels=5, n_levels=64, contrast=0.0005, mode='linear', n_processes=1, progress_bar=False)
        labels = seg.data
        central = int(labels[int(round(cy)), int(round(cx))])
        if central == 0:
            near = labels[np.hypot(X - cx, Y - cy) < 2]
            near = near[near > 0]
            if len(near):
                central = int(Counter(near).most_common(1)[0][0])
    other = (labels > 0) & (labels != central)
    other = binary_dilation(other, structure=disk(3))
    catmask = np.zeros(b.shape, bool)
    cat_count = 0
    noise = max(float(np.median(err[(R > 18) & ~bad])), 1.0)
    for n in neighbors:
        if n['id'] == row.get('id'):
            continue
        if n.get('pixc') is None or n.get('psf_flux') is None or n['psf_flux'] < 300 or (n['pixc'] < max(30, 3 * noise)):
            continue
        nx = cx + n['x_cal'] - row['x_cal']
        ny = cy + n['y_cal'] - row['y_cal']
        dist = np.hypot(nx - cx, ny - cy)
        if dist < 1 or not -12 < nx < 63 or (not -12 < ny < 63):
            continue
        radius = float(np.clip(2 + np.sqrt(2 * 1.2 ** 2 * np.log(max(n['pixc'] / noise, 1))), 4, 10))
        catmask |= np.hypot(X - nx, Y - ny) <= radius
        cat_count += 1
    return (other | catmask, labels, cat_count)

def process_arrays(row, sci, err, dq, neighbors):
    row = dict(row, algorithm_version=VERSION, accepted=False)
    bad = ~np.isfinite(sci) | ~np.isfinite(err) | (err <= 0) | ~np.isfinite(dq) | (dq != 0)
    dqmask = bad.copy()
    audit = {'SCI': sci, 'DQMASK': dqmask.astype('u1')}

    def reject(reason):
        row['reason'] = reason
        return (row, audit)
    if np.any(dqmask[R <= 1.5]):
        return reject('flagged_central_pixels')
    good = (R > 18) & ~dqmask
    if good.sum() < 100:
        return reject('insufficient_background_pixels')
    _, sky, noise = sigma_clipped_stats(sci[good], sigma=3, maxiters=5)
    b = sci - sky
    yy, xx = np.indices((9, 9))
    z = b[21:30, 21:30]
    ze = err[21:30, 21:30]
    ok = ~dqmask[21:30, 21:30]

    def model(p):
        return p[0] * np.exp(-0.5 * (((xx - p[1]) / p[3]) ** 2 + ((yy - p[2]) / p[4]) ** 2)) + p[5]
    fit = least_squares(lambda p: ((model(p) - z) / np.maximum(ze, 1))[ok], [max(float(z.max()), 1), 4, 4, 1, 1, 0], bounds=([0, 2, 2, 0.4, 0.4, -10000], [10000000.0, 6, 6, 3, 3, 10000]))
    amp, cx0, cy0, sx, sy, bg = fit.x
    cx = 21 + cx0
    cy = 21 + cy0
    dx = 25 - cx
    dy = 25 - cy
    row.update(shift_x=float(dx), shift_y=float(dy), fwhm=float(2.355 * np.sqrt(sx * sy)), axis_ratio=float(max(sx, sy) / min(sx, sy)))
    if abs(dx) > 1.5 or abs(dy) > 1.5 or row['axis_ratio'] > 1.6 or (not 0.9 < row['fwhm'] < 3.5):
        return reject('core_shape_or_center')
    neighbor, labels, ncat = segment_mask(b, err, dqmask, cx, cy, neighbors, row)
    mask = neighbor | dqmask
    background = (R > 16) & ~mask
    if background.sum() < 100:
        audit['MASKRAW'] = mask.astype('u1')
        return reject('insufficient_clean_background')
    design = np.stack([np.ones_like(X), X - 25, Y - 25], axis=-1)
    keep = background.copy()
    for _ in range(3):
        beta = np.linalg.lstsq(design[keep] / err[keep, None], sci[keep] / err[keep], rcond=None)[0]
        residual = sci - design @ beta
        keep = background & (abs(residual) < 3 * np.maximum(err, 1))
        if keep.sum() < 100:
            return reject('structured_background')
    b = sci - design @ beta
    neighbor, labels, ncat = segment_mask(b, err, dqmask, cx, cy, neighbors, row)
    mask = neighbor | dqmask
    audit.update(MASKRAW=mask.astype('u1'), NEIGHBOR=neighbor.astype('u1'), LABELS=labels.astype('i4'), BGSUB=b)
    row.update(background_plane=[float(v) for v in beta], catalog_neighbors=ncat, raw_masked_pixels=int(mask.sum()), dq_bad_pixels=int(bad.sum()))
    wx, wy = registration(dx, dy)
    outmask = (abs(wy) > 1e-08).astype(int) @ mask.astype(int) @ (abs(wx.T) > 1e-08).astype(int) > 0
    centered = wy @ np.where(mask, 0, b) @ wx.T
    variance = wy ** 2 @ np.where(mask, 0, err ** 2) @ (wx ** 2).T
    centered[outmask] = np.nan
    variance[outmask] = np.nan
    audit.update(REGISTER=centered, OUTMASK=outmask.astype('u1'), VAR=variance)
    row.update(masked_pixels=int(outmask.sum()), valid_fraction=float((~outmask).mean()))
    if outmask[20, 20] or np.mean(outmask[RR <= 3]) > 0.3 or np.mean(outmask[RR <= 6]) > 0.5:
        return reject('insufficient_clean_core')
    if (~outmask).mean() < 0.45:
        return reject('too_little_unmasked_area')
    ap = (RR <= 6) & ~outmask
    measured_flux = float(np.sum(centered[ap]))
    ap_weights = wy.T @ ap.astype(float) @ wx
    flux = float(row['psf_flux'])
    snr = measured_flux / np.sqrt(np.sum(np.where(mask, 0, err ** 2) * ap_weights ** 2))
    row['aperture_snr'] = float(snr)
    if flux <= 0 or snr < 100:
        return reject('low_aperture_snr')
    processed = centered / flux
    variance /= flux ** 2
    audit['PROC'] = processed
    audit['VAR'] = variance
    row['negative_flux_fraction'] = float(-np.nansum(np.minimum(processed, 0)))
    if row['negative_flux_fraction'] > 0.1:
        return reject('background_structure_or_noise')
    row.update(accepted=True, reason='passed_v2_screen')
    return (row, audit)

def unique_clean(rows):
    rows = [r for r in rows if r['accepted']]
    rows.sort(key=lambda r: r.get('rank_score', r.get('target_distance_px', 0) / 1000) + 2 * (1 - r['valid_fraction']) + 20 / r['aperture_snr'])
    labels = sky_groups(rows)
    seen = set()
    out = []
    for r, g in zip(rows, labels):
        if int(g) not in seen:
            out.append(dict(r, sky_group=int(g)))
            seen.add(int(g))
    return out

def stack(cube):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        clipped = np.array(cube, copy=True)
        for _ in range(2):
            reference = np.nanmedian(clipped, axis=0)
            for a in clipped:
                good = (RR <= 6) & np.isfinite(a) & np.isfinite(reference)
                scale = np.sum(a[good] * reference[good]) / np.sum(reference[good] ** 2)
                if np.isfinite(scale) and scale > 0:
                    a /= scale
        med = np.nanmedian(clipped, axis=0)
        mad = 1.4826 * np.nanmedian(abs(clipped - med), axis=0)
        bad = abs(clipped - med) > np.maximum(5 * mad, 0.0005)
        clipped[bad] = np.nan
        m = np.nanmedian(clipped, axis=0)
        n = np.sum(np.isfinite(clipped), axis=0)
        scat = 1.4826 * np.nanmedian(abs(clipped - m), axis=0)
    return (m, n, scat)

class EmpiricalBuilder:
    def __init__(self, root, targets):
        self.root = Path(root)
        self.out = self.root / 'refined_v2'
        self.targets = targets
        self.lock = threading.Lock()
        self.locks = {}
        for name in ('queries', 'cutouts', 'neighbors', 'processed', 'records', 'psfs', 'diagnostics'):
            (self.out / name).mkdir(parents=True, exist_ok=True)

    def query(self, path, filters, columns=COLS, maxpages=10):
        if path.exists():
            return json.loads(path.read_text())
        req = {'service': 'Mast.Catalogs.Filtered.Wfc3Psf.Uvis', 'params': {'columns': columns, 'filters': filters}, 'format': 'json', 'pagesize': 5000, 'page': 1}
        rows = []
        total = 0
        for page in range(1, maxpages + 1):
            req['page'] = page
            for attempt in range(3):
                try:
                    response = requests.post('https://mast.stsci.edu/api/v0/invoke', data={'request': json.dumps(req)}, timeout=90)
                    response.raise_for_status()
                    p = response.json()
                    if p.get('status') != 'COMPLETE':
                        raise RuntimeError(str(p)[:300])
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(1)
            rows.extend(p['data'])
            total = p.get('paging', {}).get('rowsFiltered', len(rows))
            if page >= p.get('paging', {}).get('pagesFiltered', 1):
                break
        result = {'request': req, 'rows': rows, 'total': total, 'truncated': len(rows) < total}
        save(path, result)
        return result

    def candidates(self, band, tier):
        filters = [val('filter', band), val('chip', 2), val('subarray', 0), val('n_sat_pixels', 0), val('mtflag', 0), rng('qfit', 1e-06, 0.1), rng('iso_index', 10, 100000), rng('pixc', 2000, 60000), rng('sky', 0, 500), rng('psf_flux', 5000, 1000000000.0), rng('exptime', 10, 30000), rng('focus', -10, 10), rng('psf_x_center', 26, 1000 if tier == 'near' else 4070), rng('psf_y_center', 26, 1000 if tier == 'near' else 2022)]
        q = self.query(self.out / 'queries' / f'{band}_{tier}.json', filters)
        q['unique_rows'] = len({r['id'] for r in q['rows']})
        q['duplicate_rows'] = len(q['rows']) - q['unique_rows']
        print(band, tier, 'catalog', len(q['rows']), 'total', q['total'], 'truncated', q['truncated'], flush=True)
        return q

    def neighbor_catalog(self, row):
        name = row['rootname'].lower()
        with self.lock:
            lock = self.locks.setdefault('neighbor:' + name, threading.Lock())
        with lock:
            path = self.out / 'neighbors' / f'{name}.json'
            q = self.query(path, [val('rootname', name), val('chip', 2)], 'id,x_cal,y_cal,pixc,qfit,psf_flux', maxpages=30)
            complete = not q['truncated'] and len({v['id'] for v in q['rows']}) == q['total']
            mode = 'complete_exposure_catalog'
            if not complete:
                path = self.out / 'neighbors' / f"{name}_{row['id']}_window.json"
                q = self.query(path, [val('rootname', name), val('chip', 2), rng('x_cal', row['x_cal'] - 38, row['x_cal'] + 38), rng('y_cal', row['y_cal'] - 38, row['y_cal'] + 38)], 'id,x_cal,y_cal,pixc,qfit,psf_flux', maxpages=1)
                mode = 'complete_local_window'
                if q['truncated'] or len({v['id'] for v in q['rows']}) != q['total']:
                    raise ValueError('incomplete_neighbor_window')
            if row['id'] not in {v['id'] for v in q['rows']}:
                raise ValueError('target_missing_from_neighbor_catalog')
        return (q['rows'], {'file': str(path.relative_to(self.out)), 'mode': mode, 'rows': len(q['rows'])})

    def plane(self, row, ext):
        name = f"{row['rootname'].lower()}_{row['id']}_{row['filter']}_ext{ext}.fits"
        path = self.out / 'cutouts' / name
        old = self.root / 'cutouts' / f"{row['rootname'].lower()}_{row['id']}_{row['filter']}_flc.fits"
        if ext == 1 and old.exists():
            path = old
        uri = f"mast:WFC3PSF/url/cgi-bin/fitscut.cgi?red={row['rootname'].lower()}_flc[{ext}]&size=51&x={row['x_cal']}&y={row['y_cal']}&badvalue=-999999&format=fits/{name}"
        if path.exists():
            payload = path.read_bytes()
        else:
            for attempt in range(3):
                try:
                    response = requests.get('https://mast.stsci.edu/api/v0.1/Download/file', params={'uri': uri}, timeout=60)
                    response.raise_for_status()
                    payload = response.content
                    with fits.open(io.BytesIO(payload)) as h:
                        h.verify('exception')
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(1)
            path.write_bytes(payload)
        with fits.open(io.BytesIO(payload)) as h:
            h.verify('exception')
            a = next((v.data.astype(float) for v in h if v.data is not None and v.data.ndim == 2))
        if a.shape != (51, 51):
            raise ValueError('wrong_cutout_shape')
        if ext == 3 and (not np.isfinite(a).all() or not np.allclose(a, np.rint(a))):
            raise ValueError('invalid_DQ_values')
        return (a, dict(path=str(path.relative_to(self.root)), uri=uri, sha256=hashlib.sha256(payload).hexdigest()))

    def rank(self, rows, band):
        positions = np.array([[e['detector_x'], e['detector_y']] for t in self.targets if t['filter'] == band for e in t['exposures']])
        out = []
        seen_ids = set()
        for r in rows:
            source_key = (r['rootname'].lower(), r['id'])
            if source_key in seen_ids:
                continue
            seen_ids.add(source_key)
            if r.get('psf_ra') is None or r.get('psf_dec') is None:
                continue
            d = float(np.min(np.linalg.norm(positions - [r['psf_x_center'], r['psf_y_center']], axis=1)))
            score = d / 1000 + r['qfit'] / 0.1 + np.sqrt(max(r['sky'], 0) + 16) / max(r['pixc'], 1) * 100
            out.append(dict(r, target_distance_px=d, rank_score=score))
        out.sort(key=lambda r: r['rank_score'])
        labels = sky_groups(out)
        ng = Counter()
        ne = Counter()
        selected = []
        for r, g in zip(out, labels):
            if ng[int(g)] >= 2 or ne[r['rootname']] >= 8:
                continue
            r['candidate_sky_group'] = int(g)
            selected.append(r)
            ng[int(g)] += 1
            ne[r['rootname']] += 1
        return selected

    def process(self, row):
        key = f"record:{row['rootname'].lower()}:{row['id']}"
        with self.lock:
            lock = self.locks.setdefault(key, threading.Lock())
        with lock:
            return self.process_locked(row)

    def process_locked(self, row):
        key = f"{row['rootname'].lower()}_{row['id']}_{row['filter']}"
        record = self.out / 'records' / f'{key}.json'
        stamp = self.out / 'processed' / f'{key}.fits'
        if record.exists():
            cached = json.loads(record.read_text())
            if cached.get('algorithm_version') == VERSION and cached.get('reason') != 'processing_error':
                return dict(cached, **{k: row[k] for k in ('rank_score', 'target_distance_px') if k in row})
        row = dict(row)
        try:
            sci, sp = self.plane(row, 1)
            err, ep = self.plane(row, 2)
            dq, dp = self.plane(row, 3)
            neighbors, neighbor_provenance = ([], {'mode': 'skipped_flagged_central_pixels'}) if np.any((dq != 0)[R <= 1.5]) else self.neighbor_catalog(row)
            result, audit = process_arrays(row, sci, err, dq, neighbors)
            result['inputs'] = {'SCI': sp, 'ERR': ep, 'DQ': dp}
            result['neighbor_catalog'] = neighbor_provenance
            header = fits.Header()
            header['FILTER'] = row['filter']
            header['SOURCEID'] = row['id']
            header['PASS'] = result['accepted']
            header['METHOD'] = VERSION
            hdus = [fits.PrimaryHDU(header=header)] + [fits.ImageHDU(a.astype('f4') if a.dtype.kind == 'f' else a, name=name) for name, a in audit.items()]
            fits.HDUList(hdus).writeto(stamp, overwrite=True, checksum=True)
            result['processed_file'] = str(stamp.relative_to(self.out))
        except Exception as e:
            result = dict(row, algorithm_version=VERSION, accepted=False, reason='processing_error', error=str(e)[:400])
        save(record, result)
        return result

    def read_processed(self, row, ext='PROC'):
        with fits.open(self.out / row['processed_file']) as h:
            return h[ext].data.astype(float)

    def finish_band(self, band, rows, tiers):
        clean = unique_clean(rows)
        info = dict(filter=band, attempted=len(rows), accepted_detections=sum((r['accepted'] for r in rows)), distinct_sky_groups=len(clean), status='insufficient_clean_stars', tiers=tiers, rejection_reasons=dict(Counter((r['reason'] for r in rows if not r['accepted']))), focus_matched=False, drizzle_matched=False)
        save(self.out / 'diagnostics' / f'{band}_selection.json', rows)
        if len(clean) < 10:
            return info
        roots = sorted({r['rootname'] for r in clean})
        held = set(roots[::5])
        train = [r for r in clean if r['rootname'] not in held]
        test = [r for r in clean if r['rootname'] in held]
        if len(train) < 8 or len(test) < 2:
            train = clean
            test = []
        cube = np.array([self.read_processed(r) for r in train])
        m, n, scatter = stack(cube)
        minimum = max(3, int(np.ceil(0.25 * len(train))))
        selected = None
        for size in [41, 35, 31, 25, 21]:
            b = (41 - size) // 2
            sl = np.s_[b:41 - b, b:41 - b]
            if np.isfinite(m[sl]).all() and n[sl].min() >= minimum:
                selected = (size, b, sl)
                break
        if selected is None:
            info['status'] = 'insufficient_wing_coverage'
            return info
        size, b, sl = selected
        total = float(m[sl].sum())
        psf = m[sl] / total
        if total <= 0:
            raise ValueError('invalid_stack_normalization')
        res = []
        metrics = []
        y, x = np.indices(psf.shape)
        rad = np.hypot(x - (size - 1) / 2, y - (size - 1) / 2)
        core = rad <= 6
        wing = (rad > 6) & (rad <= 10)
        old = fits.getdata(self.root / 'psfs' / f'MAST_UVIS2_{band}_trial.fits').astype(float)[sl] if (self.root / 'psfs' / f'MAST_UVIS2_{band}_trial.fits').exists() else None
        legacy_path = self.root / 'diagnostics' / f'{band}_selection.json'
        legacy_rows = [v for v in json.loads(legacy_path.read_text()) if v.get('stack_role') == 'training'] if legacy_path.exists() else []
        legacy_coords = SkyCoord([v['psf_ra'] for v in legacy_rows] * u.deg, [v['psf_dec'] for v in legacy_rows] * u.deg) if legacy_rows else None
        for r in test:
            a = self.read_processed(r)[sl]
            ok = core & np.isfinite(a)
            amp = float(np.sum(a[ok] * psf[ok]) / np.sum(psf[ok] ** 2))
            residual = a - amp * psf
            metric = {'id': r['id'], 'rootname': r['rootname'], 'core_l1_v2': float(np.sum(abs(residual[ok])) / np.sum(abs(a[ok])))}
            w = wing & np.isfinite(a)
            metric['wing_l1_v2'] = float(np.sum(abs(residual[w])) / np.sum(abs(a[ok]))) if w.any() else None
            overlap = legacy_coords is not None and np.min(SkyCoord(r['psf_ra'] * u.deg, r['psf_dec'] * u.deg).separation(legacy_coords).arcsec) < 0.3
            metric['overlaps_v1_training_sky_group'] = bool(overlap)
            if old is not None and (not overlap):
                amp1 = float(np.sum(a[ok] * old[ok]) / np.sum(old[ok] ** 2))
                d = a - amp1 * old
                metric['core_l1_v1'] = float(np.sum(abs(d[ok])) / np.sum(abs(a[ok])))
                metric['wing_l1_v1'] = float(np.sum(abs(d[w])) / np.sum(abs(a[ok]))) if w.any() else None
            metrics.append(metric)
            res.append(residual)
        header = fits.Header()
        header['FILTER'] = band
        header['INSTRUME'] = 'WFC3'
        header['DETECTOR'] = 'UVIS'
        header['CCDCHIP'] = 2
        header['NSTAR'] = len(train)
        header['NTEST'] = len(test)
        header['PSFSTAT'] = 'TRIAL_V2'
        header['DQMASK'] = True
        header['FOCMATCH'] = False
        header['DRIZMAT'] = False
        header['OVERSAMP'] = 1
        header['BUNIT'] = 'FRACTION'
        header['METHOD'] = 'BILINEAR'
        header['PIXSCALE'] = 0.03962
        header['HISTORY'] = 'One detection per 0.3 arcsec sky group; exposure-disjoint holdout. Native-frame trial.'
        header['HISTORY'] = 'Segmentation + catalog-neighbor + DQ masks before compact interpolation; not target matched.'
        hdus = [fits.PrimaryHDU(psf.astype('f4'), header), fits.ImageHDU(n[sl].astype('i2'), name='NCONTRIB'), fits.ImageHDU((scatter[sl] / total).astype('f4'), name='SCATTER')]
        if res:
            hdus.append(fits.ImageHDU(np.nanmedian(res, axis=0).astype('f4'), name='HOLDRES'))
        fits.HDUList(hdus).writeto(self.out / 'psfs' / f'MAST_UVIS2_{band}_v2_trial.fits', overwrite=True, checksum=True)
        for r in train:
            r['stack_role'] = 'training'
        for r in test:
            r['stack_role'] = 'heldout_exposure'
        save(self.out / 'diagnostics' / f'{band}_stack_members.json', train + test)
        save(self.out / 'diagnostics' / f'{band}_holdout.json', metrics)
        convergence = []
        for count in sorted(set([min(len(train), v) for v in [10, 20, 40, 80, len(train)]])):
            sm, _, _ = stack(cube[:count])
            sm = sm[sl]
            finite = np.isfinite(sm)
            sm /= np.nansum(sm)
            convergence.append({'n': count, 'core_l1_to_full': float(np.nansum(abs(sm[core] - psf[core])) / psf[core].sum()), 'finite_fraction': float(finite.mean())})
        save(self.out / 'diagnostics' / f'{band}_convergence.json', convergence)
        info.update(status='trial_stack' if len(train) >= 30 else 'limited_sample_trial', training=len(train), heldout=len(test), stamp_size=size, min_contributors=int(n[sl].min()), exposure_roots=len(roots), median_core_l1_v2=float(np.median([v['core_l1_v2'] for v in metrics])) if metrics else None, median_core_l1_v1=float(np.median([v['core_l1_v1'] for v in metrics if 'core_l1_v1' in v])) if any(('core_l1_v1' in v for v in metrics)) else None, median_wing_l1_v2=float(np.median([v['wing_l1_v2'] for v in metrics if v['wing_l1_v2'] is not None])) if metrics else None, max_target_distance_px=max((r.get('target_distance_px', 0) for r in clean)))
        info['paired_holdout_count'] = sum(('core_l1_v1' in v for v in metrics))
        info['paired_median_core_l1_v2'] = float(np.median([v['core_l1_v2'] for v in metrics if 'core_l1_v1' in v])) if info['paired_holdout_count'] else None
        info['paired_median_wing_l1_v1'] = float(np.median([v['wing_l1_v1'] for v in metrics if v.get('wing_l1_v1') is not None])) if info['paired_holdout_count'] else None
        info['paired_median_wing_l1_v2'] = float(np.median([v['wing_l1_v2'] for v in metrics if v.get('wing_l1_v1') is not None])) if info['paired_holdout_count'] else None
        return info

    def run(self, band, max_candidates=800, goal=100):
        rows = []
        seen = set()
        tiers = []
        for tier in ['near', 'chipwide']:
            q = self.candidates(band, tier)
            tiers.append({'tier': tier, 'catalog_rows': len(q['rows']), 'total': q['total'], 'truncated': q['truncated'], 'unique_catalog_rows': q['unique_rows'], 'duplicate_catalog_rows': q['duplicate_rows']})
            poolrows = [r for r in self.rank(q['rows'], band) if r['id'] not in seen]
            quota = min(max_candidates - len(rows), 400 if tier == 'near' else max_candidates)
            for start in range(0, min(len(poolrows), quota), 40):
                batch = poolrows[start:start + min(40, max_candidates - len(rows))]
                with ThreadPoolExecutor(max_workers=8) as pool:
                    for r in pool.map(self.process, batch):
                        rows.append(r)
                        seen.add(r['id'])
                clean = unique_clean(rows)
                save(self.out / 'diagnostics' / f'{band}_selection.json', rows)
                print(band, tier, 'attempted', len(rows), 'accepted', sum((r['accepted'] for r in rows)), 'unique', len(clean), 'rejections', dict(Counter((r['reason'] for r in rows if not r['accepted']))), flush=True)
                if len(clean) >= goal or len(rows) >= max_candidates:
                    break
            if len(unique_clean(rows)) >= goal or len(rows) >= max_candidates:
                break
        result = self.finish_band(band, rows, tiers)
        print('RESULT', json.dumps(result), flush=True)
        return result
