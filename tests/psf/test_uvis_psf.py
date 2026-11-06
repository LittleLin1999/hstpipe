"""UVIS PSF CLI/configuration and independent scientific regression checks."""
import json
from pathlib import Path
import sys
import types
import numpy as np
import pytest
import yaml
from hstpipe.uvis.psf.config import read_config, read_targets, read_records


def config(tmp_path, **overrides):
    values = dict(workdir='psf', filters=['F275W'])
    values.update(overrides)
    path = tmp_path / 'psf.yaml'
    path.write_text(yaml.safe_dump(values))
    return path


def test_config_relative_to_yaml_and_no_import_io(tmp_path, monkeypatch):
    path = config(tmp_path, targets='inputs/targets.json')
    monkeypatch.chdir(tmp_path.parent)
    c = read_config(path)
    assert c['targets'] == tmp_path / 'inputs/targets.json'
    assert c['manifest'] == tmp_path / 'psf/drizzled/input_products/manifest.json'
    assert not c['workdir'].exists()


@pytest.mark.parametrize('bad', [dict(chip=1), dict(filters='F275W'),
    dict(filters=['../../x']), dict(products=['bogus']), dict(models=[]),
    dict(goal=20, limit=10), dict(goal=True), dict(typo=True)])
def test_invalid_configuration_fails(tmp_path, bad):
    with pytest.raises(ValueError):
        read_config(config(tmp_path, **bad))


@pytest.mark.parametrize('instrument', ['uvis', 'stis'])
def test_instrument_specific_cli_dispatch(tmp_path, monkeypatch, instrument):
    from hstpipe.cli import main
    called = []
    module = types.ModuleType(f'hstpipe.{instrument}.psf.cli')
    module.run = lambda *args: called.append(args)
    monkeypatch.setitem(sys.modules, module.__name__, module)
    stage = 'empirical' if instrument == 'uvis' else 'library'
    main([instrument, 'psf', '--config', 'example.yaml', '--stage', stage])
    assert called == [('example.yaml', stage)]


def test_targets_missing_requested_filter(tmp_path):
    c = read_config(config(tmp_path))
    c['workdir'].mkdir()
    c['targets'].write_text(json.dumps([dict(filter='F336W', exposures=[dict(detector_x=1, detector_y=2)])]))
    with pytest.raises(ValueError, match='No detector positions'):
        read_targets(c)


def test_missing_science_selection_fails_before_work(tmp_path):
    c = read_config(config(tmp_path))
    c['manifest'].parent.mkdir(parents=True)
    c['manifest'].write_text(json.dumps(dict(records=[])))
    with pytest.raises(ValueError, match='No science products'):
        read_records(c)


def test_failure_receipt_and_cwd_preserved(tmp_path, monkeypatch):
    from hstpipe.uvis.psf import cli
    def fail(c):
        raise RuntimeError('insufficient stars')
    monkeypatch.setattr(cli, 'empirical', fail)
    original = Path.cwd()
    with pytest.raises(RuntimeError, match='insufficient stars'):
        cli.run(config(tmp_path), 'empirical')
    assert Path.cwd() == original
    receipt = json.loads(next((tmp_path / 'psf/runs').glob('*.json')).read_text())
    assert receipt['status'] == 'failed'
    assert 'insufficient stars' in receipt['error']


def test_insufficient_stars_propagates_and_keeps_summary(tmp_path, monkeypatch):
    from hstpipe.uvis.psf import cli
    c = read_config(config(tmp_path))
    c['workdir'].mkdir()
    c['targets'].write_text(json.dumps([dict(filter='F275W', exposures=[dict(detector_x=1, detector_y=2)])]))
    class Builder:
        def __init__(self, root, targets):
            self.out = root / 'refined_v2'
        def run(self, *args):
            return dict(filter='F275W', status='insufficient_clean_stars')
    module = types.ModuleType('hstpipe.uvis.psf.empirical')
    module.EmpiricalBuilder = Builder
    monkeypatch.setitem(sys.modules, module.__name__, module)
    with pytest.raises(RuntimeError, match='insufficient_clean_stars'):
        cli.empirical(c)
    summary = json.loads((c['workdir'] / 'refined_v2/summary.json').read_text())
    assert summary[0]['status'] == 'insufficient_clean_stars'


def test_workdir_lock_excludes_second_writer(tmp_path):
    from hstpipe.uvis.psf.cli import lock
    with lock(tmp_path):
        with pytest.raises(RuntimeError, match='Another PSF run'):
            with lock(tmp_path):
                pytest.fail('Acquired a second lock')


def test_render_cropped_psf_preserves_center_and_phase():
    pytest.importorskip('scipy')
    from hstpipe.uvis.psf.render import render_empirical
    y, x = np.indices((61, 61))
    for size in [21, 25, 31, 35, 41]:
        a, b = np.indices((size, size))
        model = np.exp(-((a-size//2)**2 + (b-size//2)**2) / 2)
        for shift in [0, .25]:
            result = render_empirical(model, x-30-shift, y-30)
            assert np.unravel_index(np.argmax(result), result.shape) == (30, 30)
            assert abs((result*x).sum()/result.sum() - (30+shift)) < .01


def empirical_module():
    pytest.importorskip('scipy')
    pytest.importorskip('photutils')
    pytest.importorskip('requests')
    from hstpipe.uvis.psf import empirical
    return empirical


def test_central_dq_rejected_before_segmentation():
    e = empirical_module()
    sci = np.ones((51, 51)); dq = np.zeros_like(sci); dq[25, 25] = 4096
    result, _ = e.process_arrays({}, sci, sci, dq, [])
    assert not result['accepted']
    assert result['reason'] == 'flagged_central_pixels'


def test_fresh_stack_without_v1_and_disjoint_holdouts(tmp_path):
    e = empirical_module()
    from astropy.io import fits
    builder = e.EmpiricalBuilder(tmp_path, [])
    y, x = np.indices((41, 41))
    psf = np.exp(-((x-20)**2 + (y-20)**2) / 2); psf /= psf.sum()
    rows = []
    for i in range(15):
        path = builder.out / 'processed' / f'{i}.fits'
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(psf, name='PROC')]).writeto(path)
        rows.append(dict(id=i, rootname=f'exposure{i:02}', accepted=True,
                         psf_ra=10+i*.01, psf_dec=2., valid_fraction=1.,
                         aperture_snr=200., target_distance_px=0.,
                         processed_file=str(path.relative_to(builder.out))))
    result = builder.finish_band('F275W', rows, [])
    assert result['training'] == 12 and result['heldout'] == 3
    output = fits.getdata(builder.out / 'psfs/MAST_UVIS2_F275W_v2_trial.fits')
    np.testing.assert_allclose(output, psf, atol=1e-7)
    members = json.loads((builder.out / 'diagnostics/F275W_stack_members.json').read_text())
    train = {r['rootname'] for r in members if r['stack_role'] == 'training'}
    held = {r['rootname'] for r in members if r['stack_role'] == 'heldout_exposure'}
    assert train.isdisjoint(held)


def test_masked_neighbor_cannot_leak_into_registered_pixels():
    e = empirical_module()
    import photutils
    if int(photutils.__version__.split('.')[0]) < 3:
        pytest.skip('Production segmentation API requires photutils >=3')
    y, x = np.indices((51, 51))
    sci = 30 + 20000*np.exp(-((x-25.3)**2 + (y-24.8)**2)/2)
    err = np.full_like(sci, 5.)
    dq = np.zeros_like(sci); dq[25, 35] = 1
    row = dict(id=1, psf_flux=20000*2*np.pi, x_cal=100., y_cal=100.)
    a, aa = e.process_arrays(row, sci, err, dq, [])
    sci[25, 35] = 1e12
    b, bb = e.process_arrays(row, sci, err, dq, [])
    assert a['accepted'] and b['accepted']
    np.testing.assert_allclose(aa['PROC'], bb['PROC'], equal_nan=True)


def test_validation_rejects_changed_product(tmp_path):
    from astropy.io import fits
    from hstpipe.uvis.psf.validation import validate, sha
    c = read_config(config(tmp_path))
    c['manifest'].parent.mkdir(parents=True)
    reference = c['manifest'].parent / 'science.fits'
    reference.write_bytes(b'science identity')
    model = tmp_path / 'model.fits'; model.write_bytes(b'model identity')
    base = 'example_F275W_hstpipe'
    folder = c['workdir'] / 'drizzled/products/example'; folder.mkdir(parents=True)
    path = folder / (base + '_psf_empirical_v2.fits')
    p = np.zeros((61, 61)); p[30, 30] = 1
    fits.HDUList([fits.PrimaryHDU(p), fits.ImageHDU(p, name='WHT'), fits.ImageHDU(p.astype('i4'), name='CTX')]).writeto(path, checksum=True)
    record = dict(target='example', filter='F275W', products={
        'hstpipe': dict(path=reference.name, cutout=dict(sha256=sha(reference)))})
    c['manifest'].write_text(json.dumps(dict(records=[record])))
    provenance = dict(validation=dict(context_equal_fraction=1., weight_relative_core_max=0., weight_relative_median=0.),
                      outputs=[dict(model='empirical_v2', path=str(path.relative_to(c['workdir']/'drizzled')),
                                    sha256=sha(path), model_source=str(model), model_sha256=sha(model))])
    (folder / (base + '_provenance.json')).write_text(json.dumps(provenance))
    assert validate(c)[0]['status'] == 'passed'
    # Existing science bundles have no skipped_models field in their provenance.
    c['models'] = ['empirical_v2', 'anderson']
    assert validate(c)[1]['status'] == 'skipped_no_library'
    c['anderson_library'].mkdir(parents=True)
    library = c['anderson_library'] / 'STDPSF_WFC3UV_F275W.fits'
    library.write_bytes(b'available model')
    with pytest.raises(ValueError, match='Missing anderson'):
        validate(c)
    c['models'] = ['empirical_v2']
    with path.open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(ValueError, match='changed'):
        validate(c)


def test_legacy_review_sample_mask_regressions():
    e = empirical_module()
    import photutils
    if int(photutils.__version__.split('.')[0]) < 3:
        pytest.skip('Production segmentation API requires photutils >=3')
    rng = np.random.default_rng(4)
    y, x = np.indices((51, 51))
    star = 30000*np.exp(-((x-25.25)**2+(y-24.7)**2)/(2*1.1**2))
    a = 30+star+rng.normal(0, 5, (51, 51)); err = np.sqrt(46+star)
    dq = np.zeros((51, 51))
    row = dict(id=1, x_cal=1000., y_cal=500., psf_flux=float(star.sum()))
    result, data = e.process_arrays(row, a, err, dq, [])
    assert result['accepted']
    assert np.unravel_index(np.nanargmax(data['PROC']), data['PROC'].shape) == (20, 20)
    b = a+300*np.exp(-((x-44)**2+(y-26)**2)/(2*1.3**2))
    neighbors = [dict(id=2, x_cal=1018.75, y_cal=501.3, pixc=300., psf_flux=3185.)]
    _, data = e.process_arrays(row, b, np.sqrt(err**2+np.maximum(b-a, 0)), dq, neighbors)
    assert data['MASKRAW'][26, 44] and data['OUTMASK'][21, 39]
    wide = 200*np.exp(-((x-44)**2+(y-26)**2)/(2*3.**2))
    _, data = e.process_arrays(row, a+wide, np.sqrt(err**2+wide), dq, [])
    assert data['MASKRAW'][26, 44] and data['MASKRAW'][26, 39]


def test_pilot_rejections_are_distinct_from_processing_errors(tmp_path, monkeypatch):
    from hstpipe.uvis.psf import cli
    c = read_config(config(tmp_path))
    c['workdir'].mkdir()
    c['review_sample'].write_text(json.dumps([dict(filter='F275W', id=1)]))
    class Builder:
        def __init__(self, root, targets):
            self.out = root / 'refined_v2'
        def process(self, row):
            return dict(row, accepted=False, reason='flagged_central_pixels')
    module = types.ModuleType('hstpipe.uvis.psf.empirical')
    module.EmpiricalBuilder = Builder
    monkeypatch.setitem(sys.modules, module.__name__, module)
    assert cli.pilot(c)[0]['reason'] == 'flagged_central_pixels'
    Builder.process = lambda self, row: dict(row, reason='processing_error')
    with pytest.raises(RuntimeError, match='processing failed'):
        cli.pilot(c)
