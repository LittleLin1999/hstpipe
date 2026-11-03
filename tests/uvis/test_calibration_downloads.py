from pathlib import Path
from types import SimpleNamespace
import json

import pytest
from astropy.io import fits

from hstpipe.common.calibration import ensure_uvis_references
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity


def setup_case(tmp_path, monkeypatch, **options):
    monkeypatch.delenv('iref', raising=False)
    raw = tmp_path / 'raw/uvis'
    raw.mkdir(parents=True)
    image = raw / 'example_flc.fits'
    fits.PrimaryHDU(header=fits.Header({'IDCTAB':'iref$example_idc.fits',
                                      'MDRIZTAB':'N/A', 'BIASFILE':'iref$unused_bias.fits'})).writeto(image)
    config = SimpleNamespace(instrument='uvis', input=str(raw), output=str(tmp_path/'redux'),
                             options=dict(calibration=dict(timeout=1, retries=0, **options)))
    paths = ReductionPaths(tmp_path/'redux', 'TEST', 'uvis', 'F336W')
    return image, config, paths


def test_download_exact_header_reference_and_reuse(tmp_path, monkeypatch):
    import hstpipe.common.calibration as c
    image, config, paths = setup_case(tmp_path, monkeypatch)
    before = identity(image)
    calls = []
    def fetch(name, directory, timeout):
        calls.append(name)
        fits.PrimaryHDU().writeto(directory/name, checksum=True)
    monkeypatch.setattr(c, '_fetch', fetch)
    first = ensure_uvis_references([image, image], config, paths)
    config.options['calibration']['offline'] = True
    second = ensure_uvis_references([image], config, paths)
    assert calls == ['example_idc.fits']
    assert first['references'][0]['status'] == 'downloaded'
    assert second['references'][0]['status'] == 'reused'
    assert identity(image) == before
    assert first['references'][0]['identity']['sha256']


def test_offline_missing_required_reference_is_reported(tmp_path, monkeypatch):
    import hstpipe.common.calibration as c
    image, config, paths = setup_case(tmp_path, monkeypatch, offline=True)
    monkeypatch.setattr(c, '_fetch', lambda *a:pytest.fail('Network attempted offline'))
    with pytest.raises(RuntimeError, match='Required calibration reference unavailable'):
        ensure_uvis_references([image], config, paths)
    report = json.loads(paths.file('logs',None,'calibration_references').read_text())
    assert report['status'] == 'failed'
    assert report['references'][0]['status'] == 'unavailable'


def test_invalid_download_is_not_installed_and_retries_are_bounded(tmp_path, monkeypatch):
    import hstpipe.common.calibration as c
    image, config, paths = setup_case(tmp_path, monkeypatch)
    config.options['calibration']['retries'] = 1
    calls = []
    def fetch(name, directory, timeout):
        calls.append(name)
        (directory/name).write_text('not FITS')
    monkeypatch.setattr(c, '_fetch', fetch)
    with pytest.raises(RuntimeError, match='Required calibration reference unavailable'):
        ensure_uvis_references([image], config, paths)
    assert len(calls) == 2
    assert not (tmp_path/'calibration/uvis/example_idc.fits').exists()
    assert not list((tmp_path/'calibration/uvis').glob('.crds-*'))


def test_optional_mdriztab_failure_does_not_block(tmp_path, monkeypatch):
    image, config, paths = setup_case(tmp_path, monkeypatch, offline=True)
    with fits.open(image, mode='update') as hd:
        hd[0].header['IDCTAB'] = 'N/A'
        hd[0].header['MDRIZTAB'] = 'iref$missing_mdz.fits'
    report = ensure_uvis_references([image], config, paths)
    assert report['status'] == 'ready_with_warnings'


def test_download_directory_cannot_be_raw(tmp_path, monkeypatch):
    image, config, paths = setup_case(tmp_path, monkeypatch)
    config.options['calibration']['directory'] = config.input
    with pytest.raises(ValueError, match='raw input'):
        ensure_uvis_references([image], config, paths)


def test_configuration_paths_and_dry_run(tmp_path, monkeypatch):
    from hstpipe.common.config import ReductionConfig
    from hstpipe.uvis.config import settings
    config = tmp_path/'config.yaml'
    config.write_text('input: raw/uvis\noutput: redux\ninstrument: uvis\noptions:\n  calibration:\n    directory: calibration/uvis\n')
    c = ReductionConfig.read(config)
    assert c.options['calibration']['directory'] == str(tmp_path/'calibration/uvis')
    assert settings(c.options)['calibration']['auto_download'] is True


def test_stis_oref_references_are_downloaded_without_calibration(tmp_path, monkeypatch):
    import hstpipe.common.calibration as c
    image, config, _ = setup_case(tmp_path, monkeypatch)
    monkeypatch.delenv('oref', raising=False)
    with fits.open(image, mode='update') as hd:
        del hd[0].header['IDCTAB']
        hd[0].header['BIASFILE']='N/A'
        hd[0].header['PHOTTAB']='oref$example_pht.fits'
        hd[0].header['DARKFILE']='N/A'
    config.instrument='stis'
    config.options['reference_dir']=str(tmp_path/'calibration/stis')
    paths=ReductionPaths(tmp_path/'redux','TEST','stis','G140L')
    before=identity(image);calls=[]
    def fetch(name, directory, timeout):
        calls.append(name);fits.PrimaryHDU().writeto(directory/name)
    monkeypatch.setattr(c,'_fetch',fetch)
    report=c.ensure_references([image],config,paths)
    assert calls==['example_pht.fits']
    assert identity(image)==before and report['status']=='ready'


def test_stis_runner_prepares_references_before_engine_and_dry_run_is_offline(tmp_path, monkeypatch):
    import hstpipe.common.calibration as c
    import hstpipe.common.runner as runner
    from hstpipe.common.config import ReductionConfig
    image, _, paths=setup_case(tmp_path,monkeypatch)
    monkeypatch.delenv('oref',raising=False)
    with fits.open(image,mode='update') as hd:
        hd[0].header['IDCTAB']='N/A';hd[0].header['BIASFILE']='N/A'
        hd[0].header['PHOTTAB']='oref$example_pht.fits'
    row=SimpleNamespace(path=image)
    monkeypatch.setattr(runner,'discover',lambda *a:[row])
    monkeypatch.setattr(runner,'groups',lambda *a:{('TEST','G140L',None):[row]})
    fetched=[];engines=[]
    def fetch(name,directory,timeout):
        fetched.append(name);fits.PrimaryHDU().writeto(directory/name)
    monkeypatch.setattr(c,'_fetch',fetch)
    def engine(config,paths,rows):
        engines.append(config.options)
        assert 'calibration' not in config.options
        assert (Path(config.options['reference_dir'])/'example_pht.fits').exists()
        paths.area('products',create=True)
        product=paths.product('s2d_drz');fits.PrimaryHDU().writeto(product)
        return [product]
    config=ReductionConfig(str(image.parent),str(tmp_path/'redux'),'stis',options={'calibration':{'timeout':1,'retries':0}},resume=True)
    runner.run(config,engine,dry_run=True)
    assert not fetched and not engines and not Path(config.output).exists()
    before=identity(image)
    result=runner.run(config,engine)
    assert len(engines)==1 and fetched==['example_pht.fits'] and identity(image)==before
    manifest=json.loads((result[0]/'manifest.json').read_text())
    assert manifest['calibration_references']['status']=='ready'
    runner.run(config,engine)
    assert len(engines)==1 and fetched==['example_pht.fits']


def test_reference_cache_discovery_does_not_require_iref_before_download(tmp_path, monkeypatch):
    import hstpipe.common.calibration as c
    from hstpipe.uvis.config import settings
    from hstpipe.uvis.storage import references
    image, config, paths=setup_case(tmp_path,monkeypatch)
    with fits.open(image,mode='update') as hd:hd[0].header['PFLTFILE']='iref$example_pfl.fits'
    directory=tmp_path/'calibration/uvis'
    options=settings({'calibration':{'directory':str(directory)},'background':{'amplifier_correction':True}})
    rows=[SimpleNamespace(path=image,metadata={'PFLTFILE':'iref$example_pfl.fits'})]
    assert references(options,rows)==[]
    config.options=options
    monkeypatch.setattr(c,'_fetch',lambda name,path,timeout:fits.PrimaryHDU().writeto(path/name))
    ensure_uvis_references([image],config,paths)
    assert {p.name for p in references(options,rows)}=={'example_idc.fits','example_pfl.fits'}
