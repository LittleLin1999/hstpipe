import json
from pathlib import Path

import pytest

from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity
from hstpipe.common.target_runner import save
from hstpipe.common.uvis_cleanup import cleanup_uvis


def group(root, *, ready=True, name='F336W'):
    paths = ReductionPaths(root, 'TARGET', 'uvis', name)
    raw = root.parent / (name + '_raw_flc.fits')
    raw.parent.mkdir(parents=True, exist_ok=True)
    if not raw.exists():
        raw.write_bytes(b'original science data')
    product = paths.product('drc')
    product.parent.mkdir(parents=True, exist_ok=True)
    product.write_bytes(b'published science data')
    state = dict(schema_version=2, target='TARGET', filters={name: dict(
        filter=name, inputs=[identity(raw)], status='complete' if ready else 'relative_astrometry_unvalidated',
        products=[identity(product)] if ready else [],
        scientific_assessment=dict(status='science_ready' if ready else 'relative_astrometry_unvalidated'))})
    save(paths.directory / 'manifest.json', state)
    comparison = paths.file('intermediate', 'drizzle', 'joint_candidate_inspection', 'fits')
    comparison.parent.mkdir(parents=True, exist_ok=True)
    comparison.write_bytes(b'diagnostic coadd')
    for kind in ('provisional', 'unvalidated_final', 'source_qa_cutouts'):
        paths.file('intermediate', 'drizzle', kind, 'fits').write_bytes(kind.encode())
    paths.cache().record('reduction', 'unchanged-key', [product, comparison])
    return paths, raw, comparison


def test_preview_and_idempotent_cleanup_preserve_science_and_receipts(tmp_path):
    paths, raw, image = group(tmp_path / 'redux')
    protected = [raw, paths.product('drc'), paths.directory / 'manifest.json',
                 paths.cache().receipt('reduction')]
    protected += [p for p in image.parent.glob('*.fits') if p != image]
    before = {p: p.read_bytes() for p in protected}
    result = cleanup_uvis(tmp_path / 'redux')
    assert result['mode'] == 'preview' and result['bytes'] == image.stat().st_size
    assert image.exists() and not paths.area('logs').exists()
    result = cleanup_uvis(tmp_path / 'redux', apply=True)
    assert result['removed'] == [str(image)] and not image.exists()
    assert len(result['receipts']) == 1
    receipt = json.loads(Path(result['receipts'][0]).read_text())
    assert receipt['files'][0]['sha256'] and receipt['removed'] == [str(image)]
    assert {p: p.read_bytes() for p in protected} == before
    assert not paths.cache().valid('reduction', 'unchanged-key')
    assert cleanup_uvis(tmp_path / 'redux', apply=True)['removed'] == []


def test_blocked_group_is_preserved_even_with_old_product(tmp_path):
    paths, raw, image = group(tmp_path / 'redux', ready=False)
    result = cleanup_uvis(paths.directory, apply=True)
    assert not result['removed'] and image.exists()
    assert result['preserved'][0]['reason']


def test_changed_publication_fails_before_deletion(tmp_path):
    paths, _, image = group(tmp_path / 'redux')
    paths.product('drc').write_bytes(b'changed')
    with pytest.raises(ValueError, match='differs from its manifest'):
        cleanup_uvis(paths.directory, apply=True)
    assert image.exists()


def test_retry_requires_same_inputs_and_removes_only_diagnostic_alias(tmp_path):
    first, raw, image = group(tmp_path / 'run/redux', ready=False)
    retry, _, retry_image = group(tmp_path / 'run/retry/redux')
    state = json.loads((retry.directory / 'manifest.json').read_text())
    state['filters']['F336W']['inputs'] = [identity(raw)]
    save(retry.directory / 'manifest.json', state)
    alias = tmp_path / 'run/view/TARGET/uvis/intermediate/drizzle' / image.name
    alias.parent.mkdir(parents=True)
    alias.symlink_to(image)
    export = tmp_path / 'run/published.fits'
    export.symlink_to(retry.product('drc'))
    result = cleanup_uvis(tmp_path / 'run', apply=True)
    assert set(result['removed']) == {str(image), str(retry_image)}
    assert not alias.is_symlink() and export.exists()


def test_symlinked_comparison_and_mismatched_inputs_are_not_deleted(tmp_path):
    first, raw, image = group(tmp_path / 'run/redux', ready=False)
    retry, _, retry_image = group(tmp_path / 'run/retry/redux')
    outside = tmp_path / 'outside.fits'
    outside.write_bytes(b'keep outside data')
    retry_image.unlink()
    retry_image.symlink_to(outside)
    result = cleanup_uvis(tmp_path / 'run', apply=True)
    assert not result['removed'] and image.exists() and outside.exists()
    assert retry_image.is_symlink()


def test_cli_defaults_to_preview_and_explicit_apply(tmp_path, capsys):
    from hstpipe.cli import main
    paths, _, image = group(tmp_path / 'redux')
    main(['uvis', 'cleanup', '--output', str(tmp_path / 'redux')])
    assert json.loads(capsys.readouterr().out)['mode'] == 'preview' and image.exists()
    main(['uvis', 'cleanup', '--output', str(tmp_path / 'redux'), '--target', 'TARGET', '--apply'])
    assert json.loads(capsys.readouterr().out)['removed'] == [str(image)]


def test_running_group_lock_blocks_cleanup(tmp_path):
    from hstpipe.common.target_runner import target_lock
    paths, _, image = group(tmp_path / 'redux')
    with target_lock(paths.directory, paths.group_key):
        with pytest.raises(RuntimeError, match='Another reduction'):
            cleanup_uvis(paths.directory, apply=True)
    assert image.exists()


def test_symlinked_drizzle_directory_cannot_delete_external_data(tmp_path):
    paths, _, image = group(tmp_path / 'redux')
    directory = image.parent
    outside = tmp_path / 'outside'
    directory.rename(outside)
    directory.symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match='escapes reduction directory'):
        cleanup_uvis(paths.directory, apply=True)
    assert (outside / image.name).exists()
