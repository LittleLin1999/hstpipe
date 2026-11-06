"""Report layout changes preserve scientific data and old bookmarks."""
import json
from pathlib import Path

from hstpipe.cli import main
from hstpipe.common.paths import ReductionPaths
from hstpipe.common.provenance import identity
from hstpipe.common.target_runner import save
from hstpipe.uvis.reporting import generate_report
from hstpipe.uvis.shallow import _destination


def test_new_report_has_one_entry_point(tmp_path):
    paths=ReductionPaths(tmp_path, 'T', 'uvis', 'F275W')
    save(paths.directory/'manifest.json', dict(target='T', filters={}))
    report=generate_report(paths.directory)
    assert report==paths.directory/'index.html'
    assert not paths.area('products').exists()
    assert not (tmp_path/'reports').exists()
    assert _destination(paths, Path('index.html'))==report


def test_refresh_preserves_science_and_old_bookmarks(tmp_path, capsys):
    paths=ReductionPaths(tmp_path, 'T', 'uvis', 'F275W')
    product=paths.product('drc')
    product.parent.mkdir(parents=True)
    product.write_bytes(b'unchanged science artifact')
    cache=paths.cache()
    cache.record('example', 'same-inputs', [product])
    receipt=cache.receipt('example')
    before=(identity(product), receipt.read_bytes())
    legacy=product.parent/'T_uvis_report.html'
    legacy.write_text('previous report')
    record=dict(filter='F275W', status='complete', products=[identity(product)],
                scientific_assessment=dict(status='science_ready'))
    save(paths.directory/'manifest.json', dict(target='T', filters={'F275W':record},
                                             report=dict(path=str(legacy))))
    for _ in range(2):
        main(['uvis', 'report', '--directory', str(paths.directory)])
        assert capsys.readouterr().out.strip()==str(paths.report())
        assert legacy.is_symlink() and legacy.readlink()==Path('../index.html')
        assert legacy.read_bytes()==paths.report().read_bytes()
        state=json.loads((paths.directory/'manifest.json').read_text())
        assert state['filters']=={'F275W':record}
        assert state['report']['path']==str(paths.report())
        assert (identity(product), receipt.read_bytes())==before
        assert cache.valid('example', 'same-inputs')


def test_report_command_rejects_non_target_directory(tmp_path):
    import pytest
    directory=tmp_path/'collected_reports'
    save(directory/'manifest.json', dict(target='T', filters={}))
    before=(directory/'manifest.json').read_bytes()
    with pytest.raises(SystemExit):
        main(['uvis', 'report', '--directory', str(directory)])
    assert (directory/'manifest.json').read_bytes()==before
    assert not (tmp_path/'T').exists()
