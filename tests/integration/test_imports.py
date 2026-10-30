import subprocess,sys

def test_public_imports_without_drizzlepac():
    subprocess.run([sys.executable,'-c','import hstpipe; from hstpipe.stis import reduce_stis; from hstpipe.uvis import reduce_uvis; import sys; assert "drizzlepac" not in sys.modules'],check=True)

def test_cli_help():
    subprocess.run(['hstpipe','--help'],check=True,capture_output=True)

def test_only_hstpipe_console_entry_point():
    from importlib.metadata import distribution
    scripts = {entry.name for entry in distribution('hstpipe').entry_points
               if entry.group == 'console_scripts'}
    assert scripts == {'hstpipe'}
