"""CLI entry point used by `hstpipe stis psf`."""
from pathlib import Path
import os,logging
from .common import config,setup,save_json

STAGES=['discover','preview','download','recalibrate','screen','reduce','extract','library','validate','report','plot']

def run(config_path,stage='all'):
    path=Path(config_path).resolve();c=config(path);original=Path.cwd();os.chdir(c['workdir']);setup()
    try:
        from . import mast,extraction,psf,validation,reporting,calibration,hstpipe_adapter,product_profiles,library_plots
        operations=dict(discover=mast.discover,preview=mast.preview,download=mast.download,recalibrate=calibration.recalibrate,screen=extraction.extract,reduce=hstpipe_adapter.reduce_stars,extract=product_profiles.measure_products,library=psf.build_library,validate=validation.validate,report=reporting.report,plot=plot_products)
        selected=STAGES if stage=='all' else [stage]
        for name in selected:
            logging.info('PSF calibration stage %s',name)
            operations[name](c)
    finally:os.chdir(original)


def plot_products(c):
    from .library_plots import library_diagnostics
    from .position_plots import position_diagnostics
    library_diagnostics(c)
    return position_diagnostics(c)
