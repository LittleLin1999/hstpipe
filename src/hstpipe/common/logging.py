from __future__ import annotations
import logging

def reduction_logger(paths):
    logger=logging.getLogger('hstpipe.'+str(paths.directory));logger.setLevel(logging.INFO)
    if not logger.handlers:
        file=logging.FileHandler(paths.area('logs',create=True)/'reduction.log')
        console=logging.StreamHandler()
        for handler in (file,console):
            handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'));logger.addHandler(handler)
    return logger
