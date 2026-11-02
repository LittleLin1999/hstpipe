"""Per-reduction logging; reopen handlers after an archived force/rebuild."""
import logging

def reduction_logger(paths):
    paths.area('logs',create=True)
    logfile=paths.file('logs',None,'reduction','log') if paths.instrument=='uvis' else paths.area('logs')/'reduction.log'
    logger=logging.getLogger('hstpipe.'+str(paths.directory));logger.setLevel(logging.INFO)
    for handler in list(logger.handlers):
        logger.removeHandler(handler);handler.close()
    for handler in (logging.FileHandler(logfile),logging.StreamHandler()):
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
        logger.addHandler(handler)
    return logger
