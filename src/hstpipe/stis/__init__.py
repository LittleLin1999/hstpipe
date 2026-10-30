"""STIS spectroscopy public API; heavy dependencies load at execution."""
def reduce_stis(config,*,dry_run=False):
    from .pipeline import reduce_stis as run
    return run(config,dry_run=dry_run)
__all__=['reduce_stis']
