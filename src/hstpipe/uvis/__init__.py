"""UVIS imaging public API; importing does not import DrizzlePac."""
def reduce_uvis(config,*,dry_run=False):
    from .pipeline import reduce_uvis as run
    return run(config,dry_run=dry_run)
__all__=['reduce_uvis']
