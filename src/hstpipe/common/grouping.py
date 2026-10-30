from __future__ import annotations
def groups(exposures,instrument):
    if instrument=='stis':
        from ..stis.grouping import group_exposures
    else:
        from ..uvis.grouping import group_exposures
    return group_exposures(exposures)
