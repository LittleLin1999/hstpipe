from __future__ import annotations
def groups(exposures,instrument,field=None):
    if instrument=='stis':
        from ..stis.grouping import group_exposures
    else:
        from ..uvis.grouping import group_exposures
    return group_exposures(exposures,field=field) if instrument=='uvis' else group_exposures(exposures)
