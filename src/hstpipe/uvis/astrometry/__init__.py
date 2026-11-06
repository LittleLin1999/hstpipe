"""Direct reference catalogs and validated two-pass UVIS astrometry."""
from .catalogs import query_gaia_dr3,query_legacy_dr10,query_footprint
__all__=['query_gaia_dr3','query_legacy_dr10','query_footprint']


def align_cross_filters(directory, selected=None, options=None):
    """Register validated provisional filter images before absolute calibration.

    Reads the target manifest and writes cross-filter match/validation records.
    Uses only HST measurements; preserves provisional and original FLC files.
    """
    from .cross_filter import align_cross_filters as align
    return align(directory, selected, options)


__all__.append('align_cross_filters')


def resume_cross_filters(directory, selected=None):
    """Reprocess saved provisional mosaics without rerunning FLC registration."""
    from .orchestration import resume_cross_filters as resume
    return resume(directory, selected)


__all__.append('resume_cross_filters')
