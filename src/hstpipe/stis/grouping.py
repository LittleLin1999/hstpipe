from collections import defaultdict

def group_exposures(exposures):
    result=defaultdict(list)
    for e in exposures:
        m=e.metadata
        if e.product!='flt':continue
        if m['DETECTOR']!='FUV-MAMA' or m['OPT_ELEM']!='G140L':raise ValueError('STIS engine currently supports FUV-MAMA G140L')
        element=m['OPT_ELEM'];extra=(m['CENWAVE'],m['APERTURE'],m['PROPAPER'],m['BINAXIS1'],m['BINAXIS2'])
        result[(m['TARGNAME'],element,extra)].append(e)
    return dict(result)
