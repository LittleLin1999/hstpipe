from collections import defaultdict

def group_exposures(exposures):
    result=defaultdict(list)
    for e in exposures:
        m=e.metadata
        if e.product!='flc':continue
        if m['DETECTOR']!='UVIS' or m['PCTECORR']!='COMPLETE':raise ValueError('UVIS requires CTE-corrected FLC observations')
        element=m['FILTER'];extra=(m['APERTURE'],m['SUBARRAY'],m['CCDAMP'],str(m['geometry']))
        result[(m['TARGNAME'],element,extra)].append(e)
    return dict(result)
