"""Inspect archived HST solutions and validate exposure registration on sources."""
from astropy.io import fits
from .matching import detect_sources,validate


def inspect_wcs(inputs):
    from stwcs.wcsutil import altwcs
    records=[]
    for p in inputs:
        with fits.open(p) as hd:
            chips=[]
            for i,h in enumerate(hd):
                if h.name!='SCI':continue
                names=altwcs.wcsnames(hd,ext=i)
                chips.append(dict(extver=h.header.get('EXTVER',1),solutions=names,primary=h.header.get('WCSNAME'),
                    gaia_solutions={k:v for k,v in names.items() if 'GAIA' in str(v).upper()},fit_metadata={k:h.header.get(k,hd[0].header.get(k)) for k in ('WCSTYPE','RMS_RA','RMS_DEC','NMATCHES','FIT_RMS','FITGEOM','CRDER1','CRDER2','CATALOG','IDCTAB','D2IMFILE','NPOLFILE')}))
            records.append(dict(exposure=p.name,chips=chips,headerlets=[dict(hdrname=h.header.get('HDRNAME'),wcsname=h.header.get('WCSNAME')) for h in hd if h.name=='HDRLET']))
    return records


def validate_relative(inputs,options,paths):
    catalogs=[detect_sources(p,options['detection']) for p in inputs];results=[]
    for p,cat in zip(inputs,catalogs):cat.write(paths.file('intermediate','alignment',p.stem+'_relative_sources','ecsv'),format='ascii.ecsv',overwrite=True)
    for p,cat in zip(inputs[1:],catalogs[1:]):
        report,matched=validate(cat,catalogs[0],options['validation'])
        report['exposure']=p.name;report['reference']=inputs[0].name;results.append(report)
        matched.write(paths.file('intermediate','alignment',p.stem+'_relative_matches','ecsv'),format='ascii.ecsv',overwrite=True)
    visits=[]
    for p in inputs:
        header=fits.getheader(p);visits.append((header.get('PROPOSID'),str(header.get('LINENUM','')).split('.')[0]))
    valid=all(r['status']=='validated' for r in results) and bool(results)
    return dict(status='single_exposure' if len(inputs)==1 else ('validated' if valid else 'unresolved_relative_registration'),
                pairs=results,multiple_visits=len(set(visits))>1,
                may_combine=valid or len(inputs)==1,
                warnings=[] if valid else ['Relative registration not fully constrained; absolute image fit cannot repair exposure smearing or misregistration'])
