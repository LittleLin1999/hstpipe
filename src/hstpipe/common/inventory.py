from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from astropy.io import fits
from .provenance import identity
from .fits import imsets

@dataclass
class Exposure:
    path: Path
    metadata: dict
    product: str
    identity: dict

FIELDS=('INSTRUME','DETECTOR','TARGNAME','OPT_ELEM','FILTER','APERTURE','PROPAPER','CENWAVE',
        'OBSTYPE','SUBARRAY','CCDAMP','BINAXIS1','BINAXIS2','POSTARG1','POSTARG2','PA_APER',
        'DARKCORR','FLATCORR','PCTECORR','EXPSTART','EXPTIME','ROOTNAME')

def discover(directory,instrument=None):
    """Flat local inventory; duplicates are identified by content, not filename."""
    if instrument=='uvis':
        from ..uvis.inventory import discover as discover_uvis
        return discover_uvis(directory)
    rows=[];seen={}
    for path in sorted(Path(directory).glob('*.fits')):
        with fits.open(path) as hd:
            h=hd[0].header;m={k:h.get(k) for k in FIELDS}
            if 'SCI' not in hd:continue
            if instrument and str(m['INSTRUME']).upper()!=('STIS' if instrument=='stis' else 'WFC3'):continue
            m['imsets']=imsets(path)
            m['geometry']=[{k:a.header.get(k) for k in ('EXTVER','CCDCHIP','NAXIS1','NAXIS2','LTV1','LTV2','LTM1_1','LTM2_2')} for a in hd if a.name=='SCI']
        ident=identity(path)
        if ident['sha256'] in seen:raise ValueError(f'Duplicate FITS contents: {path}, {seen[ident["sha256"]]}')
        seen[ident['sha256']]=path
        # Product identity is checked against detector metadata by grouping.
        product=path.stem.rsplit('_',1)[-1].lower()
        rows.append(Exposure(path.resolve(),m,product,ident))
    return rows
