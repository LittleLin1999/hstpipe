from __future__ import annotations
from pathlib import Path
from astropy.io import fits

def imsets(path):
    """Validate matched SCI/ERR/DQ shapes, supporting multiple chips/imsets."""
    with fits.open(path) as hd:
        versions=[int(h.header.get('EXTVER',1)) for h in hd if h.name=='SCI']
        if not versions:raise ValueError(f'No SCI HDU: {path}')
        for v in versions:
            shapes=[]
            for name in ('SCI','ERR','DQ'):
                hdu=hd[name,v]
                if hdu.data is not None:shapes.append(hdu.data.shape)
                elif 'PIXVALUE' in hdu.header:
                    shapes.append((hdu.header['NPIX2'],hdu.header['NPIX1']))
                else:raise ValueError(f'Empty {name} imset {v}: {path}')
            if len(set(shapes))!=1:raise ValueError(f'Mismatched imset {v}: {path}')
        return versions

def validate_product(path,required=('SCI',)):
    with fits.open(path,checksum=True) as hd:
        hd.verify('exception')
        for key in required:
            if key not in hd:raise ValueError(f'Missing {key}: {path}')
        for h in hd:
            if 'CHECKSUM' in h.header and h.verify_checksum()!=1:raise ValueError(f'Bad checksum: {path}')
    return Path(path)
