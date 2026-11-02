"""Original detector-to-provisional WCS identity shared by internal stages."""
import json
from astropy.io import fits
from astropy.wcs import WCS

def signature(image):
    import hashlib
    with fits.open(image) as hd:
        h=hd['SCI'];key=h.header.get('HPORIGK',' ')
        header=WCS(h.header,hd,key=key).to_header(relax=True,key=' ')
        for k in list(header):
            # Fit uncertainty metadata do not change the coordinate transform.
            # STWCS restoration can inherit these from the corrected primary.
            if k.startswith(('WCSNAME','CRDER','CSYER')):del header[k]
        digest=hashlib.sha256(h.data.tobytes())
        canonical={k:(format(v,'.14g') if isinstance(v,float) else v) for k,v in header.items()}
        digest.update(json.dumps(canonical,sort_keys=True).encode())
        return digest.hexdigest()



def detector_wcs(hdulist,extver):
    """Use STWCS's complete HST model; retain lookup arrays via the open HDUList."""
    primary=hdulist[0].header
    if primary.get('INSTRUME')=='WFC3' and primary.get('DETECTOR')=='UVIS':
        from stwcs.wcsutil import HSTWCS
        return HSTWCS(hdulist,ext=('SCI',int(extver)))
    return WCS(hdulist['SCI',int(extver)].header,hdulist)
