"""Header-driven FLC inventory. Raw files are opened read-only."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from ..common.inventory import Exposure
from ..common.provenance import identity
from ..common.fits import imsets

FIELDS = ('INSTRUME DETECTOR TARGNAME RA_TARG DEC_TARG ROOTNAME PROPOSID LINENUM POSTARG1 POSTARG2 PA_APER '
          'FILTER EXPTIME EXPSTART EXPEND DATE-OBS TIME-OBS APERTURE PROPAPER CCDAMP '
          'SUBARRAY OBSTYPE BIASCORR BLEVCORR DARKCORR FLSHCORR FLASHDUR FLASHLVL FLATCORR '
          'PCTECORR PHOTCORR DQICORR FGSLOCK GYROMODE EXPFLAG QUALITY '
          'FLUXCORR CCDGAIN READNSEA READNSEB READNSEC READNSED IDCTAB NPOLFILE D2IMFILE PFLTFILE BIASFILE DARKFILE PCTETAB CAL_VER').split()

def validate_flc(path):
    versions=imsets(path)
    if len(versions)!=len(set(versions)):raise ValueError('Duplicate SCI EXTVER')
    with fits.open(path) as hd:
        p=hd[0].header;m={k:p.get(k) for k in FIELDS}; warnings=[]
        if (m['INSTRUME'],m['DETECTOR'],m['OBSTYPE'])!=('WFC3','UVIS','IMAGING'):
            raise ValueError(f'Not WFC3 UVIS imaging: {path}')
        if not m['FILTER'] or not m['EXPTIME'] or m['EXPTIME']<=0:
            raise ValueError(f'Missing science identity/exposure time: {path}')
        for k in ('BIASCORR','BLEVCORR','DARKCORR','FLATCORR','PCTECORR','PHOTCORR','DQICORR'):
            if m[k]!='COMPLETE': warnings.append(f'{k}={m[k]} (expected COMPLETE)')
        if m['FLUXCORR']!='COMPLETE':warnings.append(f'FLUXCORR={m["FLUXCORR"]}; verify chip photometric normalization; no additional correction applied')
        if (m['FLASHDUR'] or 0)>0 and m['FLSHCORR']!='COMPLETE':warnings.append('Post-flash correction incomplete')
        if m['FGSLOCK'] not in ('FINE','FINE/GYRO'):warnings.append(f'Guiding requires review: FGSLOCK={m["FGSLOCK"]}')
        if m['EXPFLAG'] not in (None,'NORMAL'):warnings.append(f'EXPFLAG={m["EXPFLAG"]}')
        if m['QUALITY'] not in (None,'','OK','NORMAL'):warnings.append(f'QUALITY={m["QUALITY"]}')
        m['visit']=str(m['LINENUM']).split('.')[0] if m['LINENUM'] is not None else None
        chips=[]
        for v in versions:
            for name in ('SCI','ERR','DQ'):
                if sum(a.name==name and a.header.get('EXTVER',1)==v for a in hd)!=1:raise ValueError(f'Ambiguous {name},{v} extension')
            h=hd['SCI',v]; w=WCS(h.header,hd);shape=h.shape
            if len(shape)!=2 or min(shape)<2:raise ValueError('SCI must be a two-dimensional image')
            if not w.has_celestial:raise ValueError('Missing celestial WCS')
            corners=w.all_pix2world([[0,0],[shape[1]-1,0],[shape[1]-1,shape[0]-1],[0,shape[0]-1]],0)
            if not np.isfinite(corners).all():raise ValueError('Nonfinite WCS footprint')
            dq=hd['DQ',v].data
            if dq is not None and not np.issubdtype(dq.dtype,np.integer):raise ValueError('DQ must contain integer bit flags')
            unit=h.header.get('BUNIT','').upper()
            if unit not in ('ELECTRONS','ELECTRON'):raise ValueError(f'FLC input units must be electrons, found {unit}')
            if not w.has_distortion:warnings.append(f'SCI,{v}: no embedded distortion model; review before drizzle')
            for k in ('PHOTFLAM','PHOTPLAM'):
                if k not in h.header and k not in p:warnings.append(f'SCI,{v}: missing {k}')
            chips.append(dict(extver=v,sci=hd.index_of(('SCI',v)),err=hd.index_of(('ERR',v)),dq=hd.index_of(('DQ',v)),
                shape=list(shape),has_distortion=bool(w.has_distortion),footprint=corners.tolist(),**{k:h.header.get(k) for k in
                ('CCDCHIP','LTV1','LTV2','LTM1_1','LTM2_2','BINAXIS1','BINAXIS2','WCSNAME','BUNIT','PHOTFLAM','PHOTPLAM','PHOTBW','MDRIZSKY','PHTFLAM1','PHTFLAM2','PHTRATIO','PHOTZPT')}))
        m.update(imsets=versions,geometry=chips,warnings=warnings)
        return m

def discover(directory):
    rows=[];seen=set()
    for path in sorted(Path(directory).glob('*_flc.fits')):
        m=validate_flc(path);ident=identity(path)
        if ident['sha256'] in seen:raise ValueError(f'Duplicate FLC content: {path}')
        seen.add(ident['sha256']);rows.append(Exposure(path.resolve(),m,'flc',ident))
    return rows

def scan(directory):
    """Isolate invalid inputs by target/filter without silently dropping them."""
    import re
    rows=[];errors=[];seen=set()
    for path in sorted(Path(directory).glob('*_flc.fits')):
        try:
            m=validate_flc(path);ident=identity(path)
            if ident['sha256'] in seen:raise ValueError('Duplicate FLC content')
            seen.add(ident['sha256']);rows.append(Exposure(path.resolve(),m,'flc',ident))
        except Exception as exc:
            try:h=fits.getheader(path)
            except Exception:h={}
            target=h.get('TARGNAME') or 'unidentified_'+re.sub('[^A-Za-z0-9_.+-]','_',path.stem)
            errors.append(dict(path=str(path.resolve()),target=target,filter=h.get('FILTER') or 'UNKNOWN',error=f'{type(exc).__name__}: {exc}'))
    return rows,errors
