"""Optional UVIS amplifier equalization using detector coordinates and PFLT.

Follows STScI's UVIS amplifier equalization notebook (Sunnquist 2019):
source-masked amp median minus mean of amp medians, multiplied by PFLT.
No inference from a subarray's midpoint is allowed. Disabled by default.
"""
import os
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from .background import estimate

def detector_x(header,shape):
    if header.get('BINAXIS1',1)!=1 or header.get('BINAXIS2',1)!=1:raise ValueError('Amp equalization requires unbinned geometry')
    return (np.arange(shape[1])+1-header.get('LTV1',0))/header.get('LTM1_1',1)

def reference_path(value):
    if '$' in value:
        prefix,name=value.split('$',1)
        if prefix not in os.environ:raise ValueError(f'Missing reference environment variable {prefix}')
        return Path(os.environ[prefix])/name
    return Path(value).expanduser().resolve()

def equalize(hd,options):
    corrections={};regions=[];report=dict(applied=False,regions=[],warnings=[])
    for h in hd:
        if h.name!='SCI':continue
        v=h.header.get('EXTVER',1);chip=h.header.get('CCDCHIP')
        if chip not in (1,2):raise ValueError('Missing CCDCHIP for amp mapping')
        x=detector_x(h.header,h.shape);labels=np.where(x<=2048,0,1)
        amps='AB' if chip==1 else 'CD'
        active=set(str(hd[0].header.get('CCDAMP','')))&set(amps)
        if not active:raise ValueError('CCDAMP does not identify a readout amplifier for this chip')
        if len(active)==1:labels=np.full(x.shape,amps.index(next(iter(active))),dtype=int)
        for side in np.unique(labels):
            columns=labels==side;indices=np.flatnonzero(columns)
            section=np.s_[:,indices[0]:indices[-1]+1]
            _,mask,sky=estimate(h.data[section],hd['ERR',v].data[section],hd['DQ',v].data[section],
                WCS(h.header,hd).slice(section),[hd[0].header['RA_TARG'],hd[0].header['DEC_TARG']],dict(options,mode='none'))
            good=np.zeros(h.shape,dtype=bool);good[section]=~mask;n=int(good.sum())
            reliable=n>=options['min_pixels'] and n/max(1,(labels==side).sum()*h.shape[0])>=options['min_blank']
            level=float(np.median(h.data[good])) if n else 0.0
            error=max(sky['uncertainty'],float(np.std(h.data[good]))/np.sqrt(max(1,n))) if n else 1e30
            regions.append((v,amps[side],good,labels==side,level,error,reliable))
    if len(regions)<2:
        report['warnings'].append('Only one amplifier sampled; no relative offset correction');return corrections,report
    mean=float(np.mean([r[4] for r in regions]));significant=any(abs(r[4]-mean)>options['significance']*r[5] for r in regions)
    report['regions']=[dict(extver=v,amp=a,median=level,offset=level-mean,uncertainty=error,reliable=bool(ok)) for v,a,g,c,level,error,ok in regions]
    if not significant or not all(r[6] for r in regions):
        report['warnings'].append('Amp differences unsupported by reliable blank-sky measurements; no correction');return corrections,report
    if not options.get('amplifier_correction'):
        report['warnings'].append('Assessment only; correction disabled');return corrections,report
    flatpath=reference_path(hd[0].header['PFLTFILE'])
    from ..common.provenance import identity
    report['flat']=identity(flatpath)
    with fits.open(flatpath) as flat:
        for v,a,g,columns,level,error,ok in regions:
            h=hd['SCI',v];chip=h.header['CCDCHIP']
            matches=[f for f in flat if f.name=='SCI' and f.header.get('CCDCHIP')==chip]
            if len(matches)!=1:raise ValueError('Flat SCI chip identity ambiguous')
            f=matches[0];dx=detector_x(h.header,h.shape)
            dy=(np.arange(h.shape[0])+1-h.header.get('LTV2',0))/h.header.get('LTM2_2',1)
            fx=dx*f.header.get('LTM1_1',1)+f.header.get('LTV1',0)-1
            fy=dy*f.header.get('LTM2_2',1)+f.header.get('LTV2',0)-1
            if not np.allclose(fx,np.round(fx)) or not np.allclose(fy,np.round(fy)):raise ValueError('Nonintegral flat mapping')
            ix=np.round(fx).astype(int);iy=np.round(fy).astype(int)
            if min(ix)<0 or max(ix)>=f.shape[1] or min(iy)<0 or max(iy)>=f.shape[0]:raise ValueError('Flat does not cover science geometry')
            patch=f.data[np.ix_(iy,ix)]
            if not np.isfinite(patch).all() or np.any(patch<=0):raise ValueError('Invalid reference flat')
            correction=corrections.setdefault(v,np.zeros(h.shape,dtype='f4'))
            correction[:,columns]=(level-mean)*patch[:,columns]
    # Apply only after all geometry/reference checks have passed.
    for v,c in corrections.items():hd['SCI',v].data-=c
    report['applied']=True
    return corrections,report
