"""Conservative off-source diagnostics; AstroDrizzle owns sky subtraction.

Masks must be supplied and inspected by the caller. No object classification or
amplifier geometry is inferred from the science signal. Label maps use zero for
excluded pixels and positive integers for independently selected sky regions or
amplifiers. Statistics remain in each SCI extension's declared input units.
"""
import numpy as np
from astropy.io import fits

def assess_regions(exposure,mask):
    report=[]
    with fits.open(exposure) as science,fits.open(mask) as masks:
        for h in science:
            if h.name!='SCI':continue
            v=h.header.get('EXTVER',1);labels=masks['SCI',v].data
            if labels.shape!=h.data.shape:raise ValueError('Sky/amplifier mask shape differs from SCI')
            if not np.isfinite(labels).all() or np.any(labels<0) or not np.equal(labels,np.floor(labels)).all():
                raise ValueError('Masks must contain nonnegative integer region labels')
            dq=science['DQ',v].data
            for label in np.unique(labels[labels>0]):
                good=(labels==label)&(dq==0)&np.isfinite(h.data)
                values=h.data[good].astype(float)
                if len(values)<30:raise ValueError('Insufficient independent sky pixels for region assessment')
                report.append(dict(extver=v,region=int(label),npixel=len(values),
                    mean=float(values.mean()),rms=float(values.std(ddof=1)),
                    standard_error=float(values.std(ddof=1)/np.sqrt(len(values))),
                    unit=h.header.get('BUNIT','UNKNOWN'),
                    caveat='Conditional on supplied source-free masks; correlated noise not modeled'))
    return report
