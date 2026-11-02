"""WFC3/UVIS bit definitions and purpose-specific masks (WFC3 DHB v6, §3.2)."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from .fits import image_data

VERSION='uvis-dq-2026-01'
BITS={1:'SOFTERR',2:'DATALOST',4:'DETECTORPROB',8:'DATAMASKED_UNUSED',16:'HOTPIX',32:'CTETAIL_UNSTABLE',64:'WARMPIX',128:'BADBIAS',
      256:'FULL_WELL_SATURATION',512:'BADFLAT',1024:'TRAP_SINK',2048:'ATOD_SATURATION',4096:'ASTRODRIZZLE_CR',8192:'CALWF3_REJECT',16384:'CROSSTALK_UNUSED',32768:'RESERVED'}
SATURATION=256|2048
DRIZZLE_CR=4096
CALWF3_CR=8192
SCIENCE_GOOD_BITS=16|64

def native_compatible_dq(hd):
    """DrizzlePac 3.11 uses ``dq & ~bits``; NumPy 2 rejects that on uint16.

    Signed int32 represents all 16 UVIS bits exactly. Only caller-owned working
    HDUs are promoted; public resetbits=4096 and native algorithms stay intact.
    """
    changed=[]
    if int(np.__version__.split('.')[0])<2:return changed
    for h in hd:
        if h.name!='DQ':continue
        data=image_data(h)
        if data.dtype.kind=='u' and data.dtype.itemsize<=2:
            h.data=np.asarray(data,dtype='i4');changed.append(int(h.header.get('EXTVER',1)))
    if changed:
        hd[0].header['HPDQCAST']=('INT32','Working DQ dtype for native resetbits / NumPy 2 compatibility')
    return changed

def bit_statistics(dq):
    a=np.asarray(dq,dtype='u2')
    return {str(bit):dict(name=name,pixels=int(np.count_nonzero(a&bit)),fraction=float(np.mean((a&bit)!=0))) for bit,name in BITS.items() if np.any(a&bit)}

def ignored_bits(options,purpose='centroid',fresh_cr=False):
    fallback=4096 if purpose in ('centroid','detection') else 0
    bits=int(options.get(purpose+'_ignore_bits',options.get('dq_ignore_bits',fallback)))
    if purpose in ('centroid','detection'):bits&=65535^SATURATION
    bits&=65535^CALWF3_CR
    if fresh_cr:bits&=65535^DRIZZLE_CR
    if not 0<=bits<=65535:raise ValueError('Invalid UVIS DQ mask bits')
    return bits

def masks(hd,version,options,purpose='centroid'):
    from scipy.ndimage import binary_dilation
    a=np.asarray(hd['SCI',version].data);dq=np.asarray(image_data(hd['DQ',version]),dtype='u2') if ('DQ',version) in hd else np.zeros(a.shape,'u2')
    fresh=bool(hd[0].header.get('HPCRRUN',False));ignore=ignored_bits(options,purpose,fresh)
    invalid=~np.isfinite(a)
    if ('ERR',version) in hd:
        error=image_data(hd['ERR',version]);invalid|=~np.isfinite(error)|(error<=0)
    coverage=np.zeros(a.shape,bool)
    if 'WHT' in hd:coverage=~np.isfinite(hd['WHT'].data)|(hd['WHT'].data<=0)
    base=((dq&np.uint16(65535^ignore))!=0)|invalid|coverage
    # Grow selected CR footprints only, never every nonzero DQ pixel.
    grow_bits=int(options.get('grow_dq_bits',CALWF3_CR))|(DRIZZLE_CR if fresh else 0)
    radius=int(options.get('dq_grow_radius',1));grown=base.copy()
    if radius>0 and np.any(dq&grow_bits):grown|=binary_dilation((dq&grow_bits)!=0,iterations=radius)
    edge=int(options.get('edge_pixels',2))
    if edge>0:grown[:edge]=True;grown[-edge:]=True;grown[:,:edge]=True;grown[:,-edge:]=True
    return dq,base,grown

def policy(options=None,nimages=None):
    options=options or {}
    selected=options.get('astrometry',{}).get('relative',{}).get('detection',{})
    return dict(version=VERSION,astrometry_inherited_4096='Reassess using measured morphology/centroids and HST repeatability; never edit raw DQ',
        centroid_ignore_bits=ignored_bits(selected,'centroid'),detection_ignore_bits=ignored_bits(selected,'detection'),centroid_saturation_rejected=SATURATION,calwf3_8192='Always retained and excluded',
        background_ignore_bits=options.get('background',{}).get('dq_ignore_bits',0),
        driz_sep_bits=0 if nimages==1 else options.get('drizzle',{}).get('driz_sep_bits',SCIENCE_GOOD_BITS),final_bits=0 if nimages==1 else options.get('drizzle',{}).get('final_bits',SCIENCE_GOOD_BITS),
        science_reason='MDRIZTAB 336 minus full-well saturation 256; stable/warm flags 16+64 accepted only with multi-exposure rejection; stricter single-exposure/diagnostic mask',
        resetbits='4096 on working copies for a new CR pass; 0 when carrying newly generated masks into split final combination',grow_bits=8192,blanket_dilation=False)

def mdriztab_recommendation(image,nimages):
    from .amplifiers import reference_path
    h=fits.getheader(image);value=h.get('MDRIZTAB')
    if not value:return dict(status='not_available')
    try:
        path=reference_path(value)
        with fits.open(path) as hd:
            rows=[r for r in hd[1].data if str(r['filter']).strip() in ('ANY',h['FILTER']) and int(r['numimages'])<=nimages]
            selected=max(rows,key=lambda r:(str(r['filter']).strip()==h['FILTER'],int(r['numimages'])))
            values={k:(selected[k].item() if hasattr(selected[k],'item') else selected[k]) for k in ('filter','numimages','driz_sep_bits','final_bits','driz_sep_pixfrac','final_pixfrac')}
        return dict(status='inspected',reference=str(path),recommendation=values,applied_automatically=False)
    except (OSError,ValueError,KeyError) as exc:return dict(status='unavailable',reference=value,error=str(exc))
