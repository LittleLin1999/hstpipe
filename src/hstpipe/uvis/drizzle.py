"""Official AstroDrizzle owns distortion, median/blot CR rejection and drizzle."""
def drizzle(inputs,output,options=None):
    from drizzlepac import astrodrizzle
    opts=dict(options or {})
    if any(k in opts for k in ('input','output')):raise ValueError('Use explicit input/output arguments')
    if opts.get('group') not in (None,''):raise ValueError('UVIS AstroDrizzle must include all SCI extensions; restrictive group selection is not supported')
    opts.pop('group',None)
    opts.setdefault('driz_cr',len(inputs)>1)
    opts.setdefault('final_wht_type','EXP')
    from astropy.io import fits
    from pathlib import Path
    existing=[p for p in inputs if Path(p).is_file()]
    if any('raw' in Path(p).resolve().parts for p in existing):raise ValueError('AstroDrizzle may modify DQ; pass FLC working copies, never raw inputs')
    from .dq import native_compatible_dq
    for p in existing:
        with fits.open(p,mode='update') as hd:native_compatible_dq(hd)
    current_masks=bool(existing) and all(fits.getheader(p).get('HPCRRUN',False) for p in existing)
    reset=4096 if opts['driz_cr'] or not current_masks else 0
    opts.setdefault('driz_sep_bits',80 if opts['driz_cr'] else 0)
    opts.setdefault('final_bits',80 if opts['driz_cr'] else 0)
    if len(inputs)==1:opts['driz_sep_bits']=0;opts['final_bits']=0
    opts.update(context=True,build=True,resetbits=reset,updatewcs=False,preserve=False,
                restore=False,overwrite=False,clean=False,skysub=False,skyuser='HPSKY',
                final_units='cps',proc_unit='electrons',driz_cr_corr=False,num_cores=1)
    if any(opts.get(k) is not None for k in ('final_scale','final_rot','final_ra','final_dec','final_outnx','final_outny')):opts['final_wcs']=True
    if opts.get('driz_sep_scale') is not None:opts['driz_sep_wcs']=True
    from pathlib import Path
    def local(path):
        p=Path(path)
        return p.name if p.resolve().parent==Path.cwd().resolve() else str(p)
    astrodrizzle.AstroDrizzle([local(p) for p in inputs],output=local(output),**opts)

def annotate_product(path,weight,ninput,alignment,dq_bits=None):
    from astropy.io import fits
    with fits.open(path,mode='update') as hd:
        from astropy.nddata.bitmask import interpret_bit_flags
        hd[0].header['HPDQSCI']=(int(interpret_bit_flags(80 if dq_bits is None else dq_bits)) if ninput>1 else 0,'Effective final DQ bits treated as good')
        hd[0].header['HPNINPUT']=ninput;hd[0].header['HPALIGN']=alignment
        hd[0].header['HPQASTAT']='REVIEW';hd[0].header['HPCRCLN']=(ninput>1,'Multi-image CR rejection run; inspect HPQASTAT')
        hd[0].header.add_history('hstpipe UVIS: calibrated CTE-corrected FLC; originals preserved')
        hd[0].header.add_history('No additional PAM correction; no ERR inferred from WHT')
        hd['WHT'].header['WHTTYPE']=weight
        hd['WHT'].header['BUNIT']='SECONDS' if weight=='EXP' else 'WEIGHT'
        hd['WHT'].header.add_history('Drizzle effective weight; correlated noise not represented')
        if weight=='EXP':hd['WHT'].header.add_history('Exposure-time weights are NOT inverse variance')
        else:hd['WHT'].header.add_history('Native noise weights propagated by AstroDrizzle; no output covariance')
        for h in hd:
            if hasattr(h,'add_checksum'):h.add_checksum()


def apply_cr_masks(inputs):
    """Interpret official good=1 CR masks and OR bit 4096 on working DQ.

    driz_cr_corr=False intentionally leaves detector SCI and DQ unchanged;
    split-stage execution must explicitly carry masks into final drizzle.
    """
    from pathlib import Path
    import numpy as np
    from astropy.io import fits
    for path in inputs:
        path=Path(path)
        with fits.open(path,mode='update') as hd:
            root=path.stem.rsplit('_',1)[0]
            for h in hd:
                if h.name!='SCI':continue
                v=h.header.get('EXTVER',1)
                mask=path.parent/f'{root}_sci{v}_crmask.fits'
                if not mask.exists():raise FileNotFoundError(f'Missing official CR mask: {mask}')
                good=fits.getdata(mask)
                if good.shape!=h.shape or not np.isin(good,[0,1]).all():raise ValueError('Invalid official CR mask')
                hd['DQ',v].data[good==0]|=4096
            hd[0].header['HPCRRUN']=True
