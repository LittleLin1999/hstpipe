"""Established imaging drizzle and cosmic-ray rejection, with native DQ/ERR."""
def drizzle(inputs,output,options=None):
    from drizzlepac import astrodrizzle
    options=dict(options or {})
    if any(k in options for k in ('input','output')):raise ValueError('Use explicit input/output arguments')
    options.setdefault('driz_cr',len(inputs)>1)
    options.setdefault('final_wht_type','IVM')
    options.setdefault('context',True)
    options.setdefault('clean',False)
    astrodrizzle.AstroDrizzle([str(p) for p in inputs],output=str(output),**options)
