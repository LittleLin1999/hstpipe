"""DrizzlePac astrometry on working copies of calibrated exposures."""
def align(inputs,options=None):
    from drizzlepac import tweakreg
    options=dict(options or {})
    if 'updatehdr' in options and not options['updatehdr']:raise ValueError('Alignment must update working-copy WCS')
    options.update(updatehdr=True,interactive=False)
    tweakreg.TweakReg([str(p) for p in inputs],**options)
