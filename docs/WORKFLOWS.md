# Workflows and output layout

Public API: hstpipe.stis.reduce_stis / hstpipe.uvis.reduce_uvis accept a YAML
filename, a ReductionConfig, or an equivalent dictionary. The shared CLI uses
hstpipe {stis,uvis} inspect/reduce. --dry-run never starts processing.

STIS science stages are exposed as the functions returned by
hstpipe.stis.pipeline.stages(): prepare, background, joint_drizzle, mast_coadd,
validate, extract. Their native Config and input data models are preserved;
these stages can be called from Python with their documented explicit inputs.
The main API records per-stage hash receipts and reloads validated products.
Input, parameters, package source, dependency version and native-executable
changes invalidate affected downstream stages. --resume reuses valid stages;
--force archives an old run and starts anew. Failed stages never get a receipt.

STIS outputs:

    redux/{target}/stis/{grating}/[setup]/
      products/{target}_{grating}_{s2d_drz,x2d_coadd,x1d}.fits
      intermediate/background/           # per-FLT models and rate files
      intermediate/flt/                  # input views, native operators/coadd arrays
      intermediate/x2d/                  # individually cleaned archive spectra
      intermediate/x1d/{flt,drz,x2d}/    # per-method extraction tables
      intermediate/registration.ecsv
      diagnostics/{background,coadd,extraction,fluxcal}/
      logs/{reduction,stages}.log
      config.yaml
      config.json                       # native STIS effective configuration
      manifest.json

The consolidated X1D contains FLT/DRZ/X2D BINTABLEs and separate boxcar/optimal
errors/masks. DQ aliases DQ_BOX; select optimal with DQ_OPT. No science defaults
or units were changed for this infrastructure migration.

UVIS uses calibrated CTE-corrected FLC working copies, discovers every SCI/ERR/DQ
imset, and validates internal exposure registration before provisional drizzle.
External Gaia/Legacy catalogs may guide source identities; relative fitting uses
measured HST centroids. Accepted absolute corrections on provisional images are
propagated with STScI TweakBack before a new final AstroDrizzle pass.

    redux/{target}/uvis/products/{target}_{filter}[_setup]_drc.fits

Intermediate directories are shallow: catalogs, background, alignment, drizzle,
and working. Diagnostics and logs are flat; filters/setup IDs are filenames.
A standalone HTML report under 1 MB is always produced, including on failure.

Persistent-source QA compares the final image to independent no-CR single-image
drizzles on the same WCS. Raw CR peaks and whole-region flagged-pixel fractions
are not treated as lost astronomical flux. Failed registration, astrometry or
source/image QA blocks publication and retains review products. Stage caches
allow unchanged upstream work to be reused. See [UVIS](uvis.md) for configuration,
validation criteria, FITS semantics and remaining science limitations.
