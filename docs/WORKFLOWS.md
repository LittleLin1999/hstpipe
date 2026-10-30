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

UVIS starts from calibrated FLC working copies, with optional TweakReg alignment
and official AstroDrizzle cosmic-ray rejection, sky handling and drizzle. The
public align/drizzle/inspect_flc/validate_drizzle functions are independently
usable. backgrounds are handled by AstroDrizzle options (skysub/skymethod),
not a custom detector correction. Optional background_masks/amplifier_masks
mappings supply integer-label SCI mask HDUs for exposure-specific off-source
assessments. No automatic amplifier-offset correction is applied. These masks
must be scientifically verified by the caller. Output SCI/WHT/CTX product:

    redux/{target}/uvis/{filter}/[setup]/products/{target}_{filter}_drc.fits

Retained work lies in intermediate/alignment and intermediate/drizzle;
assessment, alignment and drizzle QA lie in diagnostics. No real UVIS reduction
was requested or validated in this STIS run. UVIS currently has whole-run cache
reuse; its underlying DrizzlePac steps retain their own intermediates.
