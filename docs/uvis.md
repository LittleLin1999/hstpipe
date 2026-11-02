# UVIS imaging and astrometric calibration

The default workflow uses calibrated, CTE-corrected individual FLC exposures:

1. Validate FLC calibration, SCI/ERR/DQ, geometry, photometry, and distortion WCS.
2. Preserve/inspect primary, alternate and headerlet WCS solutions; register and
   validate relative exposure alignment using the mandatory multi-strategy stage.
   Unresolved registration in any multi-exposure group blocks science combination.
3. Estimate source-masked backgrounds once on working copies. Use official
   AstroDrizzle separate drizzle, median, blot and CR masks to make a provisional
   distortion-corrected image.
4. Validate an existing Gaia WCS or use TweakReg to align provisional sources
   to directly queried Gaia DR3; explicitly fall back to Legacy Tractor when
   query/coverage/source/match/fit validation fails. A compatible validated optical
   anchor can transfer its reference frame across filters without a new query.
5. Archive the original provisional WCS with STWCS, validate the corrected WCS,
   and use `drizzlepac.tweakback.apply_tweak` on fresh relative-alignment FLC copies.
   Corrected per-exposure STWCS headerlets provide reproducible cache restoration.
6. Run a **new full AstroDrizzle pass on corrected individual FLC copies**, including
   CR rejection with the updated registration. Final pixels are not copied from
   the provisional image. Validate the final image against the reference catalog.

No calwf3 calibration or second PAM correction is performed. Raw files are never
opened for writing. Background, native pixel errors, detector DQ and flux units
retain their existing scientific definitions. STIS code and grating paths are unchanged.

## Output layout

```text
redux/{target}/
  stis/{grating}/                # unchanged
  uvis/
    products/{target}_{filter}[_setup]_drc.fits
    products/{target}_uvis_report.html
    intermediate/
      catalogs/                 # shared Gaia/Legacy cached cones and provenance
      background/               # rootname sky models and masks
      alignment/                # per-SCI catalogs, matches, WCS/headerlets
      drizzle/                  # provisional images, CR masks, source-QA cutouts
      working/                  # reusable relative-alignment FLC copies
    diagnostics/                # flat; target/filter/rootname-qualified filenames
    logs/                       # flat
    config.yaml
    manifest.json
```

There are no filter or mandatory setup directories. Group files and diagnostic
plots include target/filter/setup identifiers; exposure files use ROOTNAME.
Disconnected compatible-field groups have deterministic identifiers based on
sorted exposure ROOTNAME metadata. Unused stages are not created. The products
contain PRIMARY/SCI/WHT/CTX and AstroDrizzle HDRTAB when available.
SCI is electrons/s. WHTTYPE documents EXP/IVM/ERR weighting; EXP is not inverse
variance. No ERR is fabricated, and drizzle covariance is not claimed to be modeled.
`ASTSTAT`, `ASTCORR`, `ASTREF`, `ASTRELS`, `ASTMODE`, `ASTRMS`, `ASTNREF`, `ASTWCS` and `ASTEPCH`
make the final astrometric interpretation available in FITS. Unconstrained RMS
is left undefined rather than reported as zero. Image-quality review remains
separate from successful astrometric validation.

## Direct reference catalogs

`hstpipe.uvis.astrometry.catalogs.query_gaia_dr3(ra, dec, radius, observation_epoch)`
uses `astroquery.gaia.Gaia.launch_job_async` and **gaiadr3.gaia_source**. A numeric
observation epoch is a Julian year in TCB; Astropy Time or a date string is also
accepted. Gaia positions are propagated with `SkyCoord.apply_space_motion` from
Gaia's TCB reference epoch to actual HST dates. Reliable positive parallaxes
provide distances; otherwise angular propagation is used without nonphysical
negative distances. Missing radial velocities are assumed zero only where a
reliable distance is used, and that approximation is recorded. Missing proper
motions are excluded by default when an observation epoch is supplied. Per-axis
position/proper-motion covariance contributes to reported catalog uncertainty;
full perspective-error covariance is not supplied. Cached native Gaia positions
allow per-exposure epoch checks without requerying the service.

`query_legacy_dr10(ra, dec, radius)` uses **pyvo** at
`https://datalab.noirlab.edu/tap`, inspecting TAP_SCHEMA and querying
**ls_dr10.tractor** using Q3C_RADIAL_QUERY. Quality cuts check primary/unique
objects, coordinates, coverage, SNR, masks, blending and compact morphology.
Optional columns are used only when present. PHOTSYS or documented release codes
identify actual **DR10-South** and **DR9-North** rows; unknown codes remain explicitly
unresolved. Northern data are not described as newly processed DR10. Measured
motion can be applied to Gaia-associated Tractor stars with an explicit epoch;
other coadd centroid epochs remain uncertain. Legacy positions have resolution,
blending, epoch and catalog-systematic limitations and are not treated as Gaia DR3
precision. HST is never resampled onto Legacy's pixel grid.

Query cones cover actual WCS/coverage bounds plus a configurable margin, using
spherical coordinates across RA wrap and near poles. Deterministic ECSV catalogs
and JSON sidecars preserve exact queries, coordinates, radius, epoch, cuts,
endpoint, retrieval time and counts. Valid caches support offline operation.
Refresh is explicit. Bounded worker processes enforce wall-clock timeouts around
supported clients; transient errors have bounded retries. This is not a custom
HTTP client or another survey framework. No survey_query dependency is used. Empty results alone do not prove a lack of
coverage; returned objects with zero r/z observations are identified separately.
Query failure, no sources, insufficient references/matches and rejected fits are
separate attempt statuses in the alignment record. Fallback is never silent.

## Configuration and invocation

On Erebus:

```bash
cd /quasar/linxiaojing/lrd-squared/hst
source configs/environment.sh
hstpipe uvis inspect --input raw/uvis
hstpipe uvis reduce --input raw/uvis --output redux --config configs/uvis.yaml --dry-run
hstpipe uvis reduce --input raw/uvis --output redux --config configs/uvis.yaml --resume
# Add --target J1025+1402 --filter F467M for one filter.
```

One config/manifest belongs to each target. Common options and recursive
filter-specific overrides are supported:

```yaml
input: ../raw/uvis
output: ../redux
instrument: uvis
resume: true
options:
  exclude: []
  background:
    mode: auto
    target_radius: 3.0
    dilate: 8
  alignment:
    mode: auto
  drizzle:
    final_scale: 0.04
    final_pixfrac: 0.8
    final_kernel: square
    final_wht_type: EXP
  astrometry:
    method: two_pass
    catalogs:
      cache: true
      refresh: false
      offline: false
      timeout: 60
      retries: 1
      margin_arcsec: 5
      gaia_dr3:
        provider: astroquery.gaia
        apply_proper_motion: true
      legacy_dr10:
        provider: pyvo
        tap_url: https://datalab.noirlab.edu/tap
        table: ls_dr10.tractor
    relative:
      mode: auto
      require_cross_visit_validation: true
    absolute:
      preferred: gaia_dr3
      fallback: legacy_dr10
      require_validated_solution: true
      fitgeometry: shift
      min_matches: 1
    cross_filter:
      enabled: true
      reference_filter: auto
      common_grid: false
    tweakback:
      enabled: true
      preserve_original_wcs: true
    final_validation:
      enabled: true
filter_overrides:
  F275W:
    drizzle:
      final_pixfrac: 0.8
```

Reference-specific residual/offset/systematic criteria are configured under
`options.astrometry.catalogs.{gaia_dr3,legacy_dr10}.validation`, in arcseconds.
They are defaults for evaluation, not measured precision promises. Source
isolation, shape, magnitude, SNR and Gaia quality cuts are configurable under
`quality`. `shift` uses at least six matches; `rscale` requires at least ten.
Final-source detection can use angular FWHM so changing output scale does not
implicitly change the angular detection model.

The target manifest has schema_version 2 and independently records each filter
or filter/setup under `filters`, including effective settings, observing metadata,
selection, grouping, scientific assessment, provenance, artifacts and status.
A failed astrometric solution leaves a provisional image and an explicit
`astrometry_unvalidated` status. It does not create a new science-ready final DRC.
Earlier products remain recoverable and are identified as previous products.
An explicitly selected `method: relative_only` retains the legacy workflow for
comparison and does not claim validated absolute astrometry.

Background and relative-alignment checkpoints, provisional drizzle, absolute WCS,
TweakBack/headerlets, final drizzle and final QA have separate dependencies.
Final pixel-scale/pixfrac changes can reuse upstream work. Reference/solution
changes invalidate TweakBack and final drizzle while leaving background models
alone. Per-filter receipts prevent unrelated filters from being invalidated.
Astrometric cache restoration starts with the relative checkpoint; the same
transformation is never compounded on an already-corrected exposure.

## Background and science review

Existing background modes auto/constant/2d/none/user and optional geometry/flat-aware
amplifier correction are retained. Auto uses a conservative masked constant and
falls back to zero subtraction when unreliable; 2d is explicit. Protect faint
host emission with an appropriately large target mask. HPSKY/HPSKYSUB and
`skysub=False, skyuser=HPSKY` prevent double subtraction. Detector DQ is preserved;
official CR masks are interpreted as good=1 and converted to added DQ bit 4096.
Single exposures are not described as fully CR cleaned.

Absolute fit-star residuals are not independent holdout accuracy. Cross-filter
accuracy inherited from an anchor is distinguished from directly measured
absolute residuals. Significant epoch differences without validated registration,
sparse fields, morphology-dependent centroids, unresolved guiding/CTE artifacts,
source-core rejection and weight deficits require review. No PSF decomposition,
uncertainty covariance product or automatic official-MAST photometry comparison
is performed.

## Recoverable migration

```bash
hstpipe uvis migrate --output redux --dry-run
hstpipe uvis migrate --output redux
```

Migration preflights conflicts and configuration differences before changing files.
Original trees remain in `{target}/.uvis-legacy-*`. FITS files are verified
byte-for-byte and linked rather than duplicated. Subsequent writes replace
inodes atomically so recovery links remain intact. Small text/config/provenance
files are independent and paths are remapped. Valid legacy receipts are converted
to the new per-group cache layout; incomplete groups remain incomplete. Existing
relative-only outputs are preserved, but the new two-pass default must validate
absolute astrometry before claiming an absolute-calibrated final product.
Conflicting filenames, ambiguous setups, differing raw input directories or
incompatible setup-specific configurations cause explicit refusal. Inspect
`migration.json` before deleting recovery trees. No backup is deleted automatically.

## Testing

Ordinary pytest tests use synthetic FITS and mocked catalog providers. Real-service
smoke queries are separate and must report actual success/failure. The opt-in
four-FLC integration test is selected with HSTPIPE_UVIS_INTEGRATION_INPUT;
setting that variable opts into external catalog requests as well as local
processing. Set HSTPIPE_UVIS_INTEGRATION_OFFLINE=1 to test the explicit
unvalidated/provisional-only outcome without catalog access. Scientific validation status
must be inspected separately from execution status.
Erebus uses DrizzlePac 3.11.0, including the supported apply_tweak API.

See [the Erebus validation record](uvis-validation.md) for actual service outcomes,
real TweakBack/final-drizzle measurements, and unresolved scientific checks.

See [mandatory relative registration](uvis-relative.md) for strategies, precision
criteria, configuration and retained diagnostic products.

## Sparse absolute references and final source validation

Shift fitting has no built-in six-reference requirement. One verified external
source supplies only a diagnostic shift; it cannot independently validate
absolute astrometry. Multiple coherent objects must pass their catalog-specific
residual, geometry and empirical stability checks. The Gaia default residual
limit remains 0.08 arcsec and the Legacy default remains 0.20 arcsec. Eligible
cone rows, image coverage, candidate associations, consensus inliers, native-fit
sources and independently withheld validation objects are recorded separately.

The maintained STScI XYXYMatch histogram matcher precedes robust rejection and
TweakReg. The installed TweakReg wrapper cannot construct a two-source catalog
polygon; two-source absolute shifts use tweakwcs's native linear fitter and
DrizzlePac's distortion-aware WCS updater. This fixes rotation and scale. A finite
proposal that worsens an already adequate original WCS is rejected.

Final persistent-source QA uses native separate no-CR drizzles on the exact final
WCS. Diagnostic median references, empirical input/blank-aperture scatter, source
widths and WHT coverage distinguish actual source losses from rejected transient
CRs. Partial cutouts support stars near subarray edges. These diagnostics neither
replace FLC pixels nor modify production CR parameters. Publication is guarded by
accepted astrometry and source/image QA, preserves old products, and marks
HPQASTAT=VALIDATED. Image arrays, WCS and calibration do not change during QA
annotation. Unassessed or failed stages retain reports and review products.

The four-filter J104242 case passed these checks; that does not validate other
targets automatically. Precision cross-filter registration and advanced aperture,
PSF/host or covariance analysis require their own science-specific validation.


For a deliberate final-sampling-only Python resume, the existing
`astrometry.workflow.finalize(..., reuse_absolute=True)` verifies the saved
absolute-stage receipt, corrected provisional WCS and contributing relative FLC
identities before restoring headerlets and running the final pass. It does not
requery or refit absolute astrometry. Modified dependencies stop that resume;
they require a normal astrometric revalidation. Final image/source QA still runs.
