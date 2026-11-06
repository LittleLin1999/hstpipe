# Full-field UVIS astrometry

[UVIS guide](uvis.md) · [Fit models and acceptance](uvis-relative.md) · [Local association](uvis-stellar-association.md)

The standard workflow detects sources on **every available SCI extension**, then
uses their native HST centroids to register exposures. Both CCDs contribute when
usable sources exist; single-chip subarrays use the same interface.

| Stage | Measurements used |
| --- | --- |
| Identify sources | HST detections, optionally labeled with Gaia/Legacy candidates. |
| Register exposures | HST-to-HST measurements transformed through complete chip-specific WCS. |
| Align the mosaic absolutely | Gaia/Legacy positions, after relative registration passes. |
| Build the final image | Corrected working FLCs after TweakBack; a new AstroDrizzle pass. |

Catalog coordinates define search windows and identities. They never replace
measured HST centroids in the relative fit. Network failure does not prevent
independent full-field HST source detection.

## Modules and APIs

- `astrometry/full_field.py`: separately cached DAOStarFinder and Photutils
  segmentation/deblending across every SCI extension, detector-method
  deduplication and per-exposure/per-chip detection statistics.
- `centroids.py`: shared QC for detector/imported catalogs. Rejects invalid or
  clipped essential cores, masked rings/halos, spikes, saturation, severe blends,
  unstable/nonfinite centroids and poor precision. A valid galaxy centroid need
  not belong to its segmentation footprint. Rejected rows and cutouts are kept.
- `wcs.py`, `sources.py`: complete STWCS `HSTWCS` detector-to-sky transformation,
  including SIP and lookup tables, followed by a linear reference TAN inverse.
  Native SCI pixels are not resampled before centroid measurement. Identification
  WCS models load once per exposure and are not reused across WCS updates.
- `identification.py`: spherical multi-SCI footprint, Gaia/Legacy cross-identity,
  bounded association with independent full-field detections, saturation and
  ambiguity records. HST-only sources remain available.
- `common_sources.py`: external labels and blind HST pair associations. Shared
  identities seed offsets using HST measurements only. A blind histogram peak
  cannot override contradictory measured common-object positions.
- `relative_fit.py`, `relative_engine.py`: maintained HAP components
  `tweakwcs.matchutils.XYXYMatch` and `linearfit.iter_linear_fit`, TweakReg custom
  catalogs, alternate/headerlet trials, histogram matching, bounded internal
  mosaic refinement and optional masked image correlation.
- Native overlap-graph traversal supports chains without direct anchor matches.
  `drizzlepac.updatehdr.updatewcs_with_shift` updates every SCI on candidate copies.
  Finite fits, enough independent distributed matches, LOO stability, physical
  bounds, propagated uncertainties, connectivity and cycle closure are required.
- `alignment.py`, `workflow.py`: covering catalog reuse, absolute TweakReg,
  STWCS WCS archival/headerlets, `tweakback.apply_tweak`, fresh final AstroDrizzle
  and final validation. Absolute fitting does not conceal failed relative WCS.
- `pipeline.py`, `storage.py`, `config.py`: stage barriers, shared paths,
  independent cache dependencies, manifest/config and unchanged public API/CLI.

The existing 1.0-pixel scatter/current-RMS ceilings were not increased. Shift is
the first model; extra degrees of freedom require more sources and spatial
support. Training RMS alone cannot accept a solution.

No-CR coadds assess significant source emission, width, elongation, flux and
EXP core coverage. Positive noise fluctuations or weight deficits cannot approve
an apparent improvement. Before-correction centroids may be displaced; corrected
centroids must be consistent. Common-ID QA samples the connected field with one
entry per object. These moments are not a PSF model or calibrated noise covariance.

A single target-region source can produce a `local_only` native shift diagnostic.
Its independent scatter is undefined, not zero. It cannot authorize full-field
absolute calibration, TweakBack or a science-ready DRC. An optional diagnostic
no-CR coadd uses temporary copies only.

## Catalogs and failures

Gaia uses Astroquery `GaiaClass.launch_job_async` and standard VOTable output.
ESA is primary; ARI Gaia DR3 TAP is a supported alternate. Requests have bounded
timeouts/retries, and identical failed query/service pairs are remembered for
one reduction. Normalization safely converts masked integer columns before
filling floating missing values, preserves integer source IDs exactly, and
propagates positions to actual HST epochs using Astropy. Nonphysical parallax
distances are not used.

Legacy uses pyvo, `https://datalab.noirlab.edu/tap`, `ls_dr10.tractor`,
`Q3C_RADIAL_QUERY` and inspected `TAP_SCHEMA` columns. DR10-South and DR9-North
are distinguished. Legacy PSF classification alone does not establish a star.

Verified raw cache subsets can be normalized without new requests. Their original
query, retrieval date, release and coverage remain explicit. A covering cached
Legacy query is preferred. Partial query-cone coverage is distinguished from
complete FLC-plus-search-window coverage. Service failure, normalization failure,
empty results, nondetection, rejected detection and ambiguity are separate states.
Unavailable catalogs do not disable blind full-field HST registration.

Gaia-linked bright PSF objects can remain identification priors despite bad survey
photometry. This does not override HST QC or make them absolute references. A
saturated star is not replaced by a detached halo centroid. Validated saturated-
star PSF-wing fitting is not implemented.

## Configuration

Existing nested detector settings remain supported. `relative.detection:
full_field` is also accepted as an alias for default full-field settings.

```yaml
options:
  astrometry:
    method: two_pass
    workflow: relative_then_absolute
    source_identification:
      enabled: true
      catalogs: [gaia_dr3, legacy_dr10]
      use_external_positions_for_relative_fit: false
      prefer_stars: true
      allow_compact_galaxies: true
      allow_masked_stellar_priors: true
      search_radius_arcsec: 2.0
    catalogs:
      cache: true
      refresh: false
      offline: false
    relative:
      detection:
        snr: 5
        fwhm_pixels: 2.5
      kernel_widths: [3.5, 6.0]
      fit_geometry: auto
      max_scatter_pix: 1.0
      max_current_rms_pix: 1.0
      local_fallback:
        enabled: true
        diagnostic_coadd: false
    cross_filter:
      enabled: true
      reference_filter: auto
    absolute:
      preferred: gaia_dr3
      fallback: legacy_dr10
      require_validated_solution: true
    tweakback:
      enabled: true
```

Detection uses per-chip sigma-clipped background/noise with a calibrated ERR
noise floor and local centroid/pixel checks. Smoothing is detection-only;
morphology uses unsmoothed SCI. Thresholds are not lowered to force detections.
SEP was unavailable on Erebus; Photutils segmentation supplies the extended-
source path. Centroid uncertainties are approximate and include configurable
floors, not complete PSF or pixel-noise covariance propagation.

## Restart and outputs

Raw FLCs are read-only. Background, full-field detection, identification, relative
WCS, provisional drizzle, absolute WCS, TweakBack and final drizzle have separate
receipts. Catalog changes reuse unchanged detector catalogs and backgrounds.
Changed source quality or relative WCS invalidates affected descendants. Per-filter
receipts allow selected-filter reruns without reducing unrelated filters.

All paths remain under `redux/uvis/{field}/`, with no filter directories. Tables,
graphs, masks/cutouts and headerlets use `intermediate/alignment/`; provisional
and flagged inspection coadds use `intermediate/drizzle/`; diagnostics remain flat, with stage-prefixed filenames. Filter/setup IDs distinguish group files, exposure roots
distinguish detector files. A shared config/manifest records separate filters.

Only validated relative groups proceed to the normal provisional product.
Cross-filter matches must independently pass validation. Accepted absolute WCS
is propagated to fresh individual FLC copies with TweakBack, retaining relative
registration. Final AstroDrizzle resamples those FLCs, not provisional SCI.
Internal source pairs are rechecked after TweakBack. Final native SCI/WHT/CTX and
HDRTAB where available retain units, calibration and provenance. Background/CR/
photometry are preserved, with no extra PAM or fabricated ERR extension.

Unvalidated groups retain diagnostics and do not proceed to absolute fitting,
TweakBack or overwrite production DRCs. Final fit-star residuals are not described
as independent held-out absolute precision.

## References

- https://docs.astropy.org/projects/astroquery/en/stable/api/astroquery.gaia.GaiaClass.html
- https://gaia.ari.uni-heidelberg.de/tap.html
- https://github.com/spacetelescope/drizzlepac
- https://spacetelescope.github.io/hst_notebooks/notebooks/DrizzlePac/align_sparse_fields/align_sparse_fields.html
- https://github.com/gbrammer/grizli


## Direct reference catalogs

### Gaia DR3

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

### Legacy Surveys

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

### Footprints, caches, and service failures

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


## Compact references and source-selection provenance

Point-source detection (`astrometry.detection.method: dao`) remains the default
for absolute astrometry. For compact galaxies or nuclei missed by DAO detection,
`compact_stable` reuses the existing Photutils segmentation and mandatory
pixel/centroid QC at two kernel widths. A source must pass both measurements,
match uniquely on the same SCI extension, meet HST S/N and centroid-stability
requirements, and respect the configured catalog size limit. Neither catalog
coordinates nor fitted residuals decide which morphology passes these cuts.
Full-field relative detection and its target-specific fallback remain unchanged.

Example explicit filter override (replace F689M as appropriate):

```yaml
filter_overrides:
  F689M:
    astrometry:
      detection:
        method: compact_stable
        compact_kernel_widths: [3.5, 6.0]  # pixels on the image being measured
        compact_min_snr: 10
        compact_max_centroid_change_pix: 0.5
      source_identification:
        catalog_quality:
          legacy_dr10:
            min_snr: 10
            isolation_arcsec: 1.5
            max_fracflux: 0.5
            max_fracmasked: 0.2
      catalogs:
        legacy_dr10:
          quality:
            min_snr: 0
            max_position_error_mas: 200
```

`source_identification.catalog_quality` optionally supplies a separate quality
policy for each identification catalog. A supplied policy replaces that stage's
shared quality settings (provider defaults still apply). Without it, existing
shared-policy behavior is preserved. Use it when changing absolute-reference
selection so a fresh run retains the intended relative-identification policy.
The effective cuts are recorded in each query's provenance and affect its cache key.

This optional catalog policy replaces the Legacy flux-S/N cut with a finite,
positive radial positional-error limit, while retaining coverage, primary-object,
mask, blending, isolation and compact-morphology cuts. It does **not** relax the
0.20-arcsec Legacy residual requirement or bypass independent HST centroid QC.
The 200-mas example is a selection ceiling, not a precision claim. Do not enable
it automatically for all fields: faint, blended or morphology-dependent catalog
centroids still require review. Missing uncertainties fail this explicit policy.

Audit tables distinguish pixel-quality acceptance from final source selection;
rejected objects remain recorded. Absolute and final validations record
`source_quality` and `source_selection` separately. Reports describe the active
selection directly from configuration/stage results, without a target-specific
resume script. Existing manifests with the older descriptive field remain readable.

The J012930 F689M compact-only run retained two Legacy counterparts and obtained
0.032 arcsec final RMS and 0.056 arcsec leave-one-out RMS. These are sparse fit-source
checks, not independent held-out accuracy or a proof that other targets pass.
Its source/CR checks and TweakBack completed; local-relative/full-field limitations
remain explicitly distinct. No filter-specific thresholds or source IDs are
hardcoded in the package.


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

The no-CR source-QA reference resets only AstroDrizzle DQ bit 4096 on disposable
copies before separate drizzle. It preserves calwf3 bit 8192, other detector
flags and SCI values; production CR masks cannot bias both sides of this check.


## Complete configuration example

Save a YAML file such as `configs/uvis.yaml`. The settings below show the main
controls; most reductions can start with the defaults. The final drizzle scale
is in arcsec/pixel, and pixfrac is dimensionless.

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
`quality`. There is no built-in six-reference or ten-reference requirement.
Mathematical fit minima are separate from scientific acceptance; one absolute
reference is diagnostic only. See the sparse-reference criteria below.
Final-source detection can use angular FWHM so changing output scale does not
implicitly change the angular detection model.

### Processing records

The target manifest has schema_version 2 and independently records each filter
or filter/setup under `filters`, including effective settings, observing metadata,
selection, grouping, scientific assessment, provenance, artifacts and status.
A failed astrometric solution leaves a provisional image and an explicit
`astrometry_unvalidated` status. It does not create a new science-ready final DRC.
Earlier products remain recoverable and are identified as previous products.
An explicitly selected `method: relative_only` retains the legacy workflow for
comparison and does not claim validated absolute astrometry.

### What a resume reuses

Background and relative-alignment checkpoints, provisional drizzle, absolute WCS,
TweakBack/headerlets, final drizzle and final QA have separate dependencies.
Final pixel-scale/pixfrac changes can reuse upstream work. Reference/solution
changes invalidate TweakBack and final drizzle while leaving background models
alone. Per-filter receipts prevent unrelated filters from being invalidated.
Astrometric cache restoration starts with the relative checkpoint; the same
transformation is never compounded on an already-corrected exposure.

When relative registration is established through a connected overlap graph,
post-TweakBack validation rechecks its accepted exposure pairs, their original
centroid measurements, and closure consistency. It does not require new direct
connections to every exposure. Disconnected graphs or degraded accepted edges
remain unvalidated; previously rejected edges remain in the diagnostic record.
