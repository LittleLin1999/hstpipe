# Mandatory UVIS relative registration

The current full-field workflow and catalog-assisted identity stage are described
in [Full-field UVIS astrometry](uvis-full-field.md). Historical regression counts
and any earlier prohibition on catalog *identification priors* below are
superseded by that document. Relative fitting still uses only independently
measured HST centroids; external RA/Dec are never relative-fit constraints.

Every multi-exposure UVIS reduction now measures relative WCS before the normal
provisional drizzle or absolute calibration, including exposures in one visit.
An unresolved relative solution blocks Gaia/Legacy calibration, TweakBack and a
new final science product. Single exposures have no relative pair to validate;
that limitation and the lack of multi-exposure CR cleaning remain explicit.

## Implementation and APIs

`astrometry/relative_engine.py` coordinates candidate WCS and alignment strategies.
`astrometry/sources.py` provides point and segmented/deblended compact-source
catalogs, centroid-error estimates and detector overlays. Morphology is measured
on unsmoothed SCI pixels so smoothing does not turn hot pixels into valid objects.
`astrometry/relative_fit.py` implements acceptance policy around maintained STScI
matching/fitting, rather than replacing their algorithms.

Inspected DrizzlePac/HAP AlignmentTable source finding, relative/2dhist fitting,
and source selection; inspected Grizli prep.tweak_align, tweak_flt,
iterate_tweak_align, align_drizzled_image and process_direct_grism_visit. The
bounded mosaic/reference iteration follows Grizli's use of original exposure
checkpoints, avoiding accumulated corrections. Grizli is not a dependency.

Verified on Erebus: DrizzlePac 3.11.0, tweakwcs 0.9.3 and Photutils 3.0.0.
HAP's underlying `tweakwcs.matchutils.XYXYMatch` supplies the 2D-histogram and
catalog matching; `tweakwcs.linearfit.iter_linear_fit` supplies constrained,
weighted, clipped fits and leave-one-out refits. `TweakReg` is used with custom
catalogs when supported. Sparse/empty chip catalog footprints can make TweakReg
fail; the native matcher/fitter and `drizzlepac.updatehdr.updatewcs_with_shift`
then update all SCI extensions consistently, including source-empty chips.
The centered zero-origin fit adapter was checked against actual native WCS
updates on full-frame and subarray copies, including rotation/scale and all SCI
extensions. Originals remained unchanged.

References inspected:
- https://github.com/spacetelescope/drizzlepac
- https://github.com/gbrammer/grizli
- https://spacetelescope.github.io/hst_notebooks/notebooks/DrizzlePac/align_sparse_fields/align_sparse_fields.html

## Strategies and decisions

1. Measure the active WCS using common distortion-corrected tangent-plane
   coordinates. Retain it when validated. Inspect shared named alternate WCS,
   IDC/OPUS solutions and available embedded headerlets using the same detector
   source positions.
2. Try standard point-source TweakReg with the simplest supported model.
3. Detect and deblend compact sources using Photutils, testing 3.5/6 pixel kernels
   by default. Validate the original again before imposing any correction.
4. Use the maintained HAP histogram matcher and native constrained fitting.
5. Build a no-CR internal reference mosaic and try bounded mosaic/reference
   refinement. Restart fits from the original working WCS and require measured
   progress; the internal mosaic is not presumed accurate.
6. Optionally enable masked phase cross-correlation on individually drizzled
   images on a common distortion-corrected grid. Catalog-source masks and DQ/WHT
   masks constrain correlation. Its offsets only seed source matching and never
   independently authorize WCS changes. This method requires the optional
   `hstpipe[uvis-image-registration]` extra (scikit-image).

The relative reference exposure is chosen using measured usable compact-source
counts. Cross-filter absolute anchors require validated relative registration and
are ranked by measured final-reference uncertainty and match count; filters have
no hardcoded optical preference. Each filter must establish its own relative
registration before inheriting an absolute frame.

Fit models are shift, fixed-scale shift/rotation, and shift/rotation/uniform-scale.
Their mathematical source minima are 1/2/2. Rank, centroid-error-based parameter
uncertainty, clipping, leave-one-out prediction when identifiable, transformation
stability and physical bounds determine support. The shared-source adjustment
compares weighted model evidence before clipping; see Sparse common-source models
below for the current scientific acceptance and validation scopes.
Initial WCS offsets, current residual RMS, clipped post-fit scatter, centroid
measurement precision and transformation uncertainty are recorded separately.
Centroid uncertainty is an approximate source-width/SNR estimate with configurable
point/galaxy floors, not a full covariance propagation.

Acceptance is not controlled by the old universal 0.3-pixel RMS cutoff. It uses a
configurable science ceiling and a centroid-precision budget. Defaults include
1-pixel maximum scatter/current RMS, a 0.5-pixel median-offset ceiling, and
0.25-pixel transformation-stability ceiling; precise point-source measurements
normally impose a tighter scatter limit through their measurement errors.
Leave-one-out errors must also pass. These are configurable criteria, not measured
accuracy promises.

Candidate corrections additionally require improvement and a separate native
AstroDrizzle coadd comparison with CR rejection disabled. FWHM, elongation, source
flux and coverage are recorded for common eligible objects. WCS updates cannot
win this check by flagging real source pixels as CRs. Proposed changes that fail
validation are discarded; the original exposure WCS remains available.

## Configuration

```yaml
astrometry:
  method: two_pass
  relative:
    validate: true
    mode: auto
    methods:
      - existing_wcs
      - tweakreg
      - improved_catalogs
      - offset_histogram
      - iterative_mosaic
    fit_geometry: auto
    max_iterations: 3
    require_improvement: true
    kernel_widths: [3.5, 6.0]
  absolute:
    preferred: gaia_dr3
    fallback: legacy_dr10
    require_validated_solution: true
  tweakback:
    enabled: true
```

`validate: false` is rejected. `mode: original` retains existing WCS but still
requires measurement; it cannot authorize an unvalidated combination. Existing
explicit user/external transformations can initialize the working copies and
must pass the same subsequent validation. Old adapter settings remain readable;
the new relative policy uses the `astrometry.relative` precision settings.

## Outputs, caching and limitations

Per-filter/setup alignment summaries and source/match tables are under
`intermediate/alignment/`; source overlays, residual plots and no-CR shape
comparisons are under `diagnostics/alignment/`. Rejected native shifts, reference
WCS files and fit logs are retained. The target manifest records all strategies,
accepted/rejected alternatives, source counts and uncertainties. Temporary full
frame copies are removed after candidate assessment. No filter directories are
created.

Failure leaves an explicitly flagged
`intermediate/drizzle/{target}_{filter}_relative_inspection.fits` for review.
This diagnostic coadd uses unvalidated relative WCS and no CR rejection; it is not
an accepted science DRC. Prior products with unresolved relative provenance are
archived before replacement. Background models have separate dependencies and
are reused when alignment policy changes.

Validation currently requires direct common-source support with the selected
reference exposure. Chain-only overlap mosaics can remain unvalidated. Sparse
fields, morphology-dependent centroids, guiding smear and fields lacking enough
independent sources may remain unresolved even when native fitting runs.
Gaussian/moment centroid errors and coadd moments are approximate diagnostics;
no advanced PSF model or science photometry validation is implied.

## J1025+1402 regression outcome

The pre-change code was rerun first on both filters. Neither obtained an accepted
absolute solution. The integrated relative-stage trials subsequently evaluated
active WCS, common IDC/OPUS alternatives, available headerlets, TweakReg, two
segmentation kernel widths, histogram matching and iterative mosaic refinement.
Neither filter established validated relative registration. No strategy is claimed
to have solved these observations.

F467M had 12/25/32/22 usable compact-source detections in the widest-kernel anchor
selection; ifpg02e0q was selected as reference. Original-WCS histogram matching
with the improved catalogs produced 1/9/7 candidate matches at kernel 3.5 and
1/2/9 at kernel 6.0, across the three exposure/reference pairs. Some pairs had
insufficient matches and the others failed fit validation.

F621M had 87/82/78/79 measured compact sources; ifpg02ecq was selected as reference.
Original-WCS histogram matching produced 18/5/5 candidate matches at kernel 3.5
and 15/4/6 at kernel 6.0. Having enough candidate matches did not establish a
reliable fit: the fit/leave-one-out/stability policy still rejected the solutions.
These are candidate associations, not certified astronomical matches.

Both manifest entries have status `relative_astrometry_unvalidated`. Flagged
relative-inspection coadds are retained, and absolute alignment/TweakBack/new
science DRC generation are blocked. Existing production products were preserved.
The group runner checks original raw input identities after processing.

A focused re-evaluation using the saved detector catalogs confirmed that the
final original-WCS preservation policy does not silently accept either filter.
Only the focused relative-registration tests were run at completion: 17 passed.
Coverage includes known offsets versus scatter, histogram matching, geometry
counts, nonfinite/unreliable fits, same-visit failure gating, no absolute/TweakBack
on failure, hot-pixel rejection, image-quality acceptance/rejection, preservation
of an original WCS validated by better catalogs, native custom-reference dispatch,
and explicit user initialization. Full-package test repetition was stopped at the
user's request.

Records on Erebus: per-filter `intermediate/alignment/*_alignment.json`,
`*_relative_decisions.json`, source/match tables and diagnostics; execution log
`/tmp/hstpipe-j1025-release-regression.log`. The optional image-correlation path
has not been demonstrated to produce an accepted solution on these observations.
Manual verification of astronomical source identities and/or more informative
reference-assisted detector centroid measurements remains necessary before
claiming scientific success for J1025.


## Sparse common-source models (2026-11-02)

The standard verified-HST common-source path uses mathematical fit minima:
shift 1, fixed-scale rotation+shift 2 distinct objects, uniform-scale
similarity 2 distinct objects, and affine 3 non-collinear objects. These
counts establish rank only. Affine is available to the numerical assessor,
not automatically added to the UVIS model search. `tweakwcs.linearfit`
provides the weighted shift/rshift/rscale fits. Installed tweakwcs 0.9.3 can
choose a reflected solution when a two-point covariance is rank one; a
proper TAN-plane rotation is then obtained with the maintained
`scipy.spatial.transform.Rotation.align_vectors` fixed-normal constraint.
This constraint fixes geometry, not the positions of any astronomical source.

Candidate models are compared on identical measured HST objects before
outlier clipping, so a coherent rotation is not erased as shift-fit outliers.
Added parameters require configured significance (default 3 sigma), weighted
chi-square improvement, the existing minimum positional improvement, physical
bounds and predictive checks. Rotation and uniform scale are distinct models;
`rscale` is never described as fixed-scale rotation. Weighted information
matrices report nonzero formal errors, parameter correlations and rotation
uncertainty; leave-one-object-out refits are reported only when identifiable.
An exactly constrained fit cannot claim independent validation.

Common source coordinates and exposure parameters are adjusted simultaneously
with SciPy least_squares, anchored to one exposure. All independent native HST
centroids and errors contribute; absent sources add no artificial observations.
Models can differ by exposure. A one-source exposure has translation parameters
only and retains its calibrated rotation/scale. Both SCI chips contribute via
complete HSTWCS to the same undistorted TAN plane. One exposure-level matrix is
propagated through DrizzlePac updatewcs_with_shift(sciext='SCI'), preserving
chip geometry and distortion. Covariance is checked over the actual SCI
footprints, including source-free CCD regions.

Acceptance combines unchanged positional-error limits, graph/frame rank,
residual degrees of freedom, robust rejection and no-CR image morphology,
flux and coverage checks. Existing WCS is retained when measured pair offsets
are adequate. Proposed corrections require measured improvement. Sparse scopes
and unmeasured rotation are explicit; joint training residuals are not
independent validation. The existing target fallback, source detection, DQ,
absolute catalogs, TweakBack implementation and drizzle algorithms are unchanged.
The post-TweakBack relative check reuses these same HST IDs and source-fit
criteria. New models/measurements live in shallow intermediate/alignment;
no-CR comparison FITS live in intermediate/drizzle. HTML includes compact
per-exposure model/error records and retains the strict 1 MB size cap.

Configuration: astrometry.relative.model_significance_sigma (default 3) and
max_rotation_uncertainty_deg (default 0.1) supplement existing precision,
physical-bound, improvement and no-CR QA settings. Centroid uncertainty and
source lever arms, rather than a universal source count, determine whether a
rotation is measurable. A fixed rotation is not reported as a measured zero.

APIs: https://tweakwcs.readthedocs.io/en/latest/source/linearfit.html ;
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html ;
https://drizzlepac.readthedocs.io/en/latest/user_reprocessing/updatehdr.html


### Matched-pair improvement check

`relative_engine.paired_comparison` evaluates the same HST source identities,
row indices, inlier observations and exposure pairs before/after a candidate
correction, across every available overlap. Corrected positions are never
rematched for this statistic. Averaging all pre-fit exposure pairs and only
reference-exposure post-fit pairs is invalid and has been removed. Pair
measurements sharing sources/exposures are correlated, not independent samples.
The configured minimum improvement is unchanged; sparse counts alone do not
authorize corrections. The report shows both RMS values, their difference,
the criterion and the number of pair measurements used.
