# Recognizing stars and compact sources

[UVIS guide](uvis.md) · [Full-field detection](uvis-full-field.md) · [Relative fitting](uvis-relative.md)

Catalog positions delimit identification windows, never relative-fit coordinates.
The default profile assessment uses Photutils circular aperture masks on native
FLC pixels at independently measured HST centroids. It does not fit a PSF or move
the measured centroid. Optional `profile_method: psf` retains the previous
`PSFPhotometry` / `CircularGaussianPRF` path for explicit use.

## Aperture concentration

Inner and outer aperture radii default to 0.5 and 2 times the configured detection
FWHM. The ratio of background-subtracted inner to outer flux measures compactness.
It is recorded with its uncertainty, aperture coverage, S/N, local background,
size, ellipticity and the original centroid error. Flux errors propagate ERR
with fractional aperture weights, covariance between the nested apertures and
the shared local-background uncertainty. These are identification measurements,
not calibrated total photometry or a new PSF model.

No universal concentration threshold certifies a star. Existing raw-pixel
morphology/CR checks, saturation, masked-core and centroid-stability checks remain
mandatory. A high concentration can also describe a cosmic ray; a compact galaxy
can resemble a star. Incomplete apertures are flagged rather than filled or
corrected to manufacture missing core flux. Distinct compact neighbors inside
an outer aperture require deblending and cannot be resolved by simple apertures.

When several credible candidates remain, compare aperture concentration and
exposure-normalized flux against an independently measured HST anchor, using
measurement errors and configurable profile tolerances. Ties remain ambiguous;
nearest position and brightness never select a winner. The relative fitter and
its validation still establish whether associations are consistent across FLCs.

The default no longer rejects a real source for disagreeing with an analytic
Gaussian PSF (`qfit`). The optional PSF path retains its own fit diagnostics and
conservative covariance/window checks. Original measured catalogs, detector
identity, mask provenance and rejected candidates are preserved in either mode.

## Field-reference fallback

When global registration fails, common sources across the exposures are ranked
by measured quality, preferring quality-assessed compact profiles over other compact
sources. A unique-detection exposure is preferred as a local anchor, followed
by centroid precision. Science-target-only recovery remains a final option,
and configured manual positions still take precedence.

A single usable source constrains local translation only. Its result is
`local_shift_only`, with `full_field_validated=False` and no independent
post-fit residual scatter. If the measured offsets are consistent with existing
WCS and their uncertainties, no working WCS correction is applied. Candidate
local no-CR coadds may be generated on disposable copies for inspection. These
cannot authorize absolute alignment, TweakBack or a science-ready DRC. Multiple
sources still have to satisfy the existing spatial, fit and image-quality gates.

## Configuration and outputs

The existing `astrometry.source_identification.psf_resolution` namespace is retained
for configuration compatibility. `profile_method: aperture` is the default;
`aperture_radii_fwhm`, `min_aperture_valid_fraction` and
`profile_concentration_fraction` control aperture measurements. It also exposes
`enabled`, and optional-PSF parameters `fit_shape`,
`check_shape`, `xy_bounds`, `fwhm_bounds_factors`, `max_qfit`,
`max_centroid_error_pix`, `max_window_shift_pix`, `duplicate_radius_pix`,
`group_separation_pix`, and HST profile-consistency tolerances. Every candidate
in the bounded spatial window is assessed; count alone never rejects a window.
The original relative-WCS RMS ceiling and acceptance criteria are unchanged.
Profile policy/code changes invalidate association and downstream relative stages;
unchanged full-field detections and background modeling remain independently
cached. External catalogs can be reused offline.

Assessment records retain the compatible filename pattern `intermediate/alignment/*_stellar_psf.json`. Compact per-source
multi-exposure figures use flat `diagnostics/*_stellar_association_*.png` files.
The HTML report embeds a bounded selection, includes numerical centroid/local
offset diagnostics and unresolved reasons, and retains the strict <1,000,000
byte limit. No target coordinates, filters or exposure counts are hardcoded.

## Validation scope

Tests cover aperture centroid preservation without PSF fitting, masked/saturated
cores, cosmic-like pixels, neighbors, concentration-based identity discrimination,
reused image contexts and cache isolation. Existing optional-PSF checks remain.
Real-data profile acceptance does not certify relative astrometry: source
identity, transformation stability and independent image QA remain required.

API: [Photutils aperture photometry](https://docs.astropy.org/projects/photutils/en/stable/user_guide/aperture.html).

## Additional objects and partial exposure coverage

A verified native-pixel stellar measurement can seed the existing XYXYMatch
translation search. Its measured HST-to-HST displacement, including centroid
uncertainty, bounds candidate pair residuals. Other sources need not pass a
contaminated whole-catalog fit before their identities can be examined.
Repeated candidate components must have unique exposure membership, compatible
flux/width/ellipticity, no conflicting external identities, and multi-exposure
cycle support (or a shared identified two-exposure pair). These are candidate
identities, not accepted WCS corrections. External RA/Dec never enters this fit.

An unambiguous repeated component may resolve a previously ambiguous catalog
window. Partial exposure coverage is retained rather than requiring every
object in every exposure. Missing measurements remain missing. When an object
already has three measured exposures and a reliable HST guide, a bounded local
peak search may recover a missing centroid. The expected position only defines
the search window: an independently detected peak, ERR-weighted Photutils
centroid_2dg in two fitting windows, and shared pixel/centroid QC are mandatory.
No nearest/brightest choice, interpolated pixels or projected substitute
centroids are allowed. A detected but unstable center remains rejected.

The local primary-source fit uses additional objects as held-out consistency
checks across every available exposure pair, including objects absent from the
anchor exposure. Their residuals and uncertainties are reported separately;
multiple pairwise checks of the same object are not independent new sources.
They do not authorize a full-field solution or relax existing fit/QA limits.
Code changes invalidate identification/relative/downstream caches without
invalidating unchanged background or full-field detector catalogs.

### Historical regression example

These results describe the earlier implementation check, not current target
readiness. Consult each target report for the latest result.

Eight focused tests passed for repeated HST association, partial coverage,
conflicting labels, held-out checks, actual versus masked missing-centroid
measurement, original-input immutability and the existing local fallback.
The J104242 regression recovered the nucleus in all four F689M exposures and
three F621M exposures. The fourth F621M nucleus is detected but has a
window-dependent centroid. No new full-field scientific validation is claimed.
