# Image-based stellar association and local field references

Catalog positions delimit identification windows. They never supply the fitted
relative-frame coordinates. Independently detected native HST positions seed
`photutils.psf.PSFPhotometry` with a pixel-integrated `CircularGaussianPRF`.
Installed Photutils 3.0 signatures were inspected and exercised; the UVIS extra
now requires that supported version. No new matching framework is introduced.

## Association

The full-field stage already merges overlapping DAO and segmentation detections.
An ambiguous window is no longer discarded simply because several candidates
remain. Each candidate receives constrained native-pixel PRF fitting, with
simultaneous fitting for neighbors grouped by `SourceGrouper`. Two fitting
windows assess centroid stability. Width bounds, profile residuals, convergence,
covariance, positive flux and shared pixel/centroid QC determine whether a
stellar interpretation is credible. Converged duplicate centroids are recorded,
not counted as independent objects. Saturation elsewhere in a search window no
longer vetoes an unmasked source core. Saturated/masked cores and incomplete
halos remain invalid; no pixels are replaced or interpolated to obtain a fit.

The analytic Gaussian PRF is a centroid proxy, not a calibrated HST photometric
PSF. Its reduced chi-square can exceed unity for a real star. That limitation
is explicitly recorded: formal radial errors are inflated by sqrt(reduced
chi-square), and the larger of those errors, the existing centroid floor and
the fitting-window displacement is retained. Profile residual limits and
covariance/stability checks still reject unsuitable fits. A calibrated ePSF
would be needed for stronger accuracy/photometry claims.

No nearest/brightest selection is made. A single credible profile can identify
a counterpart. If several real stellar profiles remain, their HST-measured
shape and exposure-normalized flux may be compared with a unique HST anchor
from another exposure in that filter. An unresolved tie stays ambiguous.
Coordinates are not used in this profile discriminator, so an unknown relative
translation is allowed. Matched pixels are then transformed with each SCI's
full HSTWCS for the existing STScI relative fitter.

Fitted HST positions replace detector estimates only in the astrometric working
catalogs. Original detector catalogs, initial positions, fit diagnostics,
detector method/segment provenance and rejection reasons remain available.
Manual/local detection also uses the same resolver for ambiguous stellar
windows. Existing galaxy centroid behavior and masked-core checks are retained.

## Field-reference fallback

When global registration fails, common sources across the exposures are ranked
by measured quality, preferring fitted stellar profiles over other compact
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

`astrometry.source_identification.psf_resolution` exposes `enabled`, `fit_shape`,
`check_shape`, `xy_bounds`, `fwhm_bounds_factors`, `max_qfit`,
`max_centroid_error_pix`, `max_window_shift_pix`, `duplicate_radius_pix`,
`group_separation_pix`, and HST profile-consistency tolerances. Every candidate
in the bounded spatial window is assessed; count alone never rejects a window.
The original relative-WCS RMS ceiling and acceptance criteria are unchanged.
PSF policy/code changes invalidate association and downstream relative stages;
unchanged full-field detections and background modeling remain independently
cached. External catalogs can be reused offline.

Fit records use `intermediate/alignment/*_stellar_psf.json`. Compact per-source
multi-exposure figures use flat `diagnostics/*_stellar_association_*.png` files.
The HTML report embeds a bounded selection, includes numerical centroid/local
offset diagnostics and unresolved reasons, and retains the strict <1,000,000
byte limit. No target coordinates, filters or exposure counts are hardcoded.

## Focused validation

Sixteen focused tests cover stellar versus CR-like profiles, ambiguous and crowded stellar windows,
joint overlapping fits, failed validation-window fits, duplicate detector measurements, inherited-4096
reassessment, true masked/saturated cores, a saturated neighbor, manual recovery,
HST-only profile discrimination, automatic field-star fallback, stage-cache
dependencies, failed-stage propagation and HTML size enforcement. No broad
test suite was run. J104242 F621M/F689M use the standard reduction entry point
with cached catalogs/background/detections; local results do not certify global
astrometry. Original FLC and prior production-product hashes are checked.

API reference: [Photutils PSFPhotometry](https://docs.astropy.org/projects/photutils/en/stable/api/photutils.psf.PSFPhotometry.html).

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

Eight focused tests passed for repeated HST association, partial coverage,
conflicting labels, held-out checks, actual versus masked missing-centroid
measurement, original-input immutability and the existing local fallback.
The J104242 regression recovered the nucleus in all four F689M exposures and
three F621M exposures. The fourth F621M nucleus is detected but has a
window-dependent centroid. No new full-field scientific validation is claimed.
