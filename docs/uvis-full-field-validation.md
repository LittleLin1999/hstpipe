# Historical J1025 full-field registration investigation

> Historical run record. Counts, paths, thresholds and target outcomes below
> describe those runs, not current defaults or current science-readiness status.
> Use [the UVIS guide](uvis.md) for current policy and each target manifest/report
> for its latest measured outcome. This evidence is retained for reproducibility.

Implemented in the existing editable package on Erebus:
`/quasar/linxiaojing/lrd-squared/codes/hstpipe`.
The focused real-data run documented here used the standard pipeline and finished with
**relative astrometry unvalidated for both F467M and F621M**. No new science DRC,
absolute solution, or TweakBack propagation was authorized for either group.
Successful software tests do not establish successful scientific alignment.

## Implemented workflow

The default sequence is inventory/WCS inspection → catalog identity priors →
full-field HST detection and centroid QC → common measured HST associations →
relative fitting/validation → provisional AstroDrizzle → cross-filter alignment
→ absolute Gaia/Legacy alignment → TweakBack → final AstroDrizzle/QA.

External RA/Dec select approximate search locations and object identities.
Relative fitting uses independently measured native HST centroids transformed
through complete detector WCS models. Projected catalog positions are never
substituted for measured centroids or used as relative-fit constraints.
Catalog failures leave blind HST-to-HST detection/matching available.

### Modules and APIs

- `astrometry/full_field.py`: separately cached full-field DAOStarFinder and
  Photutils segmentation/deblending on every SCI extension, deduplication,
  source IDs, detector statistics and quality tables.
- `centroids.py`, `sources.py`, `matching.py`: common detector/imported-catalog
  QC, unsmoothed morphology, calibrated ERR noise floor, masks, rejection
  reasons, overlays and source cutouts. Masked or edge-clipped essential cores
  cannot pass on high integrated S/N alone. Valid irregular galaxies need not
  have their centroid inside their segmentation footprint.
- `wcs.py`: `stwcs.wcsutil.HSTWCS` plus `all_pix2world` handles detector
  distortion, including applicable SIP/lookup tables; a linear common TAN
  inverse supplies undistorted fitting coordinates. SCI data are not rectified
  before centroid measurement. Identification loads each WCS model once locally.
- `identification.py`: spherical footprint/query margin, per-exposure Gaia epoch
  handling, bounded catalog association, ambiguity/saturation/nondetection
  states and independently measured HST counterparts. A saturated expected star
  cannot be replaced by a detached halo centroid.
- `common_sources.py`: combines external identities with blind HST pair matches,
  retaining HST-only sources. Identity-seeded offsets come from HST measurements.
  Contradictory blind peaks and duplicate-exposure identity components are
  rejected; unverified associations remain explicitly provisional.
- `relative_fit.py`, `relative_engine.py`: established HAP/tweakwcs
  `XYXYMatch`, `iter_linear_fit`, TweakReg custom catalogs and native
  `updatehdr.updatewcs_with_shift`. Existing/alternate/headerlet trials,
  histogram matching, bounded internal-mosaic refinement and optional masked
  image correlation remain available. Connected overlap chains are supported;
  finite parameters, independent source support, propagated uncertainties,
  closure and image QA are required. Reference selection uses measured source
  support rather than the first exposure.
- `catalogs.py`: direct Astroquery/pyvo retrieval, normalization, metadata,
  verified cache subsets, bounded retries/timeouts and failed-query memoization.
  Fixed a real masked-integer Gaia response normalization error.
- `alignment.py`, `workflow.py`: reuse covering identification catalogs with
  absolute-quality cuts; preserve original/corrected provisional WCS; native
  `tweakback.apply_tweak`; check internal source pairs after propagation; final
  AstroDrizzle operates on individual corrected FLC working copies.
- `pipeline.py`, `storage.py`, `config.py`: stage barriers, dependency receipts,
  manifests and compatible configuration. Existing public API and target-level
  output paths remain intact.

Shift-only fitting is preferred; more complex fits require additional spatial
and precision support. The 1.0-pixel RMS/scatter ceilings were not relaxed.
Acceptance also checks outliers, centroid precision, weight concentration,
leave-one-out stability, graph connectivity/closure and no-CR image quality.
A single-source shift can be recorded as a local-only diagnostic; it cannot
authorize independently validated full-field calibration or a science-ready DRC.

Erebus versions: Python 3.11.7, DrizzlePac 3.11.0, STWCS 1.7.7,
tweakwcs 0.9.3, Photutils 3.0.0, Astropy 8.0.2, Astroquery 0.4.11, pyvo 1.9.1.
SEP and Grizli are not installed; Photutils supplies compact-galaxy segmentation.
Neither was added as a mandatory dependency. STIS science was not changed.

## Catalog retrieval and identification

Gaia uses `GaiaClass.launch_job_async`, `gaiadr3.gaia_source`, standard VOTable
output and Astropy epoch propagation. ESA is primary; a supported ARI TAP route
is available. Invalid parallax is not converted into a nonphysical distance.
Integer source identifiers remain exact; masked numeric response columns are
converted before missing-value filling.

Legacy uses `pyvo.dal.TAPService` at `https://datalab.noirlab.edu/tap`, inspected
`TAP_SCHEMA` columns and `Q3C_RADIAL_QUERY` on `ls_dr10.tractor`. Morphology,
mask, flux/uncertainty and available fit/blending columns are retained. This
field's cached measurements are DR10-South, not newly processed DR9-North.

Verified raw cache subsets preserve original query/retrieval/release metadata.
The Legacy cache covers the actual FLCs plus the association search region;
coverage of a larger optional query-margin cone is recorded separately.
Identical failed query/service pairs are not repeatedly requested in a reduction.
Empty results, service failures and HST nondetections are different statuses.
Absolute alignment can reuse a covering cache with its stricter quality cuts.

Live smoke tests: the official ESA Gaia service successfully returned one source
near the supplied bright star, with degree-valued RA/Dec, at
2026-11-01 15:58:03 UTC. Fresh full-field Gaia results also became available.
A fresh small Legacy query timed out after 30 seconds; the ARI live smoke query
timed out after 40 seconds. Neither timed-out integration is claimed successful.
Legacy real-data identification used verified previously retrieved catalogs
(retrieval 2026-10-31). Endpoint routing and error behavior passed focused mocks.

## J1025 source inventory

For each filter, the actual FLC footprints contain **5 quality-selected Gaia
objects, 19 Legacy PSF objects and 31 Legacy compact candidates**. Catalog
categories overlap; these counts must not be added as distinct stars. Legacy
PSF classification alone does not establish stellar identity.

Full-field results below are unique quality-passing detections, not confirmed
astronomical sources. DAO/segmentation counts give accepted/raw detections.
Methods overlap before deduplication.

| Filter | Exposure | Unique passing QC | DAO | Segmentation 3.5 px | Segmentation 6 px |
|---|---|---:|---:|---:|---:|
| F467M | ifpg02e4q | 3 | 0/1 | 1/5 | 3/7 |
| F467M | ifpg02e2q | 13 | 3/7 | 8/23 | 10/23 |
| F467M | ifpg02e0q | 17 | 2/9 | 11/26 | 16/32 |
| F467M | ifpg02dzq | 8 | 0/7 | 4/20 | 7/22 |
| F621M | ifpg02ecq | 52 | 2/9 | 29/53 | 40/69 |
| F621M | ifpg02eaq | 46 | 0/10 | 27/55 | 40/68 |
| F621M | ifpg02e8q | 49 | 0/4 | 29/62 | 39/74 |
| F621M | ifpg02e6q | 43 | 1/12 | 30/55 | 30/51 |

Detection uses a 5-sigma threshold, DAO FWHM 2.5 native pixels and segmentation
kernels of 3.5/6 pixels; no threshold was decreased to force detections. Rejected
counts equal raw minus accepted per method. Rejection reasons overlap and are
retained per row: incomplete/asymmetric cores, masked rings, cosmic-ray-like
morphology, unstable centroids, poor morphology, low S/N and severe blending.
All SCI extensions are searched. These FLCs have high invalid/DQ-mask fractions,
and visually bright sources may have no valid core centroid. NoDetectionsWarning
does not demonstrate an empty astronomical field.

F467M has three repeatable external identity candidates (one Gaia stellar, one
Legacy point-like, one compact galaxy) on e2q/e0q only. The unified association
counts across e4q/e2q/e0q/dzq are 0/3/3/0. This does not connect all four frames
with sufficiently constrained fits. The measured-source reference was e0q.

F621M has one repeatable Legacy point-like and one compact-galaxy identity on
ecq/eaq only; no repeatable accepted Gaia stellar centroid was established.
Unified HST association-node counts are 8/2/10/6 across ecq/eaq/e8q/e6q, including
blind candidates whose astronomical identity remains unconfirmed. The chosen
reference was e8q. No globally validated exposure graph was obtained.

### The supplied bright star

The star at 10:25:34.0761 +14:00:35.153 is present. Its nearest Legacy object,
`10995547142882769`, is PSF, r≈15.56, Gaia-linked to `3887236262361041792`.
All four original F621M FLCs have core DQ 4352 = saturation 256 + CR 4096,
with SCI peaks approximately 74,886–75,280 electrons. These are original flags.

The previous strict Legacy mask selection hid this useful identification prior.
The identification stage now retains Gaia-linked masked stellar priors and
explicitly reports the saturated HST window. It does not override HST pixel QC
or accept a ring/halo centroid. A validated saturated-star wing/PSF fitting
method is not implemented, so this star does not currently supply four reliable
relative-fit centroids. Its presence is not evidence of enough usable references.

## F621M offset hypotheses and residuals

- Old 15-pair hypothesis: the three masked-core pairs remain excluded. The earlier
  pixel-quality replay left only three old pairs, without establishing all as
  genuine matches. Current identity evidence associates two old pair locations
  with the same external objects. Identity evidence does not rehabilitate an
  invalid centroid or validate the original 15-pair fit.
- Old four-pair hypothesis: none has a shared measured external identity. A
  diagnostic shift fit gives (3.906, 14.719) common-TAN pixels, four mathematical
  inliers, zero clipped outliers, fitted RMS 0.903 pixels and leave-one-out RMS
  1.177 pixels. No-CR before/after QA found **zero eligible significant sources**
  meeting the emission/core-coverage criteria. No trustworthy FWHM or ellipticity
  improvement can be reported. This pair experiment cannot validate the group.
- Current eight-pair blind ecq/e8q hypothesis: initial candidate RMS 5.824 pixels,
  fitted shift (5.720, −0.163), fitted RMS 0.928, robust scatter 0.769,
  leave-one-out RMS 1.069, centroid precision estimate 0.354 and transformation
  stability estimate 0.107 pixels. All eight are mathematical inliers with zero
  clipped outliers; they are not eight confirmed stars. Their median native
  detector displacement is only (0.280, 1.024) pixels, versus a detector-fixed WCS
  prediction of (5.440, −0.964) TAN pixels. This is a detector-pattern warning,
  not proof of a telescope pointing error. Training RMS below one pixel is
  insufficient evidence to accept the correction.

Other attempted models and partial fits are preserved in the decision tables.
None established validated four-exposure registration. Existing WCS, IDC/OPUS
alternates, embedded headerlets, TweakReg/native HAP catalogs, histogram matching
and internal-mosaic refinement did not supply an accepted group solution.
Some sparse chips cannot supply a TweakReg catalog footprint; the native HAP
fallback handles those cases without resetting distortion or source-empty chips.

There is no accepted corrected full-group WCS to compare numerically with the
original. Original working WCS is retained after failed candidate trials.
Original-WCS no-CR inspection coadds and candidate pair diagnostics are retained;
no successful image-level improvement or final astrometric precision is claimed.
The optional single-source fallback was unavailable for both groups because no
unique valid compact source was established across all four exposures.

## Tests and integrity

Across targeted runs, **69 distinct relevant parameterized cases passed**:
11 full-field cases, 16 identification cases, 19 centroid cases, 20 relative
refinement cases, and three absolute-cache/TweakBack dispatch/resume cases.
An initial obsolete network mock failed after error classification improved; it
was corrected to model the actual catalog error and its focused rerun passed.
This was not one monolithic 69-test invocation. No unrelated broad suite ran.

Coverage includes native-pixel/complete-WCS projection and distortion, subarrays,
HST-only sources, ambiguous identities, bad histogram peaks, shared-ID fitting
that ignores catalog RA/Dec, connected chains using installed native STScI
fitting/WCS updates, masked cores/rings, irregular galaxies, saturation/blends,
edge-clipped cores, network failure memoization, cache dependencies, local-only
fallback restrictions, no-CR signal/coverage QA and safe TweakBack resume.

Both four-FLC J1025 regression cases ran through the standard pipeline. Hash
checks report original raw FLCs and existing production products unchanged;
background receipts were also unchanged. Post-relative absolute-alignment calls
were empty, as required by the failed relative gate. The invariant JSON calls
this field `external_calls`; it refers to later absolute calls, not the allowed
earlier catalog-identity queries.

## Outputs, restart and limitations

Validation outputs remain under
`/quasar/linxiaojing/lrd-squared/hst/redux/uvis-validation/J1025+1402/uvis/`.
No filter-specific directories were introduced:

```text
uvis/
  products/                       # no new DRC from these failed runs
  intermediate/alignment/         # catalogs, QC, identities, pairs, WCS decisions
    catalogs/                     # external catalog receipts/provenance
  intermediate/drizzle/           # *_relative_inspection.fits, explicitly unvalidated
  diagnostics/alignment/          # overlays, cutouts, graphs, QA and JSON evidence
  diagnostics/drizzle/
  logs/
  config.yaml
  manifest.json
```

Full-field detector receipts are independent of catalog availability; raw-FLC
identification avoids cyclic invalidation from updated working headers. Changes
to centroid rules/common identities invalidate affected relative/downstream
stages. Absolute-reference changes do not invalidate unchanged background or
full-field detection. Failed group validation cannot authorize final propagation
or overwrite old production products.

Cross-filter/absolute fitting, TweakBack and final DRC are implemented but were
blocked on these real cases. Their scientific performance on J1025 therefore
remains unvalidated. Centroid errors and no-CR moments are approximate; they do
not model saturated PSF wings, complete centroid covariance or drizzle covariance.
Blind catalog candidates and histogram peaks still require genuine astronomical
confirmation and image-level support. The remaining cause is not established as
incorrect calibrated distortion or actual guiding failure.

Repeat the relevant validation safely on Erebus:

```sh
cd /quasar/linxiaojing/lrd-squared/hst
source configs/environment.sh
hstpipe uvis reduce --input raw/uvis --output redux/uvis-validation \
  --target 'J1025+1402' --resume
# One filter only:
hstpipe uvis reduce --input raw/uvis --output redux/uvis-validation \
  --target 'J1025+1402' --filter F621M --resume
```

These commands currently retain an unvalidated status, rather than producing a
new science-ready DRC. See `J1025_uvis_fullfield_summary.json` for per-method
rejection counts, graph candidates and complete fit statistics.
