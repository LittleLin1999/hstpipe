# Full-field UVIS astrometry

The default workflow is FLC inventory/WCS inspection → catalog identity priors →
full-field native-pixel detection/QC → common HST identities → validated relative
WCS → provisional AstroDrizzle → cross-filter registration → absolute Gaia/Legacy
calibration → TweakBack → final AstroDrizzle and QA.

External coordinates select identities and search windows. Relative fits use
measured HST pixels and HST distortion models, never catalog RA/Dec constraints.

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

All paths remain under `redux/{target}/uvis/`, with no filter directories. Tables,
graphs, masks/cutouts and headerlets use `intermediate/alignment/`; provisional
and flagged inspection coadds use `intermediate/drizzle/`; diagnostics remain
stage-organized. Filter/setup IDs distinguish group files, exposure roots
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
