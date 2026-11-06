# UVIS PSF generation

`hstpipe uvis psf` packages the empirical-v2 and product-matched drizzle workflow
previously run as `refine_mast_psfs.py` and `drizzle_uvis_psfs.py` on Erebus.
It is a separate operation after science reduction, not an automatic part of
`hstpipe uvis reduce`. Existing STIS PSF commands are unchanged.

## Installation and execution

```sh
python -m pip install -e '.[uvis-psf]'
hstpipe uvis psf --config psf.yaml --stage empirical
hstpipe uvis psf --config psf.yaml --stage drizzle
hstpipe uvis psf --config psf.yaml --stage validate
# Or run all three stages:
hstpipe uvis psf --config psf.yaml
```

Example `psf.yaml` for the existing Erebus data:

```yaml
workdir: /home/linxiaojing/data/lrd-squared/hst/uvis_psf
chip: 2
filters: [F275W, F336W, F467M, F547M, F621M, F689M, F763M, F845M]
products: [hstpipe, mast]
models: [empirical_v2, anderson]
goal: 100
limit: 800
```

Defaults relative to `workdir` retain the original directory layout:

| Option | Default | Purpose |
|---|---|---|
| `targets` | `targets.json` | Science positions used to rank archival stars |
| `manifest` | `drizzled/input_products/manifest.json` | Science cutouts and their provenance |
| `raw_root` | `../raw/uvis` | Original FLCs for MAST-product replay |
| `anderson_library` | `anderson_comparison/library` | Local exact-filter `STDPSF_WFC3UV_FILTER.fits` models |

Explicit relative path values resolve against the YAML file's directory.
If omitted, `products` is `[hstpipe]` and `models` is `[empirical_v2]`.
Omitting `filters` selects the filters present in the stage's input metadata.
The empirical stage needs `targets`; drizzle/validate need `manifest`.
The optional `--stage pilot` processes the saved `review_sample` catalog
(default `workdir/review_sample_manifest.json`) into `refined_v2/pilot.json`
for the retained mask-review reports. It does not rebuild stacked PSFs.
Anderson-only runs (`models: [anderson]`) do not need empirical stacks or targets.
Library downloads are not automatic. In a combined run, absent exact-filter
Anderson models are recorded as skipped; Anderson-only runs fail if missing.

The command writes a timestamped configuration/environment/source-code receipt
under `runs/`. A workdir lock prevents simultaneous writers. Errors and
insufficient usable stars return nonzero status; details remain in stage
summaries. Cached cutouts and queries are reused. The package's processing
version differs from the old scripts, so old processed decisions are rebuilt
on the first empirical run. No old v1 PSF or diagnostics are required.
Use the drizzle stage alone to retain and reuse existing native v2 models.

## Inputs

`targets.json` is a list with exact filters and zero-based full-detector UVIS2
coordinates (including the subarray LTV offset correction):

```json
[
  {"target": "example", "filter": "F275W", "exposures": [
    {"detector_x": 500.2, "detector_y": 512.7, "chip": 2}
  ]}
]
```

The science manifest is the cutout-export format used by the existing workflow:

```json
{"records": [{
  "target": "example", "filter": "F275W",
  "center_ra_deg": 150.0, "center_dec_deg": 2.0,
  "products": {"hstpipe": {
    "path": "example_F275W_drc_cutout_10arcsec.fits",
    "source_resolved": "/data/reduced/example/uvis/products/example_F275W_drc.fits",
    "cutout": {"sha256": "SHA256_OF_THE_CUTOUT"}
  }}
}]}
```

Cutout paths resolve against the manifest directory. Each cutout must retain
the original DRC primary header, SCI/WHT/CTX and HDRTAB extensions and contain
the full 61×61 output stamp. The full detector FLCs, distortion reference data,
and hstpipe final rejection masks must remain accessible. By default hstpipe
working inputs are located at `source_resolved.parent.parent/intermediate/working`;
an absolute `working_directory` in a product record can override this after a
relocation. Final CR masks must be in its sibling `drizzle` directory.

Creating this science-cutout manifest is an input-preparation step, not a new
science reduction. The existing project `export_uvis_cutouts.py` produces it.
This release does not automatically discover arbitrary reduction layouts.

## Algorithm and products

The empirical stage queries exact-filter, full-frame UVIS2 stars from MAST,
first in the original near-detector search region (26–1000 pixels), then
chipwide. It ranks by detector distance, catalog fit quality and brightness;
the search is bounded by `goal` distinct clean stars and `limit` attempted
detections, not exhaustive. It downloads SCI/ERR/DQ, masks all nonzero DQ and
invalid pixels, segments/deblends neighbors and supplements their masks with
the exposure catalog. It fits a background plane and a Gaussian core centroid,
then uses compact bilinear registration with strict mask propagation.

Accepted stars have aperture S/N >=100, adequate unmasked coverage and acceptable
core shape. One detection per 0.3-arcsec sky group is retained. Entire exposure
roots are held out when sample size permits. Flux-scaled, clipped median stacks
are cropped to 41, 35, 31, 25 or 21 pixels according to contributor coverage.
`refined_v2/psfs/MAST_UVIS2_FILTER_v2_trial.fits` contains the unit-sum PSF,
NCONTRIB, scaled-MAD SCATTER, and optional held-out residual. Selection,
rejections, holdout metrics and convergence records remain in `diagnostics/`.

The drizzle stage reconstructs exposure WCS/SIP from HDRTAB, retains lookup
distortions and uses the original full-detector mapping grid. It renders the
point source at each target phase, then applies the science masks, exposure
weights, kernel and pixfrac through DrizzlePac. Hstpipe uses retained
`final_crmask` plus calibration DQ and HPDQSCI; MAST uses FLC DQ ignoring bits
336. Static masks are reconstructed from the current input arrays. WHT replay
must agree within the recorded tolerances and CTX must agree exactly before
PSFs are written. Unsupported weight/mask recipes fail this check.

Outputs are `drizzled/products/TARGET/TARGET_FILTER_PRODUCT_psf_MODEL.fits`
and matching provenance JSONs. Primary arrays are unit-sum 61×61 stamps;
WHT and CTX describe combination coverage. The validation stage checks source
model/cutout/product hashes, FITS checksums, normalization, shapes and replay
metrics, and writes `drizzled/validation.json` and `drizzled/README.md`.

## Scientific and implementation limits

- UVIS2 only; UVIS1 and multi-plane CTX (>32 input entries) are not supported.
- Empirical stacks are native-pixel trials, not oversampled ePSFs. Bilinear
  registration broadens their cores; cubic phase interpolation is approximate.
- Anderson models are pixel-integrated ePSFs evaluated at native pixel centers;
  no additional 4×4 pixel integration is applied. They require exact-filter files.
- Drizzle geometry is matched, but target focus, stellar color and epoch are not.
- Unit sums refer to finite stamps, not infinite-aperture photometry. Keep the
  fractional reference point in PSFCENX/Y rather than assuming (30,30).
- WHT is not a PSF uncertainty map; SCATTER is not a standard error. Heldout
  agreement is internal validation, not validation against a science galaxy.
- Science inputs are read-only. Generated PSF products and current stage
  summaries in the explicitly selected workdir are replaced when rerun.
- MAST queries require network access; the full drizzle replay requires the
  retained detector inputs and the `uvis-psf` dependencies, typically on Erebus.

## Verified reproduction on the LRD UVIS2 dataset

On 2026-11-06 the package rebuilt all eight native empirical-v2 filter stacks
from the retained MAST query/cutout caches in a separate workdir. Every numerical
FITS plane and every training/holdout member selection matched the originals
exactly. It then reproduced all 39 empirical-v2 drizzled PSFs (19 hstpipe and
20 MAST science products across seven sources), plus all 22 Anderson alternatives.
PSF, WHT and CTX arrays matched exactly. All 69 original PSF files remained
unchanged and were verified against their local copies by SHA256.

The production-environment regression suite passed 47 tests without skips.
This verifies reproducibility of these retained inputs, not additional focus or
photometric calibration. Preserve the exact query and cutout caches for repeatable
selection; a fresh query against a changing archive need not select the same stars.

Project evidence is retained under
`hst/uvis_psf/hstpipe_reproduction_20261106/` on Erebus and
`hst/analysis/reports/uvis_psf_integration/` in the local analysis workspace.
