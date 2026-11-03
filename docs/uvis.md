# UVIS user guide

[Home](../README.md) · [Installation](INSTALLATION.md) · [Code architecture](ARCHITECTURE.md)

Use this guide to reduce WFC3/UVIS FLC exposures, understand the report, and
publish an image after review. Original downloaded FLCs are never modified.

- [Run a reduction](#run-a-reduction)
- [Read the result](#read-the-result)
- [Review and publication](#review-and-publication)
- [Output layout](#output-layout)
- [Cleanup](#cleanup)
- [Configuration](#configuration-and-invocation)
- [Technical reference](#technical-reference)

## Run a reduction

Activate the installed environment and work from your observations directory.
Set `iref` to your UVIS reference directory, including the trailing `/`.
Missing header-named references are downloaded automatically through HST CRDS.
See [reference locations and offline settings](INSTALLATION.md#download-calibration-references).

```bash
export iref=/absolute/path/to/observations/calibration/uvis/
hstpipe uvis inspect --input raw/uvis
hstpipe uvis reduce --input raw/uvis --output redux --dry-run
hstpipe uvis reduce --input raw/uvis --output redux --resume
```

For one filter or a saved configuration:

```bash
hstpipe uvis reduce --input raw/uvis --output redux \
  --target YOUR_TARGET --filter F336W --resume
hstpipe uvis reduce --config configs/uvis.yaml --resume
```

`--dry-run` checks selection and prints the plan; it does not establish science
readiness. `--resume` verifies cached outputs before reusing them. YAML paths are
relative to the YAML file; command-line paths are relative to the current directory.

### What the pipeline does

| Step | Purpose |
| --- | --- |
| Validate inputs | Check calibration, geometry, and every SCI/ERR/DQ set. |
| Assess background and detect sources | Measure reliable native-pixel centroids on all available CCDs. |
| Register exposures | Validate existing relative WCS; refine when justified, with a local fallback if needed. |
| Make a provisional mosaic | Run AstroDrizzle with the accepted relative registration. |
| Establish absolute astrometry | Use Gaia DR3 or Legacy Surveys; transfer between filters only when validated. |
| Propagate and drizzle again | Use TweakBack on working FLCs, then a new final AstroDrizzle pass. |
| Validate and report | Check final WCS, source preservation, CR rejection, and coverage. |

External catalogs can identify sources before relative registration. Their RA/Dec
positions do not constrain the relative HST-to-HST fit. No calwf3 rerun, extra CTE
correction, or second pixel-area-map correction is applied.

## Read the result

Open `products/{target}_uvis_report.html` inside the target's reduction directory.
It is self-contained, works offline, and stays below **1,000,000 bytes**.
Start with the filter summary, then open only the filters needing attention.

| Report status | Meaning | Next action |
| --- | --- | --- |
| Ready | Astrometric and image checks passed. | Use the DRC with its recorded limitations. |
| Ready with sparse-source limits | Checks passed with limited field constraints. | Check whether those limits suit your science. |
| Ready — human reviewed | An explicit reviewer accepted documented field-source flux flags. | Retain the caveats with the data. |
| Review image QA | Astrometry passed; source-preservation checks need review. | Inspect the flagged-source figures below. |
| Blocked | A required scientific check did not pass. | Address the stated cause; do not override it through a flux-only approval. |
| Pending / Running / Processing error | Processing is incomplete or encountered an error. | Read the stage message and resume after resolving it. |

A file left from an earlier run does not establish that the current run passed.
Use the report and manifest together. A single-exposure image cannot claim
reliable multi-exposure cosmic-ray cleaning.

## Review and publication

**Review image QA** is different from failed astrometry. First check the report's
relative, post-TweakBack and final absolute WCS statuses. A field-source review
cannot waive a failed astrometric solution, guiding problem, target-shaped WHT
hole, or other unresolved image-quality failure.

### What to inspect

For each flagged source, use the report's visible review figure and measurements:

1. Identify the source and its distance from the science target; do not assume
   a field-source warning proves that target flux was lost.
2. Compare the individual no-CR exposure panels at the fixed reference cross.
   Look for source doubling, unstable centroids, residual CRs and masked cores.
3. Compare the no-CR median reference, final SCI, difference map and WHT. Use the
   flagged aperture radius; an agreeing larger aperture does not erase a small
   aperture warning. SCI panels use a shared display scale; residual/WHT scales
   are separate. These individually drizzled panels are not raw detector arrays.
4. Inspect detector CR masks if real-source rejection is suspected. Median versus
   weighted combination, pixel sampling and centering can also change core flux.
   Visual agreement is not independent proof of precision photometry.
5. Decide whether the remaining field-source discrepancy is acceptable for the
   intended science, and record the source IDs and scientific justification.

### If the review fails

If the issue is unacceptable or affects the science target, keep the image under
review. Correct the identified cause and rerun only affected stages, then repeat
QA. Do not multiply the science image by an arbitrary flux correction, lower an
accuracy threshold to force a pass, or edit FITS status keywords manually.

### If the review passes

If astrometry and other image checks already pass, an explicitly reviewed
**field-source flux-only exception** can be accepted without another drizzle:

```bash
# Preview the current image and flagged sources; no approval is made.
hstpipe uvis review --directory /data/redux/TARGET/uvis --filter F336W

# After inspecting the current report, accept the listed field-source flux flags.
hstpipe uvis review --directory /data/redux/TARGET/uvis --filter F336W \
  --approve --reviewer "Your name" \
  --reason "Reviewed the field-source core-flux discrepancy; acceptable for my target analysis."
```

Use the actual target reduction directory containing `manifest.json`, not a
report-only merged view. The command finds the recorded image and source IDs,
and verifies the saved file checksum internally. **You do not need to calculate
or copy a checksum.** If the image has changed or its recorded identity is
missing, approval stops; inspect and validate the new assessment first.

Approval applies to all listed field-source flux flags in this filter. The same
scientific checks still refuse failed astrometry, target-region issues, coverage
or morphology flags, and unrelated image warnings. This command does not waive
those failures or start another reduction.

The lower-level Python API `approve_source_review` still accepts
`expected_sha256` for scripts that need to bind approval to a specific image.
Ordinary users should use the command above.

### What approval changes

On approval, the existing SCI/WHT/CTX arrays, WCS, units and calibration remain
unchanged. The product is published as `products/{target}_{filter}_drc.fits` with
`HPQASTAT=REVIEWED`, `HPREVIEW=T`, and reviewer/date/reason/source provenance in
FITS HISTORY. `ASTSTAT` and `ASTREL` retain their independently established values.
The original `source_validation.status=review_required` and all measurements are
preserved. The manifest records `source_review`, and the HTML says **Ready —
human reviewed** with the figures and warnings still visible. The original
review FITS remains an immutable audit input on the reduction server.

Approval is limited to the exact image and source-QA assessment; it is not a
new default masking policy and does not certify target precision photometry.
Repeated approval of the same unchanged record is idempotent. It does not rerun
source detection, alignment, TweakBack, drizzle or photometry. A normal subsequent
reduction or changed science dependencies must be revalidated; this receipt is
not a blanket approval of future products. Accepted products can be synced to
local `products/`; remove superseded local review copies only after verifying
the new file's checksum. Keep the original review evidence on the server.

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

## Cleanup

Preview files that can be removed, then apply the same policy:

```bash
hstpipe uvis cleanup --output redux
hstpipe uvis cleanup --output redux --apply
```

`--dry-run` is an explicit synonym for preview. Add `--target YOUR_TARGET` to
limit the operation. `--output` also accepts a target's `uvis/` directory or a
run directory containing several `redux/` roots, including retries:

```bash
hstpipe uvis cleanup --output /path/to/hst/runs/RUN_NAME --apply
```

The command removes only these full-size comparison images when a matching group
and exact input set have a checksum-verified published DRC:

- `*_joint_existing_inspection.fits`
- `*_joint_candidate_inspection.fits`
- `*_relative_inspection.fits`

It keeps final DRCs, reports, source tables, numerical QA, compact plots,
provisional mosaics, working FLCs, catalogs, CR masks, and original inputs.
`*_unvalidated_final.fits` files remain because they can be inputs to recorded
QA or human approval. Diagnostics for unpublished groups remain available.

Each affected target receives a cleanup receipt in `logs/`, with paths, sizes,
and checksums. Aliases to deleted comparisons inside intermediate drizzle
folders are removed too. Scientific statuses and cache receipts are unchanged.
**A later resume may regenerate deleted diagnostics** if a receipt requires them;
other stage checkpoints remain available. This command does not delete archive
directories based on their names or age.

Do not delete an entire `runs/` directory just because it looks temporary:
published products may be symlinks into it. Check their destinations first.

## Configuration and invocation

Save this as `configs/uvis.yaml`, adjusting the paths as needed:

```yaml
input: ../raw/uvis
output: ../redux
instrument: uvis
resume: true
options:
  background:
    mode: auto
  drizzle:
    final_scale: 0.04       # arcsec/pixel
    final_pixfrac: 0.8
    final_kernel: square
    final_wht_type: EXP
  astrometry:
    method: two_pass
    catalogs:
      cache: true
      offline: false
filter_overrides:
  F275W:
    background:
      mode: none          # Example: leave the sky in this filter
```

Run `hstpipe uvis reduce --config configs/uvis.yaml --resume`. The filter override
is an example, not a recommendation to disable all UV sky subtraction. Omit it
to use the conservative default background policy.

Common settings apply to all filters; `filter_overrides` changes only the named
filter. There is one target-level config and manifest, with independent filter
records and stage receipts. See the [complete configuration reference](uvis-full-field.md#complete-configuration-example)
for catalog quality, timeouts, relative alignment, and absolute-fit options.

Changing final scale or pixfrac can reuse upstream astrometry. Changing an
absolute reference invalidates downstream WCS propagation and final drizzle,
while leaving unchanged background work reusable. Removed cleanup diagnostics
may need regeneration. No receipt is rewritten to claim unperformed validation.

## Technical reference

| Topic | Read more |
| --- | --- |
| Source detection, both CCDs, and distortion | [Full-field astrometry](uvis-full-field.md) |
| Relative shifts, rotation, and uncertainty | [Relative-registration models](uvis-relative.md) |
| Gaia/Legacy APIs, epochs, and caches | [Direct reference catalogs](uvis-full-field.md#direct-reference-catalogs) |
| Compact-galaxy settings | [Compact reference selection](uvis-full-field.md#compact-references-and-source-selection-provenance) |
| Sparse absolute fits and final source QA | [Scientific acceptance](uvis-full-field.md#sparse-absolute-references-and-final-source-validation) |
| Background, DQ, cosmic rays, and calibration | [Handbook policies](uvis-handbook.md) |
| Previous tests and target investigations | [Historical records](VALIDATION.md) |

Historical results are not current target status summaries. Use each target's
report and manifest for its latest measured outcome.

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
