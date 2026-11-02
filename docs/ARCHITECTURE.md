# Architecture

The original 15-module STIS implementation was inspected and retained. Common
infrastructure handles YAML, paths, inventory, constant-array/multi-imset FITS
validation, provenance, logging and stage receipts. Instrument grouping remains
in each instrument namespace. No legacy import aliases are retained.

STIS retains its background, calibration, geometry, resampling, coadd, extraction,
optimal, products, io, reference and diagnostic modules. pipeline.py orchestrates
the six native science functions with content-addressed per-stage receipts.
No numerical scientific method changed. Intermediate path changes, new diagnostic
routing, registration.ecsv and shared manifests are infrastructure changes.
See WORKFLOWS.md for output layout and independent stage/API access.

UVIS separates instrument processing, astrometric inference and presentation:

```text
hstpipe/uvis/
    pipeline.py, config.py, storage.py     # orchestration and stage receipts
    inventory.py, grouping.py, chips.py   # actual SCI/ERR/DQ geometry
    background.py, amplifiers.py, dq.py   # detector/background policies
    drizzle.py, qa.py, source_qa.py       # STScI combination and science QA
    alignment.py                          # retained low-level/manual API
    astrometry/
        catalogs.py                       # direct Gaia/Legacy retrieval/cache
        full_field.py, sources.py         # native HST measurements
        centroids.py, stellar.py          # pixel, centroid and profile QC
        identification.py, common_sources.py, source_recovery.py
        relative_engine.py, relative_fit.py, target_fallback.py
        spatial_residuals.py              # amplitude-based residual diagnostics
        cross_filter.py                   # validated internal anchor transfer
        alignment.py, absolute_matching.py # absolute fitting/acceptance
        workflow.py, wcs.py                # TweakBack, final validation, HSTWCS
        orchestration.py                  # target/filter stage scheduling
        relative.py, matching.py          # WCS inspection/compatibility helpers
        pair_diagnostics.py, config.py
    reporting/
        __init__.py                       # generate_report public API
        report.py                         # offline HTML assembly/size limit
        decisions.py                      # concise status and recommended action
        astrometry.py                     # source/WCS provenance tables
        alignment.py                      # per-exposure reference-cross cutouts
```

The astrometry package is a scientific subsystem, not a second pipeline. External
catalog retrieval, relative HST-to-HST fitting and absolute fitting have different
inputs and acceptance criteria; keeping them together avoids scattering these
rules across the instrument root. The reporting package reads persisted stage
results and produces presentation only. Its imports must not run scientific fits.

Both full-frame and subarray metadata are supported. Native FLC centroids use
complete chip-specific HSTWCS distortion models. Standard relative registration,
local fallback, provisional AstroDrizzle, absolute Gaia/Legacy alignment,
TweakBack and final AstroDrizzle remain separate, validated stages. Conservative
background handling avoids double subtraction; corrections are on working copies.

Content-addressed receipts cover individual stages and independent filter groups.
Report-only code is excluded from scientific cache keys. No maintenance operation
rewrites a receipt to claim that changed scientific code has been validated.
Real UVIS reductions have passed the available checks, with sparse-source and
catalog-systematic limitations recorded per target/filter; historical run notes
are evidence, not current science-readiness status. STIS scientific processing is
unchanged.

Package imports are independent of the other instrument's optional dependencies.
Adding COS/ACS requires a new instrument package with its own grouping/scientific
engine and CLI dispatch, using common config/paths/run transactions. No STIS/UVIS
scientific internals need changing. No COS/ACS pipeline is implemented.

STScI wrappers use process environment/cwd; run independent reductions in separate
processes. Force/rebuild archives prior results and reopens new log files. Cache
receipts require unchanged input/content, effective settings, source, dependency
versions and native executable identities. Failed stages never get valid receipts.
