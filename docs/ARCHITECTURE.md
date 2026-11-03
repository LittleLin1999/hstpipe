# Code architecture

[Home](../README.md) · [User workflows](WORKFLOWS.md)

This page is for developers. To run a reduction or review its output, start with
the [UVIS guide](uvis.md) or [STIS workflow](WORKFLOWS.md#stis-reduction).

## Package responsibilities

| Package | Responsibility |
| --- | --- |
| `hstpipe.common` | Configuration, paths, FITS inventory, logging, provenance, and stage receipts. |
| `hstpipe.stis` | STIS background, calibration, combination, and extraction. |
| `hstpipe.uvis` | UVIS detector processing and reduction orchestration. |
| `hstpipe.uvis.astrometry` | Source measurement, catalog retrieval, relative/absolute registration, and WCS propagation. |
| `hstpipe.uvis.reporting` | Offline HTML from saved results; no scientific fitting. |

## STIS

`pipeline.py` orchestrates six native stages: prepare, background, joint drizzle,
MAST coadd, validation, and extraction. The background, calibration, geometry,
resampling, coadd, extraction, optimal, products, I/O, reference, and diagnostic
modules retain the established scientific methods. See [stage access](WORKFLOWS.md#python-stage-apis).

## UVIS

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

`common/uvis_cleanup.py` implements explicit diagnostic pruning using shared
paths, identities, and locks. It is outside the scientific code cache inputs,
so installing the cleanup command does not invalidate reduction keys. Deleting
a diagnostic can still invalidate a receipt that requires that output.

## Caches and process boundaries

Content-addressed receipts cover individual stages and independent filter groups.
Report-only code is excluded from scientific cache keys. No maintenance operation
rewrites a receipt to claim that changed scientific code has been validated.
Historical run notes are evidence for particular datasets. The current target
manifest and report determine science readiness. STIS processing is unchanged.

Package imports are independent of the other instrument's optional dependencies.
Adding COS/ACS requires a new instrument package with its own grouping/scientific
engine and CLI dispatch, using common config/paths/run transactions. No STIS/UVIS
scientific internals need changing. No COS/ACS pipeline is implemented.

STScI wrappers use process environment/cwd; run independent reductions in separate
processes. Force/rebuild archives prior results and reopens new log files. Cache
receipts require unchanged input/content, effective settings, source, dependency
versions and native executable identities. Failed stages never get valid receipts.
