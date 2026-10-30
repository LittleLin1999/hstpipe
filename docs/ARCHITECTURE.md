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

UVIS owns grouping, off-source background/amplifier-region assessment, working-copy
alignment, official AstroDrizzle integration, and SCI/WHT/CTX QA. Sky subtraction
is delegated to AstroDrizzle, and amplifier correction is not applied blindly.
Both full-frame and subarray metadata are supported; the STIS science engine
retains its current FUV-MAMA/G140L mode restrictions. UVIS has whole-run caching;
DrizzlePac intermediates are retained. Real UVIS reduction is not yet validated.

Package imports are independent of the other instrument's optional dependencies.
Adding COS/ACS requires a new instrument package with its own grouping/scientific
engine and CLI dispatch, using common config/paths/run transactions. No STIS/UVIS
scientific internals need changing. No COS/ACS pipeline is implemented.

STScI wrappers use process environment/cwd; run independent reductions in separate
processes. Force/rebuild archives prior results and reopens new log files. Cache
receipts require unchanged input/content, effective settings, source, dependency
versions and native executable identities. Failed stages never get valid receipts.
