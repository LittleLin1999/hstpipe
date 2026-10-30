# Architecture and inspection

The source was scripts/stis_pipeline (15 scientific modules), a six-stage
workflow with JSON dataclasses, flat raw inputs, work/stage1..6, spectra/, and
three named products. There was no installable stispipe distribution in the
inspected source. Existing tests, the HSTCAL coordinate-export patch and its
license are retained. No observations, products, credentials or binaries belong
in this repository.

Science modules remain in hstpipe.stis: background, calibration, geometry,
resampling, coadd, extraction, optimal, products, io, reference, diagnostics.
workflow.py retains the six native stage functions; pipeline.py adapts external
configuration/paths. AST comparisons confirm unchanged algorithms after path
literal migration (STIS_REGRESSION.json). Existing synthetic scientific tests
are retained. No real-data reduction/regression was executed for this refactor.

common owns YAML configuration, paths, inventory, generic imset validation,
content identities, manifests, logging, and content-addressed cache receipts.
Instrument grouping and science interpretation remain in stis/uvis.

UVIS uses calibrated FLC copies and lazy imports of TweakReg/AstroDrizzle.
Alignment, drizzle, and QA are independent functions. Conservative optional sky
handling is delegated to AstroDrizzle's skysub/skymethod configuration. Custom
amplifier correction is deliberately not implemented without detector-specific
validation. No detector calibration or CTE correction is rerun.

Current cache integration is whole-reduction: inputs, engine sources, references,
controls, effective options, dependency versions and output hashes must match.
StageCache can be used directly for finer stages; automatic partial-stage resume
is not yet wired into the legacy STIS workflow. --force archives the old run.
Changing upstream data cannot silently reuse old results.

The preserved STIS intermediate filenames are migrated to intermediate/background,
intermediate/flt/registration, intermediate/flt/coadd, intermediate/x2d and
intermediate/x1d/flt. Legacy comparison tables remain in intermediate/x1d;
DRZ/X2D extraction numerical-model diagnostics are under diagnostics/extraction.
The full requested per-method directory split and registration.ecsv export are
not yet migrated. Native JSON metadata and formats are preserved.

Future ACS/COS modules should expose reduce_* and group_exposures, use common
run/config/paths, and own their detector models/calibration. Extend CLI dispatch;
do not modify existing STIS/UVIS science internals. Global cwd/environment in
STScI/legacy interfaces means reductions should run in separate processes.
