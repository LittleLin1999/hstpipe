# UVIS regression fixes — 2026-11-06

The 20 failures recorded in `field-layout-validation.md` have been resolved.
The fixes address three production defects:

- Empty external associations caused `resolve_repeated_counterparts` to access
  nonexistent columns even when measured HST common sources were available.
  Empty associations now leave those measured sources unchanged, allowing the
  existing HST-only path to continue through its usual validation.
- The compatibility alignment adapter accepted `paths=None` but unconditionally
  called `paths.scratch`. It now uses a temporary candidate directory when no
  `ReductionPaths` is supplied. Both acceptance and rejection are tested with
  actual detected synthetic stars; rejected fits preserve the input bytes.
- Migration from the old nested per-filter layout minted current cache keys
  for legacy processing without current scientific validation. Migration now
  preserves the old products and recovery receipts, records
  `requires_revalidation`, and leaves current scientific receipts unset. A
  subsequent reduction executes validation rather than silently reusing them.

Each production defect was reproduced in a failing regression before its fix.
Other failures were stale fixtures: old Gaia client mocks, per-test leakage of
failed-query memoization, superseded centroid status names, old cache paths,
missing association fields, and expectations that ignored the current mask
policy or the runner's recorded per-filter errors. The refined-catalog test now
starts with too few sources and verifies that eight measured matches validate
an unchanged WCS without fitting a correction.

No science acceptance thresholds or DQ policy were relaxed. Expanded tests
reject hot/saturated pixels, calwf3 CR flags and fresh AstroDrizzle CR flags;
existing tests retain reassessment of inherited 4096 flags. Split-drizzle tests
cover all combinations of running CR rejection and possessing fresh masks.

## Validation

Executed in the existing Erebus HST environment with the candidate source:

```bash
PYTHONPATH=src MPLCONFIGDIR=/tmp/hstpipe-bugfix-mpl \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /home/linxiaojing/data/lrd-squared/hst/.venv/bin/python -m pytest \
  tests/uvis tests/common tests/stis -q --tb=short
```

Result: **429 passed, 0 failed, 1 skipped**, in 68.57 seconds.
The skipped test requires `HSTPIPE_UVIS_INTEGRATION_INPUT` with four real FLCs;
this run did not repeat an end-to-end native reduction of science observations.

The 19 currently published UVIS products were audited: all have nonempty
external-association tables and none was produced through the affected nested
legacy migration. The standard two-pass pipeline does not use the compatibility
alignment adapter unless manual initialization is requested. No scientific
products were regenerated for this patch; the existing science quality caveats
remain unchanged. In particular, J0129 F763M remains unpublished because its
relative astrometry was not validated; a software test fix does not override
that scientific rejection.
