# Field layout and synchronization validation — 2026-11-06

New UVIS runs support explicit `field`/`--field` grouping and write to
`OUTPUT/uvis/FIELD/`. Original TARGNAME headers remain unchanged. Disconnected
footprints remain separate groups. Existing legacy reductions resume in place.
Report refreshes preserve scientific filter records and product bytes.

The candidate combines the local field/report layout changes, the cross-filter
registration work already shared with Erebus, and Erebus's existing STIS PSF
and acquisition tools. The STIS reduction directory layout is unchanged.

Executed in Erebus's existing HST Python environment:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m pytest tests/uvis tests/common tests/stis -q --tb=short
```

| Version | Passed | Failed | Skipped |
|---|---:|---:|---:|
| Pre-sync Erebus source | 385 | 22 | 1 |
| Synchronized candidate | 397 | 20 | 1 |

Every remaining failure is present in the baseline. Two older layout/cache
fixture expectations were corrected. New field/grouping, report compatibility,
publication-view protection, recovery paths, and concurrent execution checks
pass. The suite is not fully green: existing catalog-client, centroid/DQ-policy,
legacy migration/cache, and pipeline fixture failures remain. This validation
establishes no new test failures; it does not certify new scientific reductions.

## Follow-up bug fixes

The historical failures above were resolved in the subsequent
[UVIS bug-fix validation](bugfix-validation-20261106.md): 429 passed, no failures,
and one opt-in integration test skipped. The earlier table is retained as the
synchronization baseline.
