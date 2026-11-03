# Historical UVIS catalog and WCS validation

> Historical run record. Counts, paths, thresholds and target outcomes below
> describe those runs, not current defaults or current science-readiness status.
> Use [the UVIS guide](uvis.md) for current policy and each target manifest/report
> for its latest measured outcome. This evidence is retained for reproducibility.

Validation performed 2026-10-31 with Python 3.11.7, DrizzlePac 3.11.0,
STWCS 1.7.7, Astropy 8.0.2, Astroquery 0.4.11 and pyvo 1.9.1.

## Offline regression suite

126 passed, 1 opt-in real-data test skipped. Tests use synthetic FITS and mocked
network clients. Coverage includes native Gaia ADQL, motion/parallax/masked
parameters, Legacy schema/query/types/quality/releases, empty/malformed responses,
timeout cancellation/retries, offline caches and invalidation, WCS preservation,
reference fallback, cross-filter transfer, setup-safe retries, SCI extension
selection, TweakBack/headerlets, corrected-FLC final drizzle, output layout,
migration, unchanged STIS paths and raw immutability.

## Actual service queries

- NOIRLab TAP `ls_dr10.tractor`: southern cone (156.405, 14.029), radius 0.015 deg,
  returned 53 objects; 17 survived the configured quality selection. Coordinates
  and standardized RA/Dec units verified as degrees.
- Northern cone (160.6768, 37.364), radius 0.015 deg: 17 selected references,
  with release code 9011 correctly identified as DR9-North. The queried schema
  does not supply PHOTSYS; optional columns were selected from TAP_SCHEMA.
- Gaia DR3: real official Astroquery requests failed with network/read timeouts.
  A successful live Gaia response and end-to-end Gaia calibration have **not**
  been demonstrated on this host. ADQL, client dispatch, propagation, filtering
  and Gaia fitting are covered with mocked responses.

Logs: `/tmp/hstpipe-direct-catalog-smoke.log` and
`/tmp/hstpipe-direct-catalog-smoke-v2.log` on Erebus.

## Real four-FLC propagation and final drizzle

J1047+0739/F689M used exposures ifpg05l0q, ifpg05kyq, ifpg05kwq and ifpg05k2q.
An existing first-pass AstroDrizzle product was reused only after checking its
recorded input identities, relative FLC checkpoints, background/alignment/CR
parameters and provisional drizzle parameters. This avoided repeating an
unchanged first pass. External reference data came from the actual direct Legacy
query cache, with DR10-South identities preserved.

Native TweakReg retained 24 matches. STWCS archived the original provisional WCS;
native `tweakback.apply_tweak(provisional[SCI,1], ...)` updated both SCI chips in
all four fresh FLC copies. A new native AstroDrizzle pass performed separate
image drizzle, median, blot, CR rejection and final combination from these copies.
The final image contains PRIMARY/SCI/WHT/CTX/HDRTAB, SCI ELECTRONS/S, exposure-time
WHT with explicit WHTTYPE=EXP, and WCSNAME=HP-LEGACY.

Final validation retained 24 matches and rejected 3, with radial RMS 0.056765 arcsec,
robust scatter 0.033710 arcsec, median offsets (+0.003054, -0.002767) arcsec,
and uncertainty estimate 0.050471 arcsec including the configured 0.05 arcsec
catalog floor. These are fit-star statistics, **not independent holdout accuracy**.
The product is under `redux/uvis-validation/J1047+0739/uvis/products/`.

The image remains **review_required**: relative exposure registration was
unresolved (pairwise RMS 0.093–0.157 arcsec), one individual exposure/reference
check exceeded its configured 0.2 arcsec threshold, and CR diagnostics flag
bright-source regions. The final target-core weight ratio was 1.024; target-region
new CR fractions were approximately 0.39–0.85%. Successful absolute residuals do
not establish fully validated PSF/photometry or cure guiding/registration errors.

A separate actual cache-resume check applied all four saved STWCS headerlets to
fresh relative FLC copies and verified byte-identical SCI/ERR/DQ arrays and
unchanged upstream file identities. Raw FLC checksums also remained unchanged.
Logs: `/tmp/hstpipe-red-final2.log` and `/tmp/hstpipe-real-resume.log`.

J1025+1402/F467M supplied only four TweakReg matches against the configured
minimum of six; its absolute solution was rejected, its provisional image retained,
and no new science-ready DRC claimed. No fit criterion was loosened to force
success.

## Migration

Production and validation legacy filter directories were migrated to the shared
`{target}/uvis/` hierarchy with verified byte-identical FITS data and recoverable
`{target}/.uvis-legacy-*` trees. Per-target `migration.json` records mappings and
verification. Old relative-only products remain recoverable when the default
absolute workflow replaces them. STIS scientific source was not changed.

The standard `hstpipe uvis reduce` CLI was also executed on the production
J1047+0739/F689M group using the same default scientific criteria, with only
network timeout/retry controls shortened for the service smoke test. It performed
relative-fit acceptance checks, provisional drizzle, actual catalog query,
TweakReg, TweakBack and a new final drizzle. Candidate relative fits were rejected
by the configured limits; the original relative WCS was retained with explicit
review flags. The production DRC is at
`redux/J1047+0739/uvis/products/J1047+0739_F689M_drc.fits`.

WCS cache signatures exclude CRDER/CSYER fit-error metadata because STWCS archival
can retain those keywords without changing the coordinate transformation. Tests
verify that this does not hide an actual CRVAL change.

Actual stage-level resume passed after explicitly removing only the whole-run
receipt to exercise the stage caches. SCI/WHT/CTX/HDRTAB arrays were identical,
the native AstroDrizzle log and background receipt were unchanged, corrected
headerlets were restored from cache, and raw checksums remained unchanged.
Log: /tmp/hstpipe-standard-resume-fixed.log. Real cross-filter transfer and a
successful live Gaia calibration remain unverified.
