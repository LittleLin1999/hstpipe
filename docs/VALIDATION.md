# Standalone validation on erebus

2026-10-30: isolated Python 3.11 virtualenv, no inherited system site-packages.
All extras installed from declared dependencies. pip check passes. HSTCAL was
built afresh at the new software root from the pinned commit/patch, with CFITSIO
supplied by the base Conda prefix. No old observational software path is needed.

66 tests pass. Existing five low-sample background variance warnings remain.
Tests cover native STIS scientific behavior, independent masks/weights, YAML,
metadata, constant-array imsets, inventory, duplicates, cache invalidation,
partial resume, force/rebuild archiving, installed CLI, optional imports, UVIS
API dispatch, supplied sky-mask diagnostics, and subarray metadata.

Actual six-exposure J1022+0841/G140L reduction completed through both independent
2D branches and consolidated FLT/DRZ/X2D extraction. Both 2D SCI arrays equal the
prior validated outputs exactly. Boxcar fluxes equal the prior outputs; optimal
flux differences are below 6e-16 of peak (roundoff). DQ and NEXP match. Final FITS
checksums/WCS/units and original input identities pass. Only wavelength-edge
bins are masked: 980 finite FLT/X2D bins, 981 DRZ bins. The previous empirical
optimal-profile branch-flux discrepancy remains a scientific limitation; this
refactor does not claim to solve it.

Explicit stage-resume reused prepare, background, joint_drizzle, mast_coadd,
validate and extract with no numerical recomputation. Source/wheel distributions
include the support recipe, native pin/checksum, environment checker and lock.
No real UVIS astrometric or photometric reduction was performed.

All 3,883 files from the old observational HST tree were preserved by verified
hard links in the new archive. No active runtime symlink resolves into the old
tree. The old directory was not deleted automatically. Migration and scientific
validation reports reside beside the data, not in the installable package.
