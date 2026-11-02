# UVIS handbook review and operational policy

Reviewed the WFC3 Data Handbook v6.0 (June 2024), DrizzlePac Handbook v3.0
(January 2025), the current STScI amplifier-equalization notebook and the
installed DrizzlePac 3.11/STWCS APIs. This extends the existing package; STIS
scientific processing is unchanged.

## Already satisfied

Inventory discovers every SCI/ERR/DQ imset and its physical CCDCHIP metadata,
including single-SCI subarrays. Native centroid measurements from both CCDs
are projected with complete HSTWCS models into a common tangent plane. A single
exposure-level correction updates all SCI extensions through native STScI
interfaces. Native AstroDrizzle receives whole FLC filenames without a group=1
restriction. Gaia/Legacy identification priors do not supply relative-fit
coordinate constraints.

Relative WCS candidates, alternates and headerlets are tested before the normal
provisional drizzle. Target-local recovery follows failed standard registration.
Absolute fitting operates on accepted provisional images, TweakBack propagates
the correction to individual working FLCs, and final AstroDrizzle uses those
FLCs. Unvalidated groups do not authorize propagation or a new science-ready DRC.

Native integration checks on small synthetic arrays carrying actual full-frame
and subarray calibrated WCS/lookup HDUs verified both detector footprints in WHT,
two HDRTAB rows for a two-chip input and one for a subarray. Actual TweakBack
updated every available SCI, retaining distortion, with propagation errors below
0.00001 output pixels. These are API/geometry checks, not external astrometric
validation of the original observations.

The conservative background model uses source masks, significance/blank-sky
checks and explicit zero-subtraction fallback. The existing amplifier correction
uses detector geometry and PFLT information, remains disabled by default and
requires measured significance. No extra CTE calibration, chip flux normalization
or final-image PAM correction is applied. SCI remains electrons/s; WHT semantics
and correlated-noise limitations are documented without a fabricated ERR array.

## Changes required by the review

`dq.py` defines official UVIS flags and separate detector, centroid, background
and science policies. Bit 4096 identifies AstroDrizzle rejection; 8192 identifies
calwf3 combination rejection. Astrometric detection reassesses inherited 4096
flags, retaining the original flags as diagnostic provenance. Saturation and
calwf3 rejection remain excluded from precision centroids. New in-run CR masks
remain bad for subsequent source checks. Only selected CR bits are grown; the
old blanket dilation of every nonzero DQ pixel was removed. Shared masked-core,
halo, morphology and centroid-stability checks still apply to every detector.

The installed MDRIZTAB `2ck18260i_mdz.fits` was inspected. Its four-exposure
recommendation includes 336 (16+64+256), pixfrac 1.0 and unchanged output scale.
The pipeline records that recommendation without silently importing every table
parameter. Its default multi-exposure science good-bit mask is 80 (16+64),
conservatively retaining saturation as bad. A single exposure and diagnostic
coadds use stricter masks. Explicit science-mask configuration remains recorded;
it cannot make saturated cores valid astrometric centroids.

New CR passes use AstroDrizzle's resetbits=4096 on working copies. A split final
combination carrying freshly generated CR masks uses resetbits=0 so it cannot
erase those masks. Native reset_dq_bits is also used when restoring cached masks.
Only 4096 may be reassessed; all other original detector/calwf3 bits are checked
for preservation. CR summaries separate inherited, cleared, reconfirmed, newly
identified and current rejection, including source-core diagnostics.

DrizzlePac 3.11's native reset uses an unsigned-array/negative-integer operation
that raises OverflowError with installed NumPy 2. Working DQ arrays are therefore
represented as signed int32 where needed. Every 16-bit UVIS flag is represented
exactly; native resetbits and CR algorithms remain unchanged. Original DQ arrays
are never edited. Focused native tests verify this compatibility treatment on
both chips, including preservation of calwf3 and saturation flags.

Auto sky subtraction is disabled for narrowband and explicitly field-filling
emission. Low-background UV retains conservative reliability/significance checks.
HPSKY records the amount actually removed; HPSKYBG and HPSKYOK retain the measured
level and reliability even when subtraction is disabled. skysub=False with
skyuser=HPSKY follows the installed API: already-removed sky is accounted for
without another subtraction, and MDRIZSKY records it. Sky left in SCI also remains
in the reference used for CR noise. A competing skyfile subtraction is not applied.
Per-chip/exposure residual levels are compared without forcing a complex fit.

WCS inspection now retains available fit/catalog/distortion metadata. FLUXCORR
and photometric normalization keywords are checked and recorded, never reapplied.
Calibrated sky alone is not interpreted as the original total CTE-mitigating
background: postflash and dark have already been removed. Low-background/CTE,
guiding and conditional detector-artifact concerns remain explicit warnings.

An additional native catalog handoff defect was found during the real run:
TweakReg regenerated an automatic catalog for its separate reference image even
though input-image catalogs had passed pipeline QC. The pipeline now supplies
that reference explicitly from independently measured HST centroids, transformed
through every contributing chip's current HSTWCS. No external coordinates enter
this relative reference. A real two-chip TweakReg test disables automatic source
finding and verifies the resulting shared correction, preventing both redundant
reference detection and a reference-side centroid-QC bypass.

The all-target rerun exposed another actual detection defect in exploratory
no-CR mosaics: an independent-detector-pixel smoothing-noise approximation
joined millions of correlated/background-contaminated mosaic pixels into
giant segments, making deblending extremely slow. Mosaic segmentation now
measures the filtered image's own background and noise and uses the more
conservative threshold. Native FLC detection and all SCI pixels remain unchanged.
The threshold and measured filtered statistics are saved with the catalog.
A focused synthetic no-CR regression retains a real compact source, rejects the
noise-dominated segmentation behavior and verifies input immutability. A changed
detector dependency hash is never silently certified as an unchanged cache.

The final product is staged until absolute and source-core/WHT review gates pass;
questionable new results remain intermediate products while existing production
products are preserved.

## Sampling and optional measures

driz_sep_pixfrac remains 1.0. The configured final defaults are 0.04 arcsec/pixel,
a square kernel and pixfrac 0.8. Final WHT/source-core checks must verify coverage
after shrinking the drops. Commanded offsets,
effective sampling settings, WHT coverage and correlation diagnostics are saved.
An approximately 0.0333 arcsec/pixel grid is configurable after evaluating actual
dither sampling and science needs; neither smaller pixels nor smaller pixfrac is
assumed superior. CR thresholds are not relaxed to compensate for registration
or DQ-mask failures.

Amplifier correction, special masks for ghosts/scattered light/crosstalk and
saturated-star photometry remain observationally dependent. No mandatory rare-
artifact stage or saturated-core interpolation was added. Saturated PSF-wing
centroid recovery and advanced uncertainty/covariance propagation remain outside
the current implementation.

## Shallow files, restart and reports

```text
redux/{target}/uvis/
  products/          # target_filter_drc.fits; target_uvis_report.html
  intermediate/
    catalogs/        # shared Gaia/Legacy raw/normalized tables and flat receipts
    background/      # root_sky.fits, root_mask.fits, optional amplifier models
    alignment/       # root_sciN_sources.ecsv, matches/WCS/headerlets/flat receipts
    drizzle/         # provisional and essential diagnostic/reusable products
    working/         # modified FLC checkpoints
  diagnostics/       # flat, compressed report previews selected from these files
  logs/              # flat
  config.yaml
  manifest.json
```

Directories are created when files are written. Native-tool scratch workspaces
are outside this tree. Filter/setup/root/SCI identifiers distinguish filenames.
Recoverable migration uses a sibling backup and translates artifact references
and verified receipts. Layout-only changes can rebase compatible checkpoints;
they cannot certify a changed DQ/centroid/sky policy. Changed scientific policies
invalidate their dependent stages. Shared catalog coverage, query/epoch and raw
receipt validation remain mandatory; no catalog download is stored in raw/uvis.

The standalone HTML report is generated from the manifest on successful or failed
target processing. It includes per-chip counts, effective DQ policy, recovered
bright-source counts, relative/absolute residuals, sky/CR/WHT warnings and product
status. Compressed embedded figures have a fixed byte budget; the complete UTF-8
HTML is measured and must be strictly below 1,000,000 bytes. Existing products
are distinguished from products generated in the current attempt.

## Sources

- [WFC3 Data Handbook v6.0](https://hst-docs.stsci.edu/wfc3dhb)
- [UVIS calibration and DQ definitions](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-2-uvis-data-calibration-steps)
- [DrizzlePac Handbook v3.0](https://hst-docs.stsci.edu/drizzpac)
- [AstroDrizzle sky, DQ, CR and sampling guidance](https://hst-docs.stsci.edu/drizzpac/chapter-6-reprocessing-with-the-drizzlepac-package/6-3-running-astrodrizzle)
- [STScI UVIS amplifier equalization notebook](https://spacetelescope.github.io/hst_notebooks/notebooks/WFC3/uvis_amp_equalization/uvis_amp_equalization.html)
- [Pixel-based CTE correction limitations](https://hst-docs.stsci.edu/wfc3dhb/chapter-6-wfc3-uvis-charge-transfer-efficiency-cte/6-4-the-pixel-based-model)
