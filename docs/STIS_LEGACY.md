# STIS reduction scripts

`reduce_stis.py` is the single entry point. Reusable science modules live in
`hstpipe.stis`; all background fitting is in `background.py`.

## Run

From the project root, with the Python dependencies installed:

```bash
python scripts/reduce_stis.py config/stis_j1022.json
```

The supplied J1022 config points to data on **erebus**. Run where those paths
are accessible, or edit the config for your local data. For another target,
copy [the example config](../config/stis_config.example.json) and set its paths.
The output prefix is read from FITS `TARGNAME` and `OPT_ELEM`; no target redshift
or source-line windows are required.

## Workflow

1. **Prepare:** validate FLTs, official x2d controls, coordinates and calibration.
2. **Background:** estimate each exposure's residual background in count-rate space.
3. **Joint drizzle:** subtract, flux-calibrate and resample original FLTs once.
4. **MAST coadd:** subtract each MAST x2d background, combine, then assess residual sky.
5. **Validate:** generate QC plots, checksums and a provenance report.
6. **Extract 1D:** integrate the clean 2D spectrum and run official `x1d` on
   count-valued background-subtracted copies of each original FLT, then coadd.

Three final files are written under `<output_dir>/products/`:

- `<TARGNAME>_<OPT_ELEM>_s2d_drz.fits`
- `<TARGNAME>_<OPT_ELEM>_x2d_coadd.fits`
- `<TARGNAME>_<OPT_ELEM>_x1d.fits`

Original inputs are preserved. Existing final products are not overwritten;
choose a new `output_dir` for another run.

The 1D FITS has `FLT`, `DRZ`, and `X2D` BINTABLE extensions with separate
11-native-pixel boxcar and 140-native-pixel empirical-profile fluxes and uncertainties. The `spectra/` directory
retains per-method comparison tables and QC plots as intermediates:

- `<name>_from2d_x1d.fits`: finite-aperture flux from the clean 2D image.
- `<name>_individual_x1d_coadd.fits`: coadd of official per-exposure extractions.

The default boxcar height is 11 native pixels (0.2706 arcsec). Optimal extraction
uses 140 native pixels (3.444 arcsec), with its empirical template fitted over
151 native pixels (3.7146 arcsec). These are independent configurable heights.
Both apertures share a signed broadband-continuum core centroid, excluding
airglow and the strongest source feature; centering moves apertures only.
Set `extraction.center_method=header` to use the configured `center_arcsec`
without measuring a centroid. Uncertain automatic centers stop extraction. Set `extraction.enabled=false` for a 2D-only run. To add
1D spectra to existing 2D products without rerunning drizzle:

```bash
python scripts/reduce_stis.py config/stis_j1022.json --extract-only
```

The final `x1d.fits` uses `extraction.flux_mode=point` by default in all three
extensions. `FLUX_APER` appears only in an intermediate comparison table; it
intentionally retains the diffuse `x2d` convention so the two 2D extraction
methods can be checked on the same calibration basis. The diagnostic plot's
top panel compares final-product fluxes in the selected mode; its lower panel
shows the point-versus-diffuse calibration difference within the FLT branch.
Stage 2 now includes an optional **additional residual** check in `background.py`.
Primary fitting uses |y| >= 100 native pixels. Residual fitting uses |y|=84–90,
and disjoint validation uses |y|=92–98; |y| <= 80 remains protected. Geometry is
configurable through `background.source_radius`, `sky_inner`, and the residual
strip offsets. These regions must be checked for extended-source contamination.

The continuum candidate is a constant (order zero), assessed in broad 50-A
bins. Wavelength structure that a constant cannot describe is flagged, not
removed with an automatic flexible model. Each airglow interval is assessed
separately and uses unclipped means at native spectral sampling. Training must
exceed 5 screening standard errors, validation must exceed 3 with the same sign,
and both sides must agree. These thresholds are screening choices, not formal
Gaussian probabilities. No significant supported residual means no correction.

`BKGPRIMARY` and `BKGRESID` are stored separately in `work/stage2/*_rate.fits`;
`BKG` is their sum, consumed identically by drizzle and FLT extraction. `PRIMARYERR`
and `RESERR` contain correlated uncertainty bounds, added linearly in `BKGERR`.
`ERRBASE` is the unchanged pipeline statistical error. `RESDQ` flags unconstrained
or inconsistent residual regions without altering science coverage or weights.
The matching `*_residual.png` plots and `models.json` record every decision.

Official `x1d` now uses **BACKCORR=OMIT** on the total-background-subtracted count
copies, preventing duplicate subtraction. Its sensitivity, trace control and
point-source calibration remain unchanged. `background.residual_enabled=false`
disables the additional native correction; it does not enable CALSTIS clipping.
`work/stage3/primary_only.npz` retains the matched primary-only coadd for comparison.
Original data and the previous `reduction_v3` are preserved; the supplied config
writes a new `reduction_v4_residual` run.

See the [full documentation](../README.md) for dependencies, calibration-control
setup, configuration options, uncertainty limitations and tests.

## MAST x2d background treatment and extraction diagnostics

Stage 4 now fits a separate 2D background to **each original MAST x2d** before
combination, then checks the primary-cleaned coadd for an additional residual.
`mast_background` has the same configuration fields as `background`:
`enabled` controls the individual fits and `residual_enabled` controls the
coadd-level correction. All sky algorithms remain in `background.py`.

The MAST model is in the existing physical surface-brightness units. There is
no repeated calibration, count-rate conversion, dark correction, or pixel-area
scaling. Native-pixel-equivalent spatial masks preserve the same physical
source region. Observed airglow wavelengths are derived from the heliocentric
WCS using each exposure's V_HELIO. Background uncertainties acknowledge, but do
not reconstruct, pre-existing x2d interpolation covariance.

`work/stage4/*_x2d_bkgsub.fits` saves the individual cleaned images and BKG/BKGERR;
`primary_only_coadd.fits` preserves the coadd before the residual check. The
final `*_x2d_coadd.fits` includes BKGPRIMARY, PRIMARYERR, BKGRESID, RESERR,
BKG (total), BKGERR, ERRBASE, and RESDQ. `COADDSUB` records whether the additional
correction was accepted. Original MAST inputs remain unchanged. The consolidated
1D `X2D` extension now extracts this background-treated final product.

Stage 6 writes `spectra/diagnostics/{DRZ,X2D,FLT_<root>}_model.png`,
`*_profile.png`, `*_fit.fits`, and `*_fit.json`. These show the actual extraction
aperture, empirical spatial model, residual image, diagonal standardized
residuals, collapsed profiles, and chi-square per nominal degree of freedom.
The diagnostic FITS retains the data, model, residual, aperture/profile weights,
statistical and model errors, validity masks, and per-wavelength statistics.
`extraction_methods.png` compares boxcar/profile extraction with uncertainty
bands in every branch; `all_branches.png` compares all three branches.

These displays were informed by the local `unfold_jwst/extract1d.py` model/chi
panels and JWST's `save_profile`, `save_scene_model`, `save_residual_image`
outputs. The JWST fitting algorithm is **not** substituted for the STIS method.
Our profile extraction still uses equal-noise weights. Chi-square here is a
**descriptive diagonal diagnostic**, not a fit probability: the profile is
estimated from the same data, MAMA errors can be underestimated, and rectified
pixels are correlated. Pixels without usable diagnostic errors are counted and
masked only in the chi display, never excluded from the science extraction.

The supplied config writes `reduction_v5_x2d_background`, preserving earlier
runs. The MAST branch is no longer an un-subtracted-sky comparison; its saved
unsubtracted intermediate is used for the independent photometric check.

Profile-template training and automatic feature selection exclude wavelengths
outside PHOTTAB's calibrated support. This prevents blue-edge sensitivity
extrapolation from contaminating the profile. The diagnostic FITS also saves
`CHI_FLOOR`/`DIAGERR`: an alternate **display-only** scaling based on unclipped
variance in off-source |y|=100–150 native pixels over 32-column blocks. Raw
ERRBASE chi-square is retained. This floor changes neither the science ERR
arrays nor the extraction weights and cannot supply a formal fit probability.
The per-fit JSON reports both statistics and the profile centroid relative to
the resolved extraction center. `spectra/centering.json` records its offset
from the original header-based guess and wavelength-block uncertainty.
Blue dotted lines mark the boxcar; black dashed lines mark the optimal aperture.
`MODEL_FULL` and `PROFILE_FULL` show the wider trained template.
Point-mode boxcar and optimal corrections are measured with separate official
unit-SCI controls at their respective heights. Morph mode requires separate
`APFRAC_BOX` and `APFRAC_OPT` model columns, including slit losses.
The wide profile is finite and wavelength-averaged, not an infinite instrumental
PSF. Optimal extraction can miss emission with a different spatial morphology;
inspect the line-profile residuals and retain the 2D data for extended-emission work.

Point mode requires an exact PCTAB extraction-height entry. HSTCAL otherwise
chooses the nearest height (61 would use 21 for these references), not a height
interpolation. The 140-pixel default has an exact entry. Each exposure records
the chosen row height and reference checksum in its extraction provenance.
The 151-pixel template-fitting region remains inside the protected source mask.


## Masked-pixel extraction and wavelength gaps

Boxcar and optimal extraction have independent `DQ_BOX`/`DQ_OPT`,
`EXPTIME_BOX`/`EXPTIME_OPT`, and `NEXP_BOX`/`NEXP_OPT` columns. Historical
`DQ`, `EXPTIME`, and `NEXP` columns are aliases for BOX only. Select good optimal
bins with `DQ_OPT == 0`, not `DQ == 0`. The plots use the matching method masks.

Following the normalized-profile fit in local `unfold_jwst.extract_optimal`,
optimal flux is fitted using valid spatial pixels while the template remains
normalized over the full aperture. Missing pixels do not redefine the requested
flux. Default acceptance requires `PROFCOV_OPT >= 0.8` (observed model light) and
`INFOFRAC_OPT >= 0.5` (retained squared-profile weight); these thresholds are
configurable as `optimal_min_profile_fraction` and
`optimal_min_information_fraction`. `OPT_MASKED` records accepted use of a
partially masked profile. The coadd support fractions are exposure-weighted
averages of the accepted inputs, not sums of coverage from different exposures.
The fitted missing light is conditional on the empirical model. Errors propagate
the actual linear coefficients; uncertain profile shape is not included.

The reference routine uses inverse variances and a 0.4 profile-fraction threshold.
We retain equal-noise MAMA weights and more conservative coverage thresholds;
we do not adopt its negative-pixel clipping or Gaussian fallback. Its variance
expression includes a surviving-profile factor; here variance is propagated
explicitly with squared coefficients for the full-aperture estimator.

A direct boxcar still requires complete coverage of its own 11-pixel aperture:
missing boxcar pixels are neither profile-filled nor silently renormalized.
Each method is coadded independently. Point-response transfer uses calibrated
wavelength support and valid unit-SCI sensitivity, independent of science DQ.
No wavelength gaps are interpolated. Unsupported calibration/coverage stays NaN.
See `spectra/extraction_coverage.png` and `coverage_summary.json` for diagnostics.
