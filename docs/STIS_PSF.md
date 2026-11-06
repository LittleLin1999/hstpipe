# Empirical STIS G140L PSF calibration

[Home](../README.md) · [STIS reduction](WORKFLOWS.md#stis-reduction) · [Installation](INSTALLATION.md) · [UVIS PSFs](uvis-psf.md)

`hstpipe stis psf --config /path/to/workdir/config/psf.yaml --stage all` builds
calibration-star products using the existing STIS reduction stages. The working
directory contains data and configurations only. Optional dependencies are
available through the `psf` extra; the existing HST environment can be reused.

## Install and run

From the package checkout, install the STIS PSF extra into your science environment:

```sh
python -m pip install -e '.[psf]'
hstpipe stis psf --config /path/to/workdir/config/psf.yaml --stage all
```

The `psf` extra includes the STIS Python dependencies. The native reduction and
recalibration stages also require the [patched HSTCAL setup](INSTALLATION.md#manual-or-instrument-only-installation)
and appropriate calibration references.

Run one stage by replacing `all` with `discover`, `preview`, `download`,
`recalibrate`, `screen`, `reduce`, `extract`, `library`, `validate`, `report`,
or `plot`. These stages execute in that order for `all`. For example, to validate
an existing library and refresh its report:

```sh
hstpipe stis psf --config /path/to/workdir/config/psf.yaml --stage validate
hstpipe stis psf --config /path/to/workdir/config/psf.yaml --stage report
```

## Configuration and method

The PSF configuration uses the standard `input`, `output`, `instrument`, and
`options` schema. `options.psf` references selection and reduction YAML files.
Selection is metadata-first and preserves rejected candidates. `screen` performs
trace-aware PSF quality screening; `reduce` invokes the existing STIS workflow;
`extract` measures native and drizzled PSFs from those hstpipe outputs. Library,
validation, and report stages retain configuration-specific predictive uncertainty.

Individual group YAML files are ordinary reduction configurations and can be
reproduced with `hstpipe stis reduce --config ... --resume`. Archive products are
never changed. The regular and D1 apertures remain distinct. Repeller-wire or
boundary-limited PSFs use an explicitly recorded finite normalization radius;
unsupported aperture heights are not inferred. A library group containing one
star has no measured independent-star scatter. No universal slit-width correction
or uncalibrated wavelength extrapolation is provided.

Validation includes repeated visits, independent stars, wavelength, slit and
detector position, background models, native versus rectified/drizzled profiles,
leave-one-star-out null fits, and pixfrac/spatial-scale sensitivity using the
same hstpipe drizzle operation. Generic simulations are calibration diagnostics,
not science-target analysis. Inspect the resulting quality flags and reports
before using a PSF in another analysis.

## Plot the exported FITS library

```sh
hstpipe stis psf --config /path/to/workdir/config/psf.yaml --stage plot
```

This reads existing templates only. PNG figures compare slit widths versus
along-slit offset, with native/drizzled and wavelength panels. Only WING-grade
templates covering the common finite normalization radius are compared, without creating
a new stack. Predictive covariance is propagated through common-aperture
normalization. Selection identities and plotted densities are saved as ECSV.
Configure `plot_normalization_radius` (minimum template support radius; default 50 pixels),
`plot_comparison_radius` (exact shared integration radius; defaults to half a
native pixel inside the template radius), `plot_apertures`, and `plot_focus_band`
under `options.psf`. Epoch/position differences prevent a causal slit correction.

The same plot stage also reads individual native FITS profiles and plots shape,
FWHM, core concentration, wing/core ratios and asymmetry against the measured
source detector row, grouped by regular slit. The row is not a calibrated sky
offset from slit center. `plot_position_band`, `plot_position_core_height` and
`plot_position_wing_radius` configure these diagnostics. Conditional errors and
confounded source/epoch/position coverage are retained explicitly.

## Read-only acquisition and collapsed-profile diagnostics

`hstpipe stis acq --acquisition RAW... --science FLT... --output PLOT.png`
uses stistools acquisition target coordinates and science aperture headers to
plot a nominal target-anchored slit footprint. It does not measure achieved slit
centering or replace absolute astrometric registration. SCI,3 is a reference
lamp image, not the science target.

`hstpipe stis psf-collapse --library LIBRARY.fits --request REQUEST.json --output PROFILE.ecsv`
projects two explicitly selected compatible, bracketing wavelength PSFs onto
specified spatial bins with flux-conserving pixel overlap and trace shifts,
then collapses with explicit wavelength weights. No extrapolation is allowed.
Aperture mismatch requires an explicit proxy flag. The output propagates
predictive covariance with a fully correlated cross-band bound. This diagnostic
alone does not establish extension significance or a validated D1 transfer.
