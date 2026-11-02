# hstpipe

Reduce **HST/STIS spectra** and **WFC3/UVIS images** from calibrated FITS exposures.
`hstpipe` combines exposures, writes science products, and keeps the intermediate
files, diagnostics, configuration, and provenance needed to inspect each run.

| Workflow | Inputs | Main outputs | Status |
| --- | --- | --- | --- |
| STIS FUV-MAMA/G140L | Calibrated `*_flt.fits`, matching `*_x2d.fits`, and reference tables | Two combined 2D spectra and boxcar/optimal 1D spectra | Validated on a six-exposure dataset |
| WFC3/UVIS | CTE-corrected `*_flc.fits` with `PCTECORR=COMPLETE`, plus distortion references | Drizzled SCI/WHT/CTX image | Two-pass Gaia/Legacy astrometry and persistent-source image QA exercised on J104242; per-target validation required |

Other STIS modes, COS, and ACS are not currently supported. Data downloading is
separate: install the software first, then supply your observations and references.

## 1. Install with one setup command

The automated setup targets **Linux**. Use **Python 3.11**, the version used for
validation. You also need git, CMake >=3.11, a C compiler, pkg-config, and CFITSIO
headers/shared libraries; see [system prerequisites](docs/INSTALLATION.md).

```bash
git clone https://github.com/LittleLin1999/hstpipe.git
cd hstpipe
bash setup.sh "$HOME/.venvs/hstpipe"
source "$HOME/.venvs/hstpipe/bin/activate"
```

`setup.sh` creates a new isolated environment, installs both instrument dependencies
and the package, builds the pinned HSTCAL with the required coordinate-export
patch, and checks the installation. It leaves observations untouched. Internet
access is required to download dependencies and HSTCAL source.

Useful setup options:

```bash
# Check prerequisites without installing anything.
bash setup.sh --check

# Choose a Python interpreter or an existing CFITSIO installation.
PYTHON=python3.11 CFITSIO_PREFIX=/path/to/cfitsio bash setup.sh /path/to/new-env

# Reproduce the tested Linux/Python 3.11 dependency versions.
HSTPIPE_USE_LOCK=1 bash setup.sh /path/to/new-env
```

Omit the environment path to create `.venv` in the checkout. The script refuses to
modify a nonempty destination; activate an existing environment instead.

For Python-only installation, `requirements.txt` installs both instrument extras
and tests from the current checkout:

```bash
# Run inside an environment you have already created, from the repository root.
python -m pip install -r requirements.txt
```

This last command does **not** build HSTCAL. Dependency declarations live in
[pyproject.toml](pyproject.toml); the tested version lock and manual/instrument-only
installation instructions are in [the installation guide](docs/INSTALLATION.md).
The full native setup has not been validated on macOS or Windows.

## 2. Keep observations separate from code

A working directory can look like this:

```text
observations/
├── raw/
│   ├── stis/              # Original *_flt.fits and matching *_x2d.fits
│   └── uvis/              # Original *_flc.fits
├── calibration/
│   ├── stis/              # Reference tables matching the STIS exposures
│   └── uvis/              # UVIS distortion and other reference files
├── configs/
│   └── stis.yaml
└── redux/                 # Created when a reduction runs
```

Keep original MAST filenames and place inputs directly in each instrument's raw
directory. Keep calibration files in `calibration`, rather than mixing them with
science exposures. No observations or reference data are bundled with the package.

## 3. Run a STIS reduction

First inspect the exposures. The command prints a JSON inventory and compatible
target/grating groups:

```bash
hstpipe stis inspect --input /path/to/observations/raw/stis
```

Save the following as `observations/configs/stis.yaml`. Replace `J1022+0841` with
the target name in your FITS headers:

```yaml
input: ../raw/stis
output: ../redux
instrument: stis
target: J1022+0841
element: G140L
options:
  reference_dir: ../calibration/stis
  prepare_controls: true
  extraction:
    height_native_pixels: 11
    optimal_height_native_pixels: 140
    profile_height_native_pixels: 151
    flux_mode: point
```

`prepare_controls: true` generates the official x2d controls using the patched
HSTCAL installed by the setup script. To reuse existing controls, set
`options.control_dir` to their directory and `options.prepare_controls: false`.
That directory must contain `*_flt_omit.fits`, `*_flt_perform.fits`, and
`*_flt_mapping.f64` for each input FLT. Use reference tables that match the
observations; changing calibration versions changes the reduction.

Preview the selected inputs and output location, then start processing:

```bash
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml --dry-run
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml
```

`--dry-run` checks input selection and prints the plan; it does not execute science
stages or establish that all calibration dependencies are ready.

To continue a run or deliberately start it again:

```bash
# Reuse verified outputs; rerun stages whose inputs or dependencies changed.
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml --resume

# Preserve the existing run in a timestamped backup and rebuild.
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml --force
```

YAML paths are relative to the YAML file. Explicit command-line paths are relative
to your current directory. Separate incompatible instrument setups into separate
input selections; `--setup` adds an output label and does not filter exposures.

## 4. Find the results

For the example above, outputs are written under
`redux/J1022+0841/stis/G140L/`:

| Location | Contents |
| --- | --- |
| `products/*_s2d_drz.fits` | Combined spectrum from the FLT/drizzle branch |
| `products/*_x2d_coadd.fits` | Combined spectrum from the x2d branch |
| `products/*_x1d.fits` | FLT/DRZ/X2D tables containing boxcar and optimal extractions |
| `diagnostics/` | Background, combination, extraction, and flux-calibration checks |
| `intermediate/` | Per-exposure models and intermediate products |
| `logs/` | Reduction and stage logs |
| `config.yaml`, `manifest.json` | Effective configuration, input identities, and provenance |

Inspect diagnostics before interpreting the spectra. For optimal extractions,
use `DQ_OPT`; `DQ` is the boxcar mask (`DQ_BOX`). See
[workflow details](docs/WORKFLOWS.md) for the complete product layout.

## WFC3/UVIS workflow

UVIS uses DrizzlePac for alignment, cosmic-ray rejection, sky handling, and
AstroDrizzle combination. Set `iref` to the directory containing the distortion
references required by your exposures:

```bash
export iref=/path/to/observations/calibration/uvis/
hstpipe uvis inspect --input /path/to/observations/raw/uvis
hstpipe uvis reduce --input /path/to/observations/raw/uvis \
  --output /path/to/observations/redux --target YOUR_TARGET --filter F606W --dry-run
```

Remove `--dry-run` to process. Compatible groups share the target directory;
filters and any necessary setup identifiers appear in filenames:
`redux/{target}/uvis/products/{target}_{filter}[_setup]_drc.fits`.

The standard workflow validates relative registration across all available SCI
extensions, makes a provisional mosaic, obtains validated Gaia/Legacy absolute
astrometry, propagates it with TweakBack, and drizzles the corrected FLC working
copies. Unvalidated stages block publication; original FLCs remain unchanged.
An offline, self-contained report is generated even on failure and is smaller
than 1,000,000 bytes.

Configure `options.astrometry` and `options.drizzle`; use
`options.background.mode: none` when no sky subtraction is wanted. The pipeline
handles AstroDrizzle sky keywords to avoid double subtraction. Final defaults are
0.04 arcsec/pixel and pixfrac 0.8; provisional/CR reference sampling is separate.
See [UVIS configuration and scientific limitations](docs/uvis.md), including
compact-galaxy references and scoped restart behavior.

## Use from Python

The Python API accepts the same YAML configuration as the CLI:

```python
from hstpipe.stis import reduce_stis

plan = reduce_stis("/path/to/observations/configs/stis.yaml", dry_run=True)
results = reduce_stis("/path/to/observations/configs/stis.yaml")
```

`hstpipe.uvis.reduce_uvis` provides the equivalent imaging API. Both also accept a
`ReductionConfig` or dictionary. See [workflows](docs/WORKFLOWS.md) for stage APIs.

## Checks and further documentation

Use focused tests for the behavior being changed; ordinary unit tests require no
network. Historical STIS validation found exact agreement of the two 2D SCI
arrays with the prior validated reduction and 1D flux differences below
`6e-16` of peak. This covers the tested dataset, not all possible observations.
The empirical optimal-profile branch-flux discrepancy remains a known limitation.

- [Installation and troubleshooting](docs/INSTALLATION.md)
- [Workflows, restart behavior, and output layout](docs/WORKFLOWS.md)
- [Validation results and limitations](docs/VALIDATION.md)
- [Historical STIS science-method notes](docs/STIS_LEGACY.md)
- [Package architecture](docs/ARCHITECTURE.md)

UVIS imaging reduction, configuration and scientific QA: [docs/uvis.md](docs/uvis.md).
