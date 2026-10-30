# hstpipe

Installable HST/STIS spectroscopy and WFC3 UVIS imaging code. Observations and
reductions live outside the installed package. Standalone installation and the six-exposure STIS reduction have been validated on erebus.

## Installation

Python >=3.10 in a dedicated STScI-compatible environment:

```sh
pip install -e .
pip install -e '.[stis]'
pip install -e '.[uvis]'
pip install -e '.[all,test]'
```

STIS also requires HSTCAL, matching reference tables, and official x2d controls
including the coordinate export used by the original engine. See support/hstcal
and docs/STIS_LEGACY.md. Configure `oref`/CRDS according to the STScI software
environment; do not substitute new reference versions unnoticed. UVIS requires
DrizzlePac, usable distortion references and calibrated PCTECORR=COMPLETE FLCs.
Imports do not require the other instrument's optional libraries.

## Usage

```sh
hstpipe stis inspect --input /data/raw/stis
hstpipe uvis inspect --input /data/raw/uvis
hstpipe stis reduce --config /data/stis.yaml --dry-run
hstpipe uvis reduce --input /data/raw/uvis --output /data/redux --dry-run
```

Remove --dry-run to process. --target, --grating/--filter, --setup select output
identities. Incompatible setups are refused rather than silently merged.
--resume reuses only a fully verified matching run; --force archives old results.
No processing is performed on installation or import.

```python
from hstpipe import stis, uvis
from hstpipe.common.config import ReductionConfig
from hstpipe.stis import reduce_stis
from hstpipe.uvis import reduce_uvis
plan = reduce_uvis(ReductionConfig('/data/raw/uvis', '/data/redux', 'uvis'), dry_run=True)
```

## Configuration

YAML paths are relative to the configuration file. Explicit CLI paths are
relative to the caller. API paths are explicit, independent of package location.

```yaml
input: /data/raw/stis
output: /data/redux
instrument: stis
target: J1022+0841
element: G140L
options:
  reference_dir: /data/calibration_refs
  control_dir: /data/official_controls
  extraction:
    height_native_pixels: 11
    optimal_height_native_pixels: 140
    profile_height_native_pixels: 151
    flux_mode: point
```

For UVIS set instrument: uvis and element: F606W; options accepts align (default
true), alignment (TweakReg keyword mapping), drizzle (AstroDrizzle keyword mapping).
Use drizzle.skysub: false to disable additional sky subtraction. DrizzlePac owns
cosmic-ray rejection and SCI/WHT/CTX construction. Review alignment and sky masks
for each dataset; this code-only refactor did not validate real UVIS photometry.

## External directories and products

```
raw/stis/                  # flat original MAST filenames
raw/uvis/                  # flat original MAST filenames
redux/{target}/stis/{grating}/[setup]/
  products/{target}_{grating}_{s2d_drz,x2d_coadd,x1d}.fits
  intermediate/{background,flt,x2d,x1d}/
  diagnostics/             # extraction models and legacy coadd/background QC
  logs/reduction.log
  config.yaml
  manifest.json
redux/{target}/uvis/{filter}/[setup]/
  products/{target}_{filter}_drc.fits
  intermediate/{alignment,drizzle}/
  diagnostics/{alignment,drizzle}/
  logs/reduction.log
  config.yaml
  manifest.json
```

Directories are created only when used. Instrument stages use shared outer paths;
STIS retains its internal formats and per-exposure names. The final X1D retains
FLT/DRZ/X2D tables, independent DQ_BOX/DQ_OPT, full-aperture normalized masked
optimal fitting and point/diffuse/morph conventions. See docs/WORKFLOWS.md for current directory and restart behavior, and
docs/INSTALLATION.md for reproducible native/Python installation.

Use `hstpipe.stis` for STIS imports and the `hstpipe stis` CLI. The old
`stispipe` and `stis_pipeline` compatibility namespaces and `stispipe` command
have been removed. Native Config.read still accepts old JSON for callers of
legacy workflow stages.

## Validation and extension

Run `pytest`. Synthetic STIS scientific tests, shared infrastructure, import/CLI,
UVIS API dispatch and subarray metadata tests do not execute observations.
The DrizzlePac adapters require real UVIS data validation before science use.
See docs/ARCHITECTURE.md and docs/STIS_REGRESSION.json for preserved science and
remaining limitations. COS and ACS are not implemented.

Current validation: 66 tests pass in an isolated environment. The two final STIS
SCI arrays match the prior validated reduction exactly; 1D flux differences are
below 6e-16 of peak. All six scientific stages were reused in an explicit
stage-resume test. No real UVIS reduction was executed.
