# Installation support files

[User installation guide](../docs/INSTALLATION.md) · [Package README](../README.md)

This directory contains the reproducible native build recipe and installation
checks. Most users should start with `setup.sh` at the repository root.

## Requirements

Use Python 3.11 on Linux. Prerequisites: git, CMake >=3.11, a C compiler,
CFITSIO shared libraries/headers, and pkg-config. OpenMP is disabled in the
recipe. The upstream HSTCAL revision and patch checksum are pinned in
hstcal/provenance.json. No installed executable is patched silently.

## Run setup

The main entry point is `bash setup.sh /absolute/hst-env` (see the root README).
Use `bash setup.sh --check` to check prerequisites without creating an environment.
The native helper remains available directly. From the repository:

```bash
PYTHON=python3.11 CFITSIO_PREFIX=/path/to/cfitsio-prefix \
  bash support/bootstrap.sh /absolute/hst-env
export PATH=/absolute/hst-env/bin:$PATH
python support/check_environment.py --hstcal /absolute/hst-env/bin
```

CFITSIO_PREFIX is unnecessary when pkg-config locates cfitsio. On erebus CFITSIO
is supplied by the base Python Conda prefix, independently of OBSERVING/HST.
The virtualenv does not inherit system site-packages. Python dependencies,
including Drizzle 3, stistools, stsci.tools and DrizzlePac, are installed anew.

## Reproducibility and data

The `environments/erebus-python311.lock.txt` file records resolved Python versions
for the validated Linux/Python 3.11 environment; it is not a cross-platform lock.
Set HSTPIPE_USE_LOCK=1 for bootstrap.sh to install the lock before the editable
package and reproduce those versions.
Calibration reference data remain external, selected through reference_dir.
CRDS is optional when matching reference files are already present. No MAST
token is needed after observations are downloaded. No tokens or data are packaged.

## Installed files and checks

Support is included in source distributions and wheels at share/hstpipe/support.
bootstrap.sh installs a source checkout. Wheel users install extras then run the
installed hstcal/build.sh separately. `check_environment.py` checks imports and native executable presence; it does
not reduce observations or validate science outputs.

HSTCAL is installed directly into the virtualenv bin/lib directories. Build
trees use a temporary directory and are removed after installation. For a
read-only checksum/regression check use check_stis_products.py --run RUN
[--baseline PREVIOUS_RUN].
