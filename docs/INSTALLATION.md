# Installation

[Home](../README.md) · [UVIS guide](uvis.md) · [STIS workflow](WORKFLOWS.md#stis-reduction)

## Automated Linux setup

From a source checkout, run:

```bash
bash setup.sh "$HOME/.venvs/hstpipe"
source "$HOME/.venvs/hstpipe/bin/activate"
```

The script checks prerequisites, creates an isolated virtualenv, installs
`requirements.txt`, builds the pinned and patched HSTCAL into the environment's
`bin`/`lib` directories, then checks imports, native executable presence, the CLI,
and dependency consistency. Build trees are temporary. No data are downloaded or
reduced. An existing nonempty environment is never overwritten.

Run `bash setup.sh --check` for a read-only prerequisite check, or
`bash setup.sh --help` for options. `support/bootstrap.sh` remains an equivalent
entry point. Both can be called from outside the checkout.

## System prerequisites

The complete setup is validated on Linux with Python 3.11. Python package metadata
allows Python >=3.10. The setup prefers `python3.11`, falling back to `python3`;
set `PYTHON=/path/to/python` to choose explicitly.

Install these before running setup:

- Python with `venv` and `ensurepip` support (3.11 recommended).
- git, CMake >=3.11, a C compiler, and pkg-config.
- CFITSIO headers and shared library.
- Standard Linux `realpath` and `sha256sum` commands.

On Debian/Ubuntu, the development tools are available as:

```bash
sudo apt-get install git cmake gcc pkg-config libcfitsio-dev python3-venv
```

This uses your distribution's Python version; install Python 3.11 separately if
needed. The setup script does not install operating-system packages or invoke sudo.

If CFITSIO is not found by `pkg-config`, supply its installation prefix:

```bash
PYTHON=python3.11 CFITSIO_PREFIX=/path/to/cfitsio bash setup.sh /path/to/new-env
```

The prefix must contain `include/fitsio.h` and `lib/libcfitsio.so*` (or `lib64`).
On erebus, the base Conda prefix supplies CFITSIO. This is a native-library
prerequisite; the new virtualenv does not inherit base Python site-packages.

## Reproduce the tested dependency versions

```bash
HSTPIPE_USE_LOCK=1 PYTHON=python3.11 bash setup.sh /path/to/new-env
```

`support/environments/erebus-python311.lock.txt` records the validated Linux/Python
3.11 versions. It is not a cross-platform lock. Without this option, pip resolves
the dependency ranges declared in `pyproject.toml`. The generated environment's
installed versions are saved in `share/hstpipe/native/python-environment.lock.txt`.

## Manual or instrument-only installation

If you already manage your own Python environment, install from the checkout:

```bash
python -m pip install -r requirements.txt   # Both instruments and tests
# Or choose only the dependencies you need:
python -m pip install -e '.[stis]'
python -m pip install -e '.[uvis]'
```

`requirements.txt` refers to the package extras in `pyproject.toml`; it does not
build HSTCAL or provide calibration data. These Python-only commands can be tried
on other platforms, but a complete native macOS/Windows setup is not validated.

For STIS, build the required HSTCAL into your environment and add its `bin` to PATH:

```bash
bash support/hstcal/build.sh /absolute/path/to/environment
export PATH=/absolute/path/to/environment/bin:$PATH
```

The HSTCAL revision and coordinate-export patch checksum are pinned in
`support/hstcal/provenance.json`. STIS needs `cs6.e` and `cs7.e` from that build;
an unpatched stock HSTCAL does not provide the required coordinate export.
OpenMP is disabled in the supplied build recipe.

Wheels include native build support at `share/hstpipe/support` under the
installation prefix. Wheel users install the desired extras, then run that
installed `hstcal/build.sh` separately. `setup.sh` and `requirements.txt` are
source-checkout entry points and are included in source distributions.

## Verify and troubleshoot

After full setup:

```bash
python support/check_environment.py
python -m pip check
hstpipe --help
```

For code changes, run the focused tests covering the modified behavior. A full
science reduction or broad test suite is not an installation prerequisite.

The environment checker imports both instruments' dependencies and checks that
HSTCAL executables exist. It does not reduce data or validate native science
outputs. Use `--hstcal /path/to/bin` if the executables are outside `sys.prefix/bin`.

| Problem | Action |
| --- | --- |
| Destination is not empty | Activate that environment or choose a new path; setup does not overwrite it. |
| Python version or `venv`/`ensurepip` error | Select Python 3.11 with `PYTHON` and install its venv support. |
| CFITSIO not found | Install its development package or set `CFITSIO_PREFIX`. |
| Missing native build tools | Install the system prerequisites above. |
| `hstpipe` command not found | Activate the environment used for installation. |
| STIS coordinate-export error | Use the supplied pinned HSTCAL build, rather than a stock executable earlier on PATH. |
| Missing calibration reference | Supply the matching external tables and configure `reference_dir` (STIS) or `iref` (UVIS). |

Calibration references and observations remain external. CRDS is optional when
matching references are already present. A MAST token is not needed to reduce
already-downloaded observations. For a product/regression check, use:

```bash
python support/check_stis_products.py --run /redux/target/stis/grating
# Add --baseline /previous/reduction to compare with an available prior run.
```
