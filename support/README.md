# Standalone installation support

Use Python 3.11 on Linux. Prerequisites: git, CMake >=3.11, a C compiler,
CFITSIO shared libraries/headers, and pkg-config. OpenMP is disabled in the
recipe. The upstream HSTCAL revision and patch checksum are pinned in
hstcal/provenance.json. No installed executable is patched silently.

From the repository:

    PYTHON=python3.11 CFITSIO_PREFIX=/path/to/cfitsio-prefix \
      bash support/bootstrap.sh /absolute/hst-env /absolute/hst-software
    export PATH=/absolute/hst-env/bin:/absolute/hst-software/hstcal-install/bin:$PATH
    python support/check_environment.py --hstcal /absolute/hst-software/hstcal-install/bin

CFITSIO_PREFIX is unnecessary when pkg-config locates cfitsio. On erebus CFITSIO
is supplied by the base Python Conda prefix, independently of OBSERVING/HST.
The virtualenv does not inherit system site-packages. Python dependencies,
including Drizzle 3, stistools, stsci.tools and DrizzlePac, are installed anew.

The environments/erebus-python311.lock.txt file records resolved Python versions
for the validated Linux/Python 3.11 environment; it is not a cross-platform lock.
Set HSTPIPE_USE_LOCK=1 for bootstrap.sh to install the lock before the editable
package and reproduce those versions.
Calibration reference data remain external, selected through reference_dir.
CRDS is optional when matching reference files are already present. No MAST
token is needed after observations are downloaded. No tokens or data are packaged.

Support is included in source distributions and wheels at share/hstpipe/support.
bootstrap.sh installs a source checkout. Wheel users install extras then run the
installed hstcal/build.sh separately. check_environment.py performs imports only.
