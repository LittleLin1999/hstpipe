# Standalone setup

Code repository: any checkout of hstpipe. Data and Python environment roots
are explicit and need not be inside the checkout. Native HSTCAL is installed
inside the Python environment; temporary build trees are discarded.

Prerequisites: Python 3.11 (validated Linux version), git, CMake >=3.11,
C compiler, pkg-config, CFITSIO shared library and headers. On Debian these are
provided by cmake, gcc, pkg-config and libcfitsio-dev. A Conda cfitsio prefix
can be supplied with CFITSIO_PREFIX. OpenMP is disabled for the pinned build.

    PYTHON=python3.11 CFITSIO_PREFIX=/native/prefix \
      bash support/bootstrap.sh /environment/hstpipe
    export PATH=/environment/hstpipe/bin:$PATH
    python support/check_environment.py --hstcal /environment/hstpipe/bin
    python -m pip check
    hstpipe stis reduce --config /observations/configs/stis.yaml

The helper installs all extras; minimal users may install .[stis] or .[uvis]
instead. HSTCAL cs6/cs7 are needed for STIS; the coordinate export is a pinned
local patch with checksum/provenance. It is checked against ordinary x2d output
before science rectification. The script never modifies observations.

The tested Linux/Python3.11 version lock is support/environments/erebus-python311.lock.txt.
It contains no hstpipe checkout path; install those dependencies then the desired
hstpipe checkout. A wheel includes support at share/hstpipe/support, including
native build instructions and patch; no C binaries or calibration data are
bundled. The source-distribution manifest includes the same support files.

Reference files stay external and their hashes are recorded in the run manifest.
CRDS can obtain appropriate references but is not required if these are already
present. The reduction tested on erebus uses only new hst/raw, hst/calibration,
hst/.venv paths. Its active installation does not depend on
OBSERVING/HST or on base Python site-packages.

HSTCAL executables and libraries now share the Python virtualenv prefix. Build
source and CMake trees are temporary and removed after installation. Setup logs
and validation records belong to each reduction, not a top-level administration
directory. The portable product checker is support/check_stis_products.py:

    python support/check_stis_products.py --run /redux/target/stis/grating --baseline /previous/reduction
