#!/usr/bin/env bash
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
usage() {
    cat <<'HELP'
Usage: bash setup.sh [--check] [ENV_DIR]

Create a new isolated environment for STIS and UVIS, build the pinned HSTCAL,
and verify the installation. ENV_DIR defaults to .venv in this checkout.
Run --check to check prerequisites without creating an environment.

Requires Linux, Python >=3.10 (3.11 recommended), git, CMake >=3.11,
a C compiler, pkg-config, and CFITSIO headers/shared library.
The installation needs internet access. Observations are not downloaded.

Environment variables:
  PYTHON=python3.11       Interpreter to use (preferred default: python3.11)
  CFITSIO_PREFIX=/prefix CFITSIO prefix, if not found by pkg-config
  HSTPIPE_USE_LOCK=1      Use tested Linux/Python 3.11 dependency versions

Examples:
  bash setup.sh "$HOME/.venvs/hstpipe"
  bash setup.sh --check
  HSTPIPE_USE_LOCK=1 bash setup.sh /data/hst/.venv
HELP
}
fail() { printf 'Setup: %s\n' "$*" >&2; exit 1; }
check_only=0
environment_path=''
while (($#)); do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --check) check_only=1 ;;
        -*) fail "Unknown option: $1 (see --help)" ;;
        *) [[ -z "$environment_path" ]] || fail 'Supply only one environment directory.'
           environment_path=$1 ;;
    esac
    shift
done
[[ $(uname -s) == Linux ]] || fail 'The native setup recipe is validated on Linux. See docs/INSTALLATION.md for manual Python installation on other platforms.'
if [[ -n "${PYTHON:-}" ]]; then
    python=$PYTHON
elif command -v python3.11 >/dev/null 2>&1; then
    python=python3.11
else
    python=python3
fi
command -v "$python" >/dev/null 2>&1 || fail "Python not found: $python. Install Python 3.11 or set PYTHON."
"$python" - <<'CHECK'
import sys
if sys.version_info < (3, 10):
    raise SystemExit('Setup requires Python >=3.10; Python 3.11 is recommended.')
import venv, ensurepip
CHECK
if [[ ${HSTPIPE_USE_LOCK:-0} == 1 ]]; then
    "$python" -c 'import sys; assert sys.version_info[:2] == (3,11), "The supplied version lock requires Python 3.11."'
fi
venv=$("$python" -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "${environment_path:-$repo/.venv}")
[[ ! -e "$venv" ]] || [[ -d "$venv" && -z "$(ls -A "$venv")" ]] || fail "Destination is not empty: $venv. Activate the existing environment or choose a new path."
for tool in git cmake pkg-config realpath sha256sum "${CC:-cc}"; do
    command -v "$tool" >/dev/null 2>&1 || fail "Missing prerequisite: $tool. See docs/INSTALLATION.md."
done
"$python" - <<'CHECK'
import re, subprocess
output = subprocess.check_output(['cmake', '--version'], text=True)
version = re.search(r'version (\d+)\.(\d+)', output)
if version is None or tuple(map(int, version.groups())) < (3, 11):
    raise SystemExit('Setup requires CMake >=3.11.')
CHECK
if [[ -n "${CFITSIO_PREFIX:-}" ]]; then
    "$python" - "$CFITSIO_PREFIX" <<'CHECK'
from pathlib import Path
import sys
prefix = Path(sys.argv[1])
libraries = [p for directory in ('lib', 'lib64') for p in (prefix / directory).glob('libcfitsio.so*')]
if not (prefix / 'include/fitsio.h').is_file() or not libraries:
    raise SystemExit('CFITSIO_PREFIX must contain include/fitsio.h and lib[64]/libcfitsio.so. See docs/INSTALLATION.md.')
CHECK
else
    pkg-config --exists cfitsio || fail 'CFITSIO was not found. Install its development package or set CFITSIO_PREFIX (see docs/INSTALLATION.md).'
fi
printf 'Prerequisites passed. Environment: %s\n' "$venv"
if ((check_only)); then exit 0; fi
printf '\n[1/4] Creating isolated Python environment\n'
"$python" -m venv "$venv"
export PATH="$venv/bin:$PATH"
printf '\n[2/4] Installing Python dependencies and hstpipe\n'
"$venv/bin/python" -m pip install --upgrade pip setuptools wheel
if [[ ${HSTPIPE_USE_LOCK:-0} == 1 ]]; then
    "$venv/bin/python" -m pip install -r "$repo/support/environments/erebus-python311.lock.txt"
fi
(cd "$repo" && "$venv/bin/python" -m pip install -r requirements.txt)
printf '\n[3/4] Building pinned HSTCAL into the environment\n'
bash "$repo/support/hstcal/build.sh" "$venv"
printf '\n[4/4] Checking the installed environment\n'
"$venv/bin/python" "$repo/support/check_environment.py" --hstcal "$venv/bin"
"$venv/bin/python" -m pip check
"$venv/bin/hstpipe" --help
"$venv/bin/python" -m pip freeze > "$venv/share/hstpipe/native/python-environment.lock.txt"
printf '\nSetup complete. Activate with:\n  source %q/bin/activate\n' "$venv"
