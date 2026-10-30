#!/usr/bin/env bash
set -euo pipefail
# Use an explicit environment and external software root; never activate old data software.
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
venv=${1:?Supply an absolute virtualenv directory}
software=${2:?Supply an absolute native software directory}
python=${PYTHON:-python3}
"$python" -m venv "$venv"
"$venv/bin/python" -m pip install --upgrade pip setuptools wheel
if [ "${HSTPIPE_USE_LOCK:-0}" = 1 ]; then
 "$venv/bin/python" -m pip install -r "$repo/support/environments/erebus-python311.lock.txt"
fi
"$venv/bin/python" -m pip install -e "$repo[all,test]"
bash "$repo/support/hstcal/build.sh" "$software"
"$venv/bin/python" -m pip check
"$venv/bin/python" -m pip freeze > "$software/python-environment.lock.txt"
