#!/usr/bin/env bash
set -euo pipefail
# Install pinned HSTCAL into a Python environment; discard temporary build files.
prefix=$(realpath -m "${1:?Supply the Python environment installation prefix}")
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python3 - "$script_dir" <<'CHECK'
import pathlib,json,hashlib,sys
root=pathlib.Path(sys.argv[1]);meta=json.loads((root/'provenance.json').read_text())
if hashlib.sha256((root/'hstcal_coordinate_export.patch').read_bytes()).hexdigest()!=meta['patch_sha256']:
    raise SystemExit('Coordinate-export patch checksum differs from pinned provenance')
CHECK
commit=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["commit"])' "$script_dir/provenance.json")
work=$(mktemp -d "${TMPDIR:-/tmp}/hstpipe-hstcal.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
mkdir -p "$prefix"
git clone https://github.com/spacetelescope/hstcal.git "$work/source"
git -C "$work/source" checkout --detach "$commit"
git -C "$work/source" apply --check "$script_dir/hstcal_coordinate_export.patch"
git -C "$work/source" apply "$script_dir/hstcal_coordinate_export.patch"
cd "$work/source"
cmake -S "$work/source" -B "$work/build" -DCMAKE_INSTALL_PREFIX="$prefix" -DENABLE_OPENMP=OFF ${CFITSIO_PREFIX:+-DWITH_CFITSIO="$CFITSIO_PREFIX"}
cmake --build "$work/build" --parallel 4
cmake --install "$work/build"
mkdir -p "$prefix/share/hstpipe/native"
sha256sum "$script_dir/hstcal_coordinate_export.patch" > "$prefix/share/hstpipe/native/hstcal-patch.sha256"
git -C "$work/source" rev-parse HEAD > "$prefix/share/hstpipe/native/hstcal-commit.txt"
