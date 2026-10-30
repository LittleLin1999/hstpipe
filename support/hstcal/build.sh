#!/usr/bin/env bash
set -euo pipefail
# Usage: build.sh /absolute/software/root; requires git, CMake, C compiler, CFITSIO.
root=$(realpath -m "${1:?Supply an installation root outside observational data}")
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python3 - "$script_dir" <<'CHECK'
import pathlib,json,hashlib,sys
root=pathlib.Path(sys.argv[1]);meta=json.loads((root/'provenance.json').read_text())
if hashlib.sha256((root/'hstcal_coordinate_export.patch').read_bytes()).hexdigest()!=meta['patch_sha256']:
    raise SystemExit('Coordinate-export patch checksum differs from pinned provenance')
CHECK
commit=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["commit"])' "$script_dir/provenance.json")
mkdir -p "$root"
if [ ! -d "$root/hstcal/.git" ]; then git clone https://github.com/spacetelescope/hstcal.git "$root/hstcal"; fi
test "$(git -C "$root/hstcal" rev-parse HEAD)" = "$commit" || git -C "$root/hstcal" checkout --detach "$commit"
if ! git -C "$root/hstcal" apply --reverse --check "$script_dir/hstcal_coordinate_export.patch" 2>/dev/null; then
 git -C "$root/hstcal" apply --check "$script_dir/hstcal_coordinate_export.patch"
 git -C "$root/hstcal" apply "$script_dir/hstcal_coordinate_export.patch"
fi
cd "$root/hstcal"
cmake -S "$root/hstcal" -B "$root/hstcal-build" -DCMAKE_INSTALL_PREFIX="$root/hstcal-install" -DENABLE_OPENMP=OFF ${CFITSIO_PREFIX:+-DWITH_CFITSIO="$CFITSIO_PREFIX"}
cmake --build "$root/hstcal-build" --parallel 4
cmake --install "$root/hstcal-build"
sha256sum "$script_dir/hstcal_coordinate_export.patch" > "$root/hstcal-patch.sha256"
git -C "$root/hstcal" rev-parse HEAD > "$root/hstcal-commit.txt"
