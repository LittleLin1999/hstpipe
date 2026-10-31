#!/usr/bin/env bash
# Run from any directory; setup logic lives in the packaged support helper.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec bash "$repo/support/bootstrap.sh" "$@"
