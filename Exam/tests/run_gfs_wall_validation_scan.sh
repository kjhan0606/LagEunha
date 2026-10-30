#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "$0")" && pwd)
test_bin="/tmp/gfs_wall_validation_scan_test_$$"
trap 'rm -f "$test_bin"' EXIT
icx -std=c11 -O2 -qopenmp -Wall -Wextra -Werror \
  "$here/test_gfs_wall_validation_scan.c" -lm -o "$test_bin"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" "$test_bin"
