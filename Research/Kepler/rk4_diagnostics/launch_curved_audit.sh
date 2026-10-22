#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
test "$#" -ge 1
runner=$1
shift
test -f "$runner"
envelope=$(mktemp -d "$root/rk4_diagnostics/curved_launch_XXXXXXXX")
cp -- "$runner" "$envelope/runner.sh"
cp -- "$root/run_env.sh" "$envelope/run_env.sh"
sha256sum "$envelope/runner.sh" "$envelope/run_env.sh" > "$envelope/source.sha256"
printf '%q ' "$runner" "$@" > "$envelope/arguments.txt"
printf '\n' >> "$envelope/arguments.txt"
chmod a-w "$envelope/runner.sh" "$envelope/run_env.sh"
export GFS_CURVED_FROZEN_RUNNER=1
export GFS_CURVED_FROZEN_ENV=$envelope/run_env.sh
export GFS_CURVED_LAUNCH_ENVELOPE=$envelope
printf 'Frozen runner: %s\n' "$envelope/runner.sh"
# exec ensures the mutable launcher is never read again after the task starts.
exec bash "$envelope/runner.sh" "$@"
