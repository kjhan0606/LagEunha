#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
if test "${RADIAL_FROZEN_RUNNER:-0}" != 1; then
    case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
    envelope=$(mktemp -d "$root/radial_launch_XXXXXXXX")
    cp -- "$0" "$envelope/runner.sh"
    cp -- "$root/prepare_radial_resolution_v1.py" "$root/audit_snapshot.py" "$root/snapshot_offsets.json" "$envelope/"
    sha256sum "$envelope/"* > "$envelope/source.sha256"
    chmod a-w "$envelope/"*
    export RADIAL_FROZEN_RUNNER=1 RADIAL_ENVELOPE=$envelope
    exec bash "$envelope/runner.sh" "$@"
fi
version=${1:-v2}
case "$version" in v[0-9]*) ;; *) exit 2;; esac
dest=$root/radial_resolution_20261021_$version
test ! -e "$dest"
sha256sum "$root/core_refinement/coarse_integrated/khout.000001.dat" "$root/core_refinement/coarse_integrated/params.dat" > "$RADIAL_ENVELOPE/inputs.sha256"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
set +e
python3 "$RADIAL_ENVELOPE/prepare_radial_resolution_v1.py" --historical-root "$root" --output "$dest" > "$RADIAL_ENVELOPE/preparation.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$RADIAL_ENVELOPE/exit_code.txt"
printf '%s\n' "$RADIAL_ENVELOPE" > "$dest/preparation_envelope.txt"
sha256sum -c "$RADIAL_ENVELOPE/source.sha256" > "$RADIAL_ENVELOPE/final_source_integrity.txt"
exit "$result"
