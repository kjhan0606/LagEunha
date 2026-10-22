#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
guard=${2:-$root/rk4_diagnostics/curved_topology_v1}
dest=$root/rk4_diagnostics/curved_adaptive_transfer_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
python3 - "$guard" <<'PY'
import json,pathlib,sys
p=pathlib.Path(sys.argv[1])
assert (p/'status.txt').read_text().strip()=='AUDIT_COMPLETE'
d=json.loads((p/'contract.json').read_text())
assert d['all_event_bounds_pass'], 'radial event guard failed; preserve witness and do not promote'
PY
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=2 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for h in gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h; do cp "$guard/$h" "$dest/"; done
for f in material_curved_packet.h material_curved_bounds.h material_curved_roots.h material_curved_trace_bounds.h material_curved_partition.h material_curved_topology.h; do cp "$guard/$f" "$dest/"; done
if test -f "$guard/material_curved_compiled.h"; then cp "$guard/material_curved_compiled.h" "$dest/";fi
for f in test_curved_transfer_v1.cpp test_curved_adaptive_transfer_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
sha256sum -c "$guard/inputs.sha256" > "$dest/guard_inputs_unchanged.txt"
sha256sum "$guard/contract.json" "$guard"/radial_panels_rank*.csv > "$dest/guard.sha256"
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_adaptive_transfer_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
mpirun -np 4 "$dest/contract.exe" "$packet" "$sites" "$guard" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
