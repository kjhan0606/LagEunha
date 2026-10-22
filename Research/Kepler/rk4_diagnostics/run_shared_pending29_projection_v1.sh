#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
diag=$root/rk4_diagnostics
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
dest=$diag/shared_pending29_projection_${1:-v1}_mpi1
mkdir "$dest"
source "$GFS_CURVED_FROZEN_ENV"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
witness=$diag/shared_atlas_witness_v1_mpi4
repair=$diag/shared_atlas_repairs_v1_mpi1
registry=$diag/shared_atlas_adoption_v5
actual=$diag/shared_registry_domains_v1_mpi4
for input in "$witness/"*.h "$witness/"*.cpp; do cp "$input" "$dest/"; done
cp "$diag/test_shared_pending29_projection_v1.cpp" "$diag/SHARED_PENDING29_PLAN_20261021.md" "$registry/remaining29_atlas_slots.csv" "$dest/"
sha256sum -c "$witness/source.sha256" > "$dest/previous_source_integrity.txt"
sha256sum -c "$witness/output.sha256" > "$dest/previous_output_integrity.txt"
sha256sum -c "$registry/output.sha256" > "$dest/registry_output_integrity.txt"
python3 - "$diag/shared_oblique_apex_v1/contract.json" "$dest" <<'PY'
from pathlib import Path
from fractions import Fraction as F
import csv,json,sys
source=Path(sys.argv[1]);dest=Path(sys.argv[2]);c=json.loads(source.read_text())
assert c['proved_strict_inner_apex_neighborhoods']==14 and c['outside_apex_count']==14
records=[]
for r in c['chart_results']:
    bound=F(r['inner_radial_cut_exact']);cut=F(1,2)
    while cut>=bound: cut/=2
    assert 0<cut<bound
    records.append(dict(slot_id=r['slot_id'],exact_binary_cut=format(float(cut),'.30f'),exact_proved_upper=str(bound),exact_apex_margin=r['strict_inner_proof_margin_exact']))
with (dest/'oblique_excluded_cuts.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
PY
packet=$diag/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$diag/shared_trace_rhs_v1_mpi4
sites=$diag/curved_map_algebra_v2/advected_sites.csv
weights=$diag/curved_events_v9/site_weights.csv
retained=$diag/shared_target_v3_mpi4
bounds=$diag/shared_real_bound_v1_mpi4
pressure=$diag/geometric_pressure_v2_mpi4
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$retained/source_target_contributions.csv" "$bounds/source_real_bounds.csv" "$witness/failed_charts_rank"*.bin "$actual/domain_nodes_rank"*.csv "$actual/domain_boundaries_rank"*.csv "$repair/oblique_direction_proofs.csv" "$pressure/geometric_cuda.o" "$dest/remaining29_atlas_slots.csv" "$dest/oblique_excluded_cuts.csv" "$dest/SHARED_PENDING29_PLAN_20261021.md" "$diag/shared_oblique_apex_v1/contract.json" "$diag/SHARED_ATLAS_BUDGET_POLICY_20261021.json" "$registry/contract.json" > "$dest/input.sha256"
printf '%s\n' "$GFS_CURVED_LAUNCH_ENVELOPE" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_pending29_projection_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 1 "$dest/contract.exe" "$packet" "$saved" "$sites" "$weights" "$retained" "$bounds" "$actual" "$witness" "$repair" "$dest/remaining29_atlas_slots.csv" "$dest/oblique_excluded_cuts.csv" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv "$dest/"*.bin > "$dest/output.sha256"
printf 'PENDING29_REFERENCE_PROJECTION_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
