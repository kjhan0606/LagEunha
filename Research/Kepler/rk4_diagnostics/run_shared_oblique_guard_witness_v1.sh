#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
diag=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
dest=$diag/shared_oblique_guard_witness_${1:-v1}_mpi1
mkdir "$dest"
source "$GFS_CURVED_FROZEN_ENV"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
previous=$diag/shared_oblique_panels_v1_mpi1
for input in "$previous/"*.h "$previous/"*.cpp; do cp "$input" "$dest/"; done
cp "$diag/test_shared_oblique_guard_witness_v1.cpp" "$dest/"
sha256sum -c "$previous/source.sha256" > "$dest/previous_source_integrity.txt"
sha256sum -c "$previous/output.sha256" > "$dest/previous_output_integrity.txt"
python3 - "$dest" <<'PY'
from pathlib import Path
import difflib,sys
p=Path(sys.argv[1])/'material_shared_oblique_point.h';old=p.read_text();new=old
changes=[('need(sf::lower(out.s)>0,"collected oblique regular-Duffy radial pole excluded");','if(!(sf::lower(out.s)>0))throw sd::DomainFailure("collected oblique regular-Duffy radial pole excluded",out.s,Q(0));')]
for phase,expression in [(70,'os::graph(m,p,f,u0,lo,hi,real)'),(71,'face(p,f,l,u0,S(root.n))'),(72,'face(p,f,l,u,S(root.n))'),(73,'face(p,f,l,u,S(predictor.n,sf::add(predictor.e,prior)),&local_s)'),(74,'sd::evaluate(m,f.original,c.s,c.t,&c.angular)')]:
 changes.append((expression,'[&]{sd::TracePhase observer_phase('+str(phase)+');return '+expression+';}()'))
for before,after in changes:
 assert new.count(before)==1,before;new=new.replace(before,after)
reverse=new
for before,after in reversed(changes):
 assert reverse.count(after)==1;reverse=reverse.replace(after,before)
assert reverse==old;p.write_text(new)
(p.parent/'READ_ONLY_GUARD_OBSERVER.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='original/'+p.name,tofile='observer/'+p.name)))
PY
packet=$diag/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$diag/shared_trace_rhs_v1_mpi4
sites=$diag/curved_map_algebra_v2/advected_sites.csv
weights=$diag/curved_events_v9/site_weights.csv
pressure=$diag/geometric_pressure_v2_mpi4
sha256sum "$dest/"*.h "$dest/"*.cpp "$dest/READ_ONLY_GUARD_OBSERVER.patch" > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$previous/candidate_oblique_panels_rank0.bin" "$previous/candidate_axis_proofs_rank0.bin" "$previous/unchanged_oblique_frames_rank0.bin" "$previous/unchanged_actual_slabs_rank0.bin" "$previous/same_chart_linear_coordinates_rank0.bin" "$pressure/geometric_cuda.o" > "$dest/input.sha256"
printf '%s\n' "$GFS_CURVED_LAUNCH_ENVELOPE" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_oblique_guard_witness_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 1 "$dest/contract.exe" "$packet" "$saved" "$sites" "$weights" "$previous" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv "$dest/"*.bin > "$dest/output.sha256"
printf 'OBLIQUE_FIRST_GUARD_WITNESS_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
