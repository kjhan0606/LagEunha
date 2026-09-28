#!/bin/bash
# Sourced by Exam/Tests3D/*/run.slurm. Mirrors the 2D launch style in
# code_review_grok.md: LAGFORCE_* job names, one directory per run under
# /gpfs/kjhan/LagForce, binary copied into the run directory, params file,
# environment flags exported by the script, stdout+stderr to ./log,
# "EXIT:<rc>" appended to the log.
# Difference from the 2D scripts: the batch script exits with mpirun's return
# code, so sacct no longer reports COMPLETED 0:0 for a crashed run.
#
# Nothing here submits or cancels jobs.

# Binary: eunha2.exe dispatches "Simulation Model = Hydro3D" to the 3D GFS
# path (Exam/exam3d_gfs.c, linked through libexam.a). The same source also
# builds standalone without MPI/FFTW:
#   cd $LAGEUNHA_SRC/Exam && make lag3d.exe      (then LAG3D_BIN=$LAGEUNHA_SRC/Exam/lag3d.exe)
: "${LAGEUNHA_SRC:=$HOME/LagEunha}"        # checkout that holds eunha2.exe
: "${LAG3D_BIN:=$LAGEUNHA_SRC/eunha2.exe}"
: "${LAG3D_ROOT:=/gpfs/kjhan/LagForce}"
: "${LAG3D_PY:=python3}"                    # needs numpy+scipy (IC generators)
# The 3D path is one MPI rank + OpenMP (extra ranks would only idle), so the
# default is 1 rank; threads come from OMP_NUM_THREADS (flags_base.env).
: "${NRANK:=1}"
# LAG3D_LAUNCH=mpirun (default) or direct (run the binary without mpirun;
# fine for lag3d.exe, and for eunha2.exe where singleton MPI init works).
: "${LAG3D_LAUNCH:=mpirun}"
TESTS3D="$LAGEUNHA_SRC/Exam/Tests3D"

# A binary that has the 3D GFS driver embeds the string LAGEUNHA_3D_GFS_V1
# (Exam/exam3d_gfs.c) and reads "define Hydro3D IC file = <file>" from
# params3d.dat. Anything else (an old eunha2.exe, a typo in LAG3D_BIN) stops
# here instead of burning an allocation.
lag3d_preflight() {
	local bin="$1"
	if [ ! -x "$bin" ]; then
		echo "[3D] binary $bin not found or not executable"; return 3
	fi
	if ! grep -a -q "LAGEUNHA_3D_GFS_V1" "$bin"; then
		echo "[3D] $bin has no 3D GFS driver (marker LAGEUNHA_3D_GFS_V1 absent)."
		echo "[3D] See Exam/Tests3D/README.md. Not running."
		return 3
	fi
	return 0
}

# lag3d_setup <rundir> : make the run dir, copy binary, record provenance
lag3d_setup() {
	local run="$1"
	mkdir -p "$run" || return 1
	cp "$LAG3D_BIN" "$run/eunha2.exe" || return 1
	( cd "$LAGEUNHA_SRC" && git rev-parse HEAD 2>/dev/null ) > "$run/COMMIT"
	env | grep -E '^(SEDOV_PHASE1|GFS_|GFS3D_|EUNHA_|HYDRO_|LAG3D_|OMP_NUM_THREADS)' | sort > "$run/ENV"
	echo "[3D] run dir $run  commit $(cat "$run/COMMIT")  ranks $NRANK"
}

# lag3d_params <template> <out> KEY=VALUE ... : fill @KEY@ placeholders
lag3d_params() {
	local tpl="$1" out="$2"; shift 2
	cp "$tpl" "$out"
	local kv
	for kv in "$@"; do sed -i "s|@${kv%%=*}@|${kv#*=}|g" "$out"; done
	if grep -q '@[A-Z_]*@' "$out"; then
		echo "[3D] unfilled placeholders in $out:"; grep '@[A-Z_]*@' "$out"; return 1
	fi
}

# lag3d_mpirun : run in the current directory, log, return the launcher's rc.
# One rank with many OpenMP threads: switch off the launcher's core binding,
# or all threads share one core. Measured on the box (Open MPI 5.0.7, 8
# threads, SoundWave3D N=16): 12.2 s bound vs 1.8 s with binding off.
#   Open MPI 5 (PRRTE): PRTE_MCA_hwloc_default_binding_policy=none
#   Open MPI 4:         OMPI_MCA_hwloc_base_binding_policy=none
#   Intel MPI:          I_MPI_PIN_DOMAIN=omp (domain = OMP_NUM_THREADS cores)
# LAG3D_MPIRUN_OPTS adds launcher options (e.g. "--bind-to none").
lag3d_mpirun() {
	local rc
	export PRTE_MCA_hwloc_default_binding_policy=${PRTE_MCA_hwloc_default_binding_policy:-none}
	export OMPI_MCA_hwloc_base_binding_policy=${OMPI_MCA_hwloc_base_binding_policy:-none}
	export I_MPI_PIN_DOMAIN=${I_MPI_PIN_DOMAIN:-omp}
	echo "[3D] launch=$LAG3D_LAUNCH ranks=$NRANK OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}" > log
	if [ "$LAG3D_LAUNCH" = direct ]; then
		./eunha2.exe params3d.dat >> log 2>&1
	else
		mpirun -np "$NRANK" ${LAG3D_MPIRUN_OPTS:-} ./eunha2.exe params3d.dat >> log 2>&1
	fi
	rc=$?
	echo "EXIT:$rc" >> log
	return $rc
}
