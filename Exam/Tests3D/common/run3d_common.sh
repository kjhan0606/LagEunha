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

: "${LAGEUNHA_SRC:=$HOME/LagEunha}"        # checkout that holds eunha2.exe
: "${LAG3D_BIN:=$LAGEUNHA_SRC/eunha2.exe}"
: "${LAG3D_ROOT:=/gpfs/kjhan/LagForce}"
: "${LAG3D_PY:=python3}"                    # needs numpy+scipy (IC generators)
: "${NRANK:=${SLURM_NTASKS:-4}}"
TESTS3D="$LAGEUNHA_SRC/Exam/Tests3D"

# The 3D GFS driver does not exist at the commit that added this suite (see
# Exam/Tests3D/README.md, "3D path audit"). A binary that has one must embed
# the string LAGEUNHA_3D_GFS_V1, e.g.
#   __attribute__((used)) static const char lag3d_marker[] = "LAGEUNHA_3D_GFS_V1";
# and accept "define Hydro3D IC file = <file>" in params3d.dat.
# Until then every run script stops here instead of burning an allocation.
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
	env | grep -E '^(SEDOV_PHASE1|GFS_|EUNHA_|HYDRO_|LAG3D_|OMP_NUM_THREADS)' | sort > "$run/ENV"
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

# lag3d_mpirun : run in the current directory, log, return mpirun's rc
lag3d_mpirun() {
	mpirun -np "$NRANK" ./eunha2.exe params3d.dat > log 2>&1
	local rc=$?
	echo "EXIT:$rc" >> log
	return $rc
}
