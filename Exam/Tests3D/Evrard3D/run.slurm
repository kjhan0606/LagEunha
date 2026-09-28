#!/bin/bash
#SBATCH -J LAGFORCE_EVR3D
#SBATCH -N 1
#SBATCH --ntasks=1                # 3D path: 1 rank + OpenMP (NRANK=1 in run3d_common.sh)
#SBATCH --cpus-per-task=32        # = OMP_NUM_THREADS (flags_base.env)
#SBATCH -t 24:00:00
#SBATCH -o slurm-%x-%j.out
##SBATCH -p <cpu-partition>
# Example (not submitted by Grokbot):
#   sbatch --export=ALL,N=40 run.slurm            # ~33k particles, smoke run
#   sbatch --export=ALL,N=64 run.slurm            # ~137k particles, the pass/fail run
# Self-gravity: Barnes-Hut octree + Plummer softening in Exam/exam3d_gfs.c
# (LAG3D_THETA, default 0.5; LAG3D_GRAV_DIRECT=1 for direct summation).
set -u
source "${LAGEUNHA_SRC:-$HOME/LagEunha}/Exam/Tests3D/common/run3d_common.sh"
source "$TESTS3D/Evrard3D/config.env"
N=${N:-64}; BG=${BACKGROUND:-$LAG3D_BACKGROUND}
lag3d_preflight "$LAG3D_BIN" || exit 3
RUN=$LAG3D_ROOT/Evrard3D_gfs_n${N}_${BG}
lag3d_setup "$RUN" || exit 1
cd "$RUN" || exit 1
$LAG3D_PY "$TESTS3D/Evrard3D/make_ic.py" --n "$N" --eps "$LAG3D_EPS" --background "$BG" --tend "$LAG3D_TEND" -o ic.l3d > ic.log 2>&1 || { cat ic.log; exit 1; }
lag3d_params "$TESTS3D/common/params3d.dat.template" params3d.dat \
	MODEL=Hydro3D GAMMA=1.666667 COURANT=$LAG3D_COURANT NX=$N NY=$N NZ=$N \
	ICFILE=ic.l3d TEND=$LAG3D_TEND DUMPDT=$LAG3D_DUMPDT GRAVITY=1 EPS=$LAG3D_EPS || exit 1
export HYDRO_TSTOP=$LAG3D_TEND EUNHA_DUMP_DT=$LAG3D_DUMPDT
lag3d_mpirun; rc=$?
echo "[3D] Evrard3D n=$N rc=$rc; analyse with:"
echo "  $LAG3D_PY $TESTS3D/Evrard3D/analyze.py --log $RUN/log --snap $RUN/<snapshot at t=0.8>.l3d"
exit $rc
