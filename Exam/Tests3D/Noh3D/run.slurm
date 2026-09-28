#!/bin/bash
#SBATCH -J LAGFORCE_NOH3D
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH -t 12:00:00
#SBATCH -o slurm-%x-%j.out
##SBATCH -p <cpu-partition>
# Example (not submitted by Grokbot):
#   sbatch --export=ALL,N=64 run.slurm            # smoke run (R = 7 dx)
#   sbatch --export=ALL,N=128 run.slurm           # the pass/fail run (R = 14 dx)
set -u
source "${LAGEUNHA_SRC:-$HOME/LagEunha}/Exam/Tests3D/common/run3d_common.sh"
source "$TESTS3D/Noh3D/config.env"
N=${N:-128}
lag3d_preflight "$LAG3D_BIN" || exit 3
RUN=$LAG3D_ROOT/Noh3D_gfs_N${N}_L${LAG3D_L}
lag3d_setup "$RUN" || exit 1
cd "$RUN" || exit 1
$LAG3D_PY "$TESTS3D/Noh3D/make_ic.py" --n "$N" --L "$LAG3D_L" --tend "$LAG3D_TEND" -o ic.l3d > ic.log 2>&1 || { cat ic.log; exit 1; }
lag3d_params "$TESTS3D/common/params3d.dat.template" params3d.dat \
	MODEL=Hydro3D GAMMA=1.666667 COURANT=$LAG3D_COURANT NX=$N NY=$N NZ=$N \
	ICFILE=ic.l3d TEND=$LAG3D_TEND DUMPDT=$LAG3D_DUMPDT GRAVITY=0 EPS=0 || exit 1
export HYDRO_TSTOP=$LAG3D_TEND EUNHA_DUMP_DT=$LAG3D_DUMPDT
lag3d_mpirun; rc=$?
echo "[3D] Noh3D N=$N rc=$rc; analyse with:"
echo "  $LAG3D_PY $TESTS3D/Noh3D/analyze.py --snap $(ls "$RUN"/snap_*.l3d 2>/dev/null | tail -1) --ic $RUN/ic.l3d"
exit $rc
