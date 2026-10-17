#!/bin/bash
#SBATCH -J LAGFORCE_KH3D
#SBATCH -N 1
#SBATCH --ntasks=1                # 3D path: 1 rank + OpenMP (NRANK=1 in run3d_common.sh)
#SBATCH --cpus-per-task=32        # = OMP_NUM_THREADS (flags_base.env)
#SBATCH -t 12:00:00
#SBATCH -o slurm-%x-%j.out
##SBATCH -p <cpu-partition>
# Example (not submitted by Grokbot):
#   sbatch --export=ALL,N=64 run.slurm             # 64 x 64 x 16
#   sbatch --export=ALL,N=128 run.slurm            # 128 x 128 x 32, the pass/fail run
#   sbatch --export=ALL,N=256 -t 48:00:00 run.slurm
set -u
source "${LAGEUNHA_SRC:-$HOME/LagEunha}/Exam/Tests3D/common/run3d_common.sh"
source "$TESTS3D/KH3D/config.env"
N=${N:-128}
NZ=$(awk -v n="$N" -v l="$LAG3D_LZ" 'BEGIN{v=int(n*l+0.5); if(v<1)v=1; print v}')
lag3d_preflight "$LAG3D_BIN" || exit 3
RUN=$LAG3D_ROOT/KH3D_gfs_N${N}_nz${NZ}
lag3d_setup "$RUN" || exit 1
cd "$RUN" || exit 1
$LAG3D_PY "$TESTS3D/KH3D/make_ic.py" --n "$N" --lz "$LAG3D_LZ" --tend "$LAG3D_TEND" -o ic.l3d > ic.log 2>&1 || { cat ic.log; exit 1; }
lag3d_params "$TESTS3D/common/params3d.dat.template" params3d.dat \
	MODEL=Hydro3D GAMMA=1.666667 COURANT=$LAG3D_COURANT NX=$N NY=$N NZ=$NZ \
	ICFILE=ic.l3d TEND=$LAG3D_TEND DUMPDT=$LAG3D_DUMPDT GRAVITY=0 EPS=0 || exit 1
export HYDRO_TSTOP=$LAG3D_TEND EUNHA_DUMP_DT=$LAG3D_DUMPDT
lag3d_mpirun; rc=$?
echo "[3D] KH3D N=$N rc=$rc; analyse with:"
echo "  $LAG3D_PY $TESTS3D/KH3D/analyze.py --snaps '$RUN/snap_*.l3d' --ic $RUN/ic.l3d"
exit $rc
