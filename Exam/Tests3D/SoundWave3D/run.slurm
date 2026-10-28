#!/bin/bash
#SBATCH -J LAGFORCE_SW3D
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH -t 04:00:00
#SBATCH -o slurm-%x-%j.out
##SBATCH -p <cpu-partition>      # the 3D path is CPU-only; do not hold GPUs for it
# Example (not submitted by Grokbot):
#   sbatch --export=ALL,N=32 run.slurm
#   sbatch --export=ALL,N=64 run.slurm
#   sbatch --export=ALL,N=128 -t 12:00:00 run.slurm
#   sbatch --export=ALL,N=64,BOOST=1 run.slurm        # Galilean-invariance run
set -u
source "${LAGEUNHA_SRC:-$HOME/LagEunha}/Exam/Tests3D/common/run3d_common.sh"
source "$TESTS3D/SoundWave3D/config.env"
N=${N:-64}; BOOST=${BOOST:-0}
lag3d_preflight "$LAG3D_BIN" || exit 3
TAG=""; BARG=""
if [ "$BOOST" != "0" ]; then TAG="_boost$BOOST"; BARG="--boost $BOOST"; fi
RUN=$LAG3D_ROOT/SoundWave3D_gfs_N${N}${TAG}
lag3d_setup "$RUN" || exit 1
cd "$RUN" || exit 1
$LAG3D_PY "$TESTS3D/SoundWave3D/make_ic.py" --n "$N" --amp "$LAG3D_AMP" $BARG -o ic.l3d > ic.log 2>&1 || { cat ic.log; exit 1; }
lag3d_params "$TESTS3D/common/params3d.dat.template" params3d.dat \
	MODEL=Hydro3D GAMMA=1.666667 COURANT=$LAG3D_COURANT NX=$N NY=$N NZ=$N \
	ICFILE=ic.l3d TEND=$LAG3D_TEND DUMPDT=$LAG3D_DUMPDT GRAVITY=0 EPS=0 || exit 1
export HYDRO_TSTOP=$LAG3D_TEND EUNHA_DUMP_DT=$LAG3D_DUMPDT
lag3d_mpirun; rc=$?
echo "[3D] SoundWave3D N=$N boost=$BOOST rc=$rc; analyse with:"
echo "  $LAG3D_PY $TESTS3D/SoundWave3D/analyze.py --snaps <run_N32>/<last snap_NNNNNN.l3d> <run_N64>/... <run_N128>/... --boosted <run_N64_boost1>/<last snap>"
exit $rc
