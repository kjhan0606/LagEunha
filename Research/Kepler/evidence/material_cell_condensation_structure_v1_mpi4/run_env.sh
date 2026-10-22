module purge
module load intel/tbb/2022.3
module load intel/compiler-rt/2025.3.0
module load intel/umf/1.0.2
module load intel/compiler/2025.3.0
module load intel/mpi/2021.17
module load gnu13/13.2.0
module load cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
# An inherited OMP_PLACES=threads otherwise makes libgomp bind every unpinned
# MPI rank to the same first two hardware threads. Preserve explicit choices.
export OMP_PROC_BIND=${OMP_PROC_BIND:-false}
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export EUNHA_KEPLER_CX=12 EUNHA_KEPLER_CY=12
export EUNHA_KEPLER_RIN=2 EUNHA_KEPLER_ROUT=8
export EUNHA_KEPLER_EPS=0.047 EUNHA_KEPLER_DTETA=0.06
export HYDRO_TSTOP=0.2 EUNHA_DUMP_DT=0.2
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1
export GFS_GPU_FACE_DUMP_CALL=100 GFS_KEPLER_HALO_CELLS=8
unset GFS_GPU_SYNC_HALO_STATE
unset GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_VOL_AUDIT GFS_FACE_AUDIT GFS_FACE_CHARGE_LIMIT GFS_GEOM_FACE_VEL
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT
unset SEDOV_LAGVOL SEDOV_WSMOOTH EUNHA_KEPLER_HR EUNHA_KEPLER_PROFILE
