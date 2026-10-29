#!/bin/bash
set -euo pipefail
cd /gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029
unset SLURM_JOB_ID SLURM_CPUS_PER_TASK I_MPI_HYDRA_BOOTSTRAP
export CUDA_VISIBLE_DEVICES=GPU-42b3367f-1df7-2b2d-7961-7781c87b022d
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
sha256sum -c chain.sha256
/home/kjhan/miniconda3/bin/python3 /gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/run_kepler_q1_orbital_cap_chain_20261029.py /gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029
