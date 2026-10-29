#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES=GPU-42b3367f-1df7-2b2d-7961-7781c87b022d
bash /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/rhs_t10_selected_v1/run.sh
