#!/bin/bash
set -euo pipefail
cd /gpfs/kjhan/LagForce/Kepler_SlipWall_Integration_20261029/project
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
make -C Exam -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 libexam.a > ../build_const_cache_exam.log 2>&1
make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 this > ../build_const_cache_link.log 2>&1
sha256sum eunha2.exe Exam/exam.c Exam/gfs_thermal64.h > ../build_const_cache_outputs.sha256
