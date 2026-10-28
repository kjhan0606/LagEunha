#!/usr/bin/env bash
set -euo pipefail

baseline=/tmp/Kepler_sparse_fixedpoly_20261027
candidate=/tmp/Kepler_pressure_force_off_20261027
test "$(hostname)" = syntax
test "$(sha256sum "$baseline/eunha2.exe" | cut -d' ' -f1)" = \
  3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e
test "$(sha256sum "$candidate/Exam/exam_gpu.cu" | cut -d' ' -f1)" = \
  b17e681f4ef72e68450fd9bc95ec29c0bf606bc33331a03568a0d3410dd988c8
test "$(sha256sum "$candidate/Exam/Makefile" | cut -d' ' -f1)" = \
  4f45f553c2221fc253a83fa8600086003d3cc09f999569d81ba2ce81027c04b5

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
make -C "$candidate/Exam" -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 \
  GFS_DIAG_PRESSURE_FORCE_OFF=1 libexam.a
make -C "$candidate" -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 \
  GFS_DIAG_PRESSURE_FORCE_OFF=1 this
sha256sum "$baseline/eunha2.exe" "$candidate/Exam/exam_gpu.o" \
          "$candidate/eunha2.exe"
