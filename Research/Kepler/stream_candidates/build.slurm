#!/usr/bin/env bash
#SBATCH --job-name=LAGFORCE_stream_build
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/build-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/build-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3
module load intel/compiler-rt/2025.3.0
module load intel/umf/1.0.2
module load intel/compiler/2025.3.0
module load intel/mpi/2021.17
module load gnu13/13.2.0
module load cuda/13.0.2

cd /gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
make -C Exam -j1 USE_CUDA=1 libexam.a > /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/build_exam.log 2>&1
make -B -j1 USE_CUDA=1 eunha2.exe > /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/build_link.log 2>&1
test -x eunha2.exe
cp eunha2.exe /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/eunha2_stream.exe
sha256sum /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/eunha2_stream.exe > /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/binary.sha256
ldd /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/eunha2_stream.exe > /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/ldd.txt
if grep -q 'not found' /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/ldd.txt; then exit 2; fi
