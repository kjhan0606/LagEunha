#!/usr/bin/env bash
#SBATCH --job-name=LagQ2SpinCompile
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
g++ -std=c++14 -O2 -Wall -Wextra -Werror -I . host_reconstruct.cpp -o host_reconstruct
g++ -std=c++14 -O2 -Wall -Wextra -Werror -I . host_exchange.cpp -o host_exchange
./host_reconstruct
./host_exchange
nvcc -std=c++14 -O2 -arch=sm_80 gpu_reconstruct.cu -o gpu_reconstruct
nvcc -std=c++14 -O2 -arch=sm_80 gpu_exchange.cu -o gpu_exchange
sha256sum gpu_reconstruct gpu_exchange > binaries.sha256
