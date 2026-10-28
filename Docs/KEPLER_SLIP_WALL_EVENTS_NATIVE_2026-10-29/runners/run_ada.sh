#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
if [[ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]];then echo 'GPU occupied' >&2;exit 31;fi
(cd /gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/native_events_v1/flag_off;sha256sum -c inputs.sha256;bash run.sh)
(cd /gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/native_events_v1/events_no_hit;sha256sum -c inputs.sha256;bash run.sh)
(cd /gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/native_events_v1/forced_retry;sha256sum -c inputs.sha256;bash run.sh)
(cd /gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/native_events_v1/forced_failure;sha256sum -c inputs.sha256;bash run.sh)
