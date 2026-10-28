#!/usr/bin/env bash
#SBATCH --job-name=LagQ2AffineSpin65
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 - <<'PY_PARENT'
import json
from pathlib import Path
import shutil
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

parent=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1')
source=parent/'audit_t65.json'
a=json.loads(Path('w0_receipt.json').read_text());b=json.loads(source.read_text())
if (not all(a['checks'].values()) or len(b['checks'])!=21 or not all(b['checks'].values())
        or a['initial_sha256']!=b['initial_sha256'] or abs(a['time']-6.5)>1e-12
        or abs(b['time']-6.5)>1e-12
        or b['binary_sha256']!='feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b'):
    raise ValueError('Parent t=6.5 physics gate or common initial state failed')
for key in ('initial','endpoint'):
    if digest(Path(b[key]))!=b[key+'_sha256']:raise ValueError(f'Changed weighted {key}')
if abs(snapshot_time(Path(b['endpoint']),55351)-6.5)>1e-10:
    raise ValueError('Incomplete weighted t=6.5 checkpoint')
target=Path('weighted_receipt.json')
if target.exists():raise FileExistsError(target)
shutil.copy2(source,target)
if digest(target)!=digest(source):raise ValueError('Copied receipt changed')
PY_PARENT
/home/kjhan/miniconda3/bin/python3 test_kepler_affine_spin_energy_20261027.py
/home/kjhan/miniconda3/bin/python3 audit_kepler_affine_spin_energy_20261027.py w0_receipt.json --output w0.json
/home/kjhan/miniconda3/bin/python3 audit_kepler_affine_spin_energy_20261027.py weighted_receipt.json --output weighted.json
