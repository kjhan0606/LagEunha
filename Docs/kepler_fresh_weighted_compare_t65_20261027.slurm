#!/usr/bin/env bash
#SBATCH --job-name=LagQ2WeightedCompare65
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 - <<'PY_GATE'
import json
from pathlib import Path
import shutil
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

parent=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1')
source=parent/'audit_t65.json'
if not source.exists():raise ValueError('No complete t=6.5 audit; inspect predecessor failure')
a=json.loads(Path('w0_receipt.json').read_text());b=json.loads(source.read_text())
failed={k for k,v in b['checks'].items() if not v}
if (not all(a['checks'].values()) or len(b['checks'])!=21 or
        failed not in (set(),{'disk_density_equilibrium'}) or
        a['initial_sha256']!=b['initial_sha256'] or abs(a['time']-6.5)>1e-12 or
        abs(b['time']-6.5)>1e-12 or
        b['binary_sha256']!='feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b'):
    raise ValueError(f'Unexpected weighted failure or lineage: {sorted(failed)}')
for key in ('initial','endpoint'):
    if digest(Path(b[key]))!=b[key+'_sha256']:raise ValueError(f'Changed weighted {key}')
if abs(snapshot_time(Path(b['endpoint']),55351)-6.5)>1e-10:
    raise ValueError('Incomplete t=6.5 weighted endpoint')
target=Path('weighted_receipt.json')
if target.exists():raise FileExistsError(target)
shutil.copy2(source,target)
if digest(target)!=digest(source):raise ValueError('Copied receipt changed')
with Path('gate_decision.json').open('x') as out:
    json.dump({'weighted_failed_checks':sorted(failed),'weighted_physics_pass':not failed},out,indent=2)
    out.write('\n')
PY_GATE
failed=$(/home/kjhan/miniconda3/bin/python3 -c 'import json; print(len(json.load(open("gate_decision.json"))["weighted_failed_checks"]))')
/home/kjhan/miniconda3/bin/python3 audit_kepler_bounded_entropy_sign_20261027.py w0_receipt.json --output w0_entropy.json
if [[ "$failed" == 1 ]]; then
  /home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_entropy_sign_20261027.py weighted_receipt.json --allow-failed-check disk_density_equilibrium --output weighted_entropy.json
  /home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_vs_w0_fresh_20261027.py w0_receipt.json weighted_receipt.json w0_entropy.json weighted_entropy.json --allow-weighted-failed-check disk_density_equilibrium --output comparison.json
else
  /home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_entropy_sign_20261027.py weighted_receipt.json --output weighted_entropy.json
  /home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_vs_w0_fresh_20261027.py w0_receipt.json weighted_receipt.json w0_entropy.json weighted_entropy.json --output comparison.json
fi
