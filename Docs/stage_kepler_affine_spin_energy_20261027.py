"""Stage affine spin/strain diagnostic from matching accepted snapshots."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/affine_spin_energy_t1651_v1')


def main():
    directory=DOCS/'KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27'
    names={'job.slurm':DOCS/'kepler_affine_spin_energy_20261027.slurm',
           'w0_receipt.json':directory/'w0_jitter_step427_audit.json',
           'weighted_receipt.json':directory/'live_step427_audit.json'}
    for name in ('audit_kepler_affine_spin_energy_20261027.py',
                 'test_kepler_affine_spin_energy_20261027.py',
                 'audit_kepler_polygon_moment_native_faces_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        names[name]=DOCS/name
    a=json.loads(names['w0_receipt.json'].read_text());b=json.loads(names['weighted_receipt.json'].read_text())
    if (not all(a['checks'].values()) or not all(b['checks'].values()) or
            a['initial_sha256']!=b['initial_sha256'] or abs(a['time']-b['time'])>1e-6):
        raise ValueError('Expected matching accepted native checkpoints')
    if RUN.exists():raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(names['job.slurm'])],check=True)
    RUN.mkdir();hashes={}
    for name,source in names.items():
        shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Input changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)


if __name__=='__main__':main()
