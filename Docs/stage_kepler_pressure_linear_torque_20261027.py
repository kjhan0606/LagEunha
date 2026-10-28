"""Stage native-polygon midpoint-versus-linear-pressure torque analysis."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/linear_pressure_torque_v1')


def main():
    sources={'job.slurm':DOCS/'kepler_pressure_linear_torque_20261027.slurm',
             'w0_receipt.json':DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json',
             'weighted_receipt.json':DOCS/'KEPLER_WEIGHTED_FRESH_CAPACITY_FAILURE_2026-10-27/live_step1703_audit.json'}
    for name in ('audit_kepler_pressure_linear_torque_20261027.py',
                 'test_kepler_pressure_linear_torque_20261027.py',
                 'audit_kepler_affine_spin_energy_20261027.py',
                 'audit_kepler_polygon_moment_native_faces_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        sources[name]=DOCS/name
    for name in ('w0_receipt.json','weighted_receipt.json'):
        r=json.loads(sources[name].read_text())
        if not all(r['checks'].values()):raise ValueError(f'Failed native parent: {name}')
    if RUN.exists():raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir();hashes={}
    for name,source in sources.items():
        shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)


if __name__=='__main__':main()
