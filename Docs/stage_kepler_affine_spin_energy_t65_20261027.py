"""Stage a conditional matching W0/weighted t=6.5 spin diagnostic."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
PARENT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/affine_spin_energy_t1651_v1')
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/affine_spin_energy_t65_v1')


def main():
    parent_hashes={}
    for line in (PARENT/'inputs.sha256').read_text().splitlines():
        value,name=line.split('  ',1);parent_hashes[name]=value
    names={'job.slurm':DOCS/'kepler_affine_spin_energy_t65_20261027.slurm',
           'w0_receipt.json':DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json'}
    for name in ('audit_kepler_affine_spin_energy_20261027.py',
                 'test_kepler_affine_spin_energy_20261027.py',
                 'audit_kepler_polygon_moment_native_faces_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        if digest(PARENT/name)!=parent_hashes[name]:raise ValueError(f'Prior diagnostic changed: {name}')
        names[name]=PARENT/name
    w=json.loads(names['w0_receipt.json'].read_text())
    if not all(w['checks'].values()) or abs(w['time']-6.5)>1e-12:
        raise ValueError('W0 t=6.5 receipt failed')
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
