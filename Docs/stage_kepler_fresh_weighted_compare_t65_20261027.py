"""Stage afterany analysis which preserves explicit density-gate failure."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/fresh_weighted_compare_t65_v1')


def main():
    sources={'job.slurm':DOCS/'kepler_fresh_weighted_compare_t65_20261027.slurm',
             'w0_receipt.json':DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json'}
    for name in ('audit_kepler_bounded_entropy_sign_20261027.py',
                 'audit_kepler_weighted_entropy_sign_20261027.py',
                 'audit_kepler_weighted_vs_w0_fresh_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        sources[name]=DOCS/name
    a=json.loads(sources['w0_receipt.json'].read_text())
    if not all(a['checks'].values()) or abs(a['time']-6.5)>1e-12:
        raise ValueError('W0 t=6.5 receipt failed')
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
