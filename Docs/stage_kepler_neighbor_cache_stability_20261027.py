"""Stage adjacency persistence on two already completed native replays."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS=Path(__file__).resolve().parent
ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PAIR=ROOT/'gpu_energy_neutral_seedcert_w0_weighted_pair_v1'
RUN=ROOT/'neighbor_cache_stability_v1'


def main():
    gate=json.loads((PAIR/'pair_audit.json').read_text())
    if not all(gate['checks'].values()):raise ValueError('Native device pair gate failed')
    sources={'job.slurm':DOCS/'kepler_neighbor_cache_stability_20261027.slurm'}
    for name in ('audit_kepler_neighbor_cache_stability_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        sources[name]=DOCS/name
    if RUN.exists():raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir();hashes={}
    for name,source in sources.items():
        shutil.copy2(source,RUN/name);hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Input changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    print(RUN)


if __name__=='__main__':main()
