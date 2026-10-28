"""Prestage an independent half-timestep H200 run from accepted t=6.5."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
PARENT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/jitter005_to65_a10_v1')
AUDIT = DOCS / 'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json'
FILES = (
    'kepler_q2_jitter_halfdt_to10_h200_20261027.slurm',
    'prepare_kepler_bounded_followon_20261027.py',
    'prepare_kepler_sparse_area_restart_20261027.py',
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    d = json.loads(AUDIT.read_text())
    endpoint = Path(d['endpoint'])
    if (d['time'] != 6.5 or not all(d['checks'].values())
            or endpoint.parent != PARENT
            or digest(endpoint) != d['endpoint_sha256']
            or abs(snapshot_time(endpoint,55351)-6.5)>1e-10):
        raise ValueError('Accepted t=6.5 jitter parent changed')
    runner = DOCS / FILES[0]
    subprocess.run(['bash','-n',str(runner)],check=True)
    text = runner.read_text()
    if ('#SBATCH --gres=gpu:H200:1' not in text
            or 'EUNHA_KEPLER_DTETA=0.03' not in text
            or 'HYDRO_TSTOP=10' not in text):
        raise ValueError('Unexpected half-dt H200 runner')
    dest.mkdir()
    hashes = {}
    for name in FILES:
        shutil.copy2(DOCS / name, dest / name)
        hashes[name] = digest(dest / name)
    with (dest / 'prestage_manifest.json').open('x') as out:
        json.dump({'run':str(dest),'parent':str(PARENT),
                   'parent_audit_sha256':digest(AUDIT),
                   'parent_endpoint_sha256':d['endpoint_sha256'],
                   'time_start':6.5,'time_end':10.0,
                   'delta':{'EUNHA_KEPLER_DTETA':[0.06,0.03]},
                   'prestaged_sha256':hashes,
                   'scope':'One-rank H200 temporal-convergence control; same accepted W0 jitter t6.5 state, no t0 restart.'},
                  out,indent=2)
        out.write('\n')
    print(json.dumps({'run':str(dest),'inputs':hashes},indent=2))


if __name__=='__main__':
    main()
