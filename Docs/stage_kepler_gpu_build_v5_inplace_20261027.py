"""Stage an isolated CPU-node CUDA build of the in-place raw clipper."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT = ROOT / 'gpu_energy_neutral_native_build_v4_seed3'
FILES = {
    'source.tar': PARENT / 'source.tar',
    'patch_guard.py': PARENT / 'patch_guard.py',
    'patch_iterations.py': PARENT / 'patch_iterations.py',
    'patch_inplace.py': DOCS / 'patch_kepler_inplace_raw_clipping_v5_20261027.py',
    'job.slurm': DOCS / 'kepler_gpu_energy_neutral_native_build_v5_inplace_20261027.slurm',
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(['bash', '-n', str(FILES['job.slurm'])], check=True)
    parent_inputs = {name: checksum for checksum, name in
                     (line.split(maxsplit=1) for line in
                      (PARENT / 'inputs.sha256').read_text().splitlines())}
    if digest(FILES['source.tar']) != parent_inputs.get('source.tar'):
        raise ValueError('Source tar changed')
    for name in ('patch_guard.py', 'patch_iterations.py'):
        if digest(FILES[name]) != parent_inputs.get(name):
            raise ValueError(f'{name} changed')
    dest.mkdir()
    hashes = {}
    for name, source in FILES.items():
        shutil.copy2(source, dest / name)
        hashes[name] = digest(dest / name)
    with (dest / 'inputs.sha256').open('x') as stream:
        for name, value in sorted(hashes.items()):
            stream.write(f'{value}  {name}\n')
    with (dest / 'stage_manifest.json').open('x') as stream:
        json.dump({'run': str(dest), 'inputs': hashes,
                   'parent_build': str(PARENT),
                   'scope': 'Isolated experimental CUDA build; no physics or speed acceptance until matched native run.'},
                  stream, indent=2)
        stream.write('\n')
    print(json.dumps({'destination': str(dest), 'inputs': hashes}, indent=2))


if __name__ == '__main__':
    main()
