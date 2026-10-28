"""Stage hash-guarded limiter specialization build using immutable source links."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT = ROOT / 'gpu_energy_neutral_native_build_v4_seed3'
LINKS = {
    'source.tar': PARENT / 'source.tar',
    'patch_guard.py': PARENT / 'patch_guard.py',
    'patch_iterations.py': PARENT / 'patch_iterations.py',
}
COPIES = {
    'patch_limiter.py': DOCS / 'patch_kepler_limiter_only_v7_20261027.py',
    'runner.sh': DOCS / 'build_kepler_limiter_only_v7_lageunha_20261027.sh',
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(['bash', '-n', str(COPIES['runner.sh'])], check=True)
    expected = {name: checksum for checksum, name in
                (line.split(maxsplit=1) for line in
                 (PARENT / 'inputs.sha256').read_text().splitlines())}
    for name, path in LINKS.items():
        if digest(path) != expected[name]:
            raise ValueError(f'Parent input changed: {name}')
    dest.mkdir()
    hashes = {}
    for name, path in LINKS.items():
        (dest / name).symlink_to(path.resolve())
        hashes[name] = digest(dest / name)
    for name, path in COPIES.items():
        shutil.copy2(path, dest / name)
        hashes[name] = digest(dest / name)
    with (dest / 'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    with (dest / 'stage_manifest.json').open('x') as out:
        json.dump({'run': str(dest), 'parent': str(PARENT),
                   'inputs': hashes,
                   'scope': 'Experimental original clipper plus limiter=1 specialization; build only.'},
                  out, indent=2)
        out.write('\n')
    print(json.dumps({'run': str(dest), 'inputs': hashes}, indent=2))


if __name__ == '__main__':
    main()
