"""Stage the accepted bounded Kepler t=.825 checkpoint for H200 timing."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time


SOURCE = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to65_near32_a10_v1')
LOCAL = Path(__file__).resolve().parent
FILES = (
    'eunha2.exe', 'params.dat', 'glass.bin', 'areas.bin',
    'audit_kepler_bounded_native_evolution_20261027.py',
    'audit_kepler_sparse_area_t0825_20261027.py',
    'audit_kepler_sparse_native_ballistic_20261027.py',
    'audit_kepler_sparse_state_consistency_20261027.py',
    'prepare_kepler_sparse_area_restart_20261027.py',
)
RUNNER = 'kepler_q2_bounded_h200_speed_to165137_20261027.slurm'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(f'Refusing existing stage: {dest}')
    if abs(snapshot_time(SOURCE / 'khout.000001.dat', 55351) - .825) > 1e-10:
        raise ValueError('Source is not the complete accepted t=.825 checkpoint')
    if digest(SOURCE / 'eunha2.exe') != '3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e':
        raise ValueError('Source executable differs from the audited A10 baseline')
    subprocess.run(['bash', '-n', str(LOCAL / RUNNER)], check=True)
    dest.mkdir(parents=True)
    for name in FILES:
        shutil.copy2(SOURCE / name, dest / name)
    shutil.copy2(SOURCE / 'khout.000001.dat', dest / 'khout.000001.dat')
    shutil.copy2(LOCAL / RUNNER, dest / RUNNER)
    with (dest / 'inputs.sha256').open('x') as stream:
        for name in (*FILES, 'khout.000001.dat', RUNNER):
            stream.write(f'{digest(dest / name)}  {name}\n')
    manifest = {
        'source': str(SOURCE),
        'destination': str(dest),
        'time_range': [.825, 1.651374332206546],
        'purpose': 'same-binary A10/H200 timing and physics comparison',
        'runner': RUNNER,
        'source_checkpoint_sha256': digest(SOURCE / 'khout.000001.dat'),
        'binary_sha256': digest(dest / 'eunha2.exe'),
    }
    with (dest / 'stage_manifest.json').open('x') as stream:
        json.dump(manifest, stream, indent=2)
        stream.write('\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
