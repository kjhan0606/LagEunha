"""Stage two exact native restarts for the RK4 timing-only A10 comparison."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time


ROOT = Path(__file__).resolve().parent
PARENT_AUDIT = ROOT / 'KEPLER_Q2_BOUNDED_NATIVE_T165137_2026-10-27.json'
PARENT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/to65_near32_a10_v1')
PROFILE = Path('/tmp/Kepler_bounded_rk_gap_20261027/eunha2.exe')
BASE_SHA = '3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e'
PROFILE_SHA = 'e10cca50e3e158baeca28fea142b9f061ddafd17aade3756a8bc42242ae9761c'
RUNNER = 'kepler_bounded_rk4_gap_short_a10_20261027.slurm'
AUDITOR = 'audit_kepler_bounded_rk4_gap_20261027.py'
READER = 'prepare_kepler_sparse_area_restart_20261027.py'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(f'Refusing existing stage: {dest}')
    receipt = json.loads(PARENT_AUDIT.read_text())
    state = Path(receipt['endpoint']).resolve()
    if (abs(receipt['time'] - 1.651374332206546) > 1e-10
            or not all(receipt['checks'].values()) or state.parent != PARENT
            or digest(state) != receipt['endpoint_sha256']
            or abs(snapshot_time(state, 55351) - receipt['time']) > 1e-10):
        raise ValueError('Parent checkpoint is not the accepted native t=1.651374 state')
    if digest(PARENT / 'eunha2.exe') != BASE_SHA or digest(PROFILE) != PROFILE_SHA:
        raise ValueError('Base/profile binary identity changed')
    subprocess.run(['bash', '-n', str(ROOT / RUNNER)], check=True)
    dest.mkdir(parents=True)
    for name in (RUNNER, AUDITOR, READER):
        shutil.copy2(ROOT / name, dest / name)
    for kind, binary in (('base', PARENT / 'eunha2.exe'), ('profile', PROFILE)):
        target = dest / kind
        target.mkdir()
        for name in ('params.dat', 'glass.bin', 'areas.bin'):
            shutil.copy2(PARENT / name, target / name)
        shutil.copy2(state, target / 'khout.000001.dat')
        shutil.copy2(binary, target / 'eunha2.exe')
    files = [RUNNER, AUDITOR, READER]
    files.extend(f'{kind}/{name}' for kind in ('base', 'profile')
                 for name in ('eunha2.exe', 'params.dat', 'glass.bin', 'areas.bin', 'khout.000001.dat'))
    with (dest / 'inputs.sha256').open('x') as stream:
        for name in files:
            stream.write(f'{digest(dest / name)}  {name}\n')
    manifest = {
        'source_audit': str(PARENT_AUDIT),
        'source_audit_sha256': digest(PARENT_AUDIT),
        'source_checkpoint': str(state),
        'source_checkpoint_sha256': receipt['endpoint_sha256'],
        'base_binary_sha256': BASE_SHA,
        'profile_binary_sha256': PROFILE_SHA,
        'profile_source_sha256': digest(Path('/tmp/Kepler_bounded_rk_gap_20261027/Exam/exam.c')),
        'profile_patch_sha256': digest(ROOT / 'KEPLER_RK4_GAP_PROFILE_2026-10-27.patch'),
        'target_time': 1.6785,
        'purpose': 'matched 14-step byte-equality and exclusive RK4 wall-time partition',
    }
    with (dest / 'stage_manifest.json').open('x') as stream:
        json.dump(manifest, stream, indent=2)
        stream.write('\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
