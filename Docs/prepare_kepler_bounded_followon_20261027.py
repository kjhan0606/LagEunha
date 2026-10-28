"""Stage a bounded Kepler restart only after its parent native audit passes."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil

from prepare_kepler_sparse_area_restart_20261027 import snapshot_time


EXPECTED_BINARY = '3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e'
COPY_FILES = (
    'eunha2.exe', 'params.dat', 'glass.bin', 'areas.bin',
    'audit_kepler_bounded_native_evolution_20261027.py',
    'audit_kepler_sparse_area_t0825_20261027.py',
    'audit_kepler_sparse_native_ballistic_20261027.py',
    'audit_kepler_sparse_state_consistency_20261027.py',
    'prepare_kepler_sparse_area_restart_20261027.py',
)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('parent', type=Path)
    ap.add_argument('parent_audit')
    ap.add_argument('parent_time', type=float)
    ap.add_argument('runner')
    args = ap.parse_args()
    parent = args.parent.resolve()
    dest = Path.cwd().resolve()
    if (dest / 'khout.000001.dat').exists() or (dest / 'inputs.sha256').exists():
        raise FileExistsError('Follow-on stage already has restart input')
    audit_path = parent / args.parent_audit
    audit = json.loads(audit_path.read_text())
    if abs(audit['time'] - args.parent_time) > 1e-10 or not all(audit['checks'].values()):
        raise ValueError('Parent native physics audit did not accept the requested time')
    endpoint = Path(audit['endpoint']).resolve()
    if endpoint.parent != parent or not endpoint.is_file():
        raise ValueError('Parent endpoint does not belong to the parent run')
    if digest(endpoint) != audit['endpoint_sha256']:
        raise ValueError('Parent endpoint SHA differs from audit')
    if abs(snapshot_time(endpoint, 55351) - args.parent_time) > 1e-10:
        raise ValueError('Parent endpoint time is wrong or incomplete')
    initial = Path(audit['initial']).resolve()
    if digest(initial) != audit['initial_sha256']:
        raise ValueError('Initial state SHA differs from parent audit')
    if digest(parent / 'eunha2.exe') != EXPECTED_BINARY:
        raise ValueError('Parent binary differs from the accepted bounded binary')
    if not (dest / args.runner).is_file():
        raise FileNotFoundError('Submitted runner is absent from follow-on stage')
    for name in COPY_FILES:
        if (dest / name).exists():
            if name != 'prepare_kepler_sparse_area_restart_20261027.py' or digest(dest / name) != digest(parent / name):
                raise FileExistsError(f'Refusing to overwrite {dest / name}')
        else:
            shutil.copy2(parent / name, dest / name)
    shutil.copy2(endpoint, dest / 'khout.000001.dat')
    files = (*COPY_FILES, 'khout.000001.dat', args.runner,
             'prepare_kepler_bounded_followon_20261027.py')
    with (dest / 'inputs.sha256').open('x') as stream:
        for name in files:
            stream.write(f'{digest(dest / name)}  {name}\n')
    manifest = {
        'parent': str(parent), 'parent_audit': str(audit_path),
        'parent_audit_sha256': digest(audit_path),
        'accepted_parent_time': args.parent_time,
        'parent_endpoint': str(endpoint),
        'parent_endpoint_sha256': audit['endpoint_sha256'],
        'initial': str(initial), 'initial_sha256': audit['initial_sha256'],
        'binary_sha256': EXPECTED_BINARY,
        'runner': args.runner,
    }
    with (dest / 'restart_origin.json').open('x') as stream:
        json.dump(manifest, stream, indent=2)
        stream.write('\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
