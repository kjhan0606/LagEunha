"""Promote only a fully audited 54k-site profile IC t=.001 checkpoint."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import struct


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('preflight', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    receipt = json.loads((args.preflight/'audit_t001.json').read_text())
    if (not all(receipt['checks'].values()) or abs(receipt['t1']-.001)>1e-12
            or receipt['n'] != 54312):
        raise ValueError('Preflight audit did not pass')
    source = Path(receipt['evolved'])
    if source.parent.resolve() != args.preflight.resolve():
        raise ValueError('Receipt checkpoint is outside preflight run')
    match = re.fullmatch(r'khout\.(\d+)\.dat', source.name)
    if not match or int(match.group(1)) <= 0:
        raise ValueError('Unexpected complete checkpoint name')
    with source.open('rb') as stream:
        magic, n, stride, stress = struct.unpack('<4i', stream.read(16))
    if ((magic, n, stride, stress) != (-1, 54312, 576, 1) or
            source.stat().st_size != 16+n*stride+52):
        raise ValueError('Incomplete checkpoint')
    if sha256(source) != receipt['sha256'][source.name]:
        raise ValueError('Checkpoint hash changed after audit')
    if args.destination.exists() and any(args.destination.glob('khout.*.dat')):
        raise ValueError('Destination already contains a checkpoint')
    args.destination.mkdir(parents=True, exist_ok=True)
    for name in (source.name, 'eunha2.exe', 'params.dat'):
        original = source if name == source.name else args.preflight/name
        target = args.destination/name
        if target.exists():
            raise ValueError(f'Destination file already exists: {target}')
        shutil.copy2(original, target)
    inputs = {}
    for name in (source.name, 'eunha2.exe', 'params.dat'):
        inputs[name] = sha256(args.destination/name)
    plan = dict(preflight_receipt=str(args.preflight/'audit_t001.json'),
                preflight_checkpoint=str(source), time=.001,
                step=int(match.group(1)), target_time=.825,
                inputs_sha256=inputs)
    with (args.destination/'restart_plan.json').open('x') as stream:
        json.dump(plan, stream, indent=2)
        stream.write('\n')
    print(json.dumps(plan, indent=2))


if __name__ == '__main__':
    main()
