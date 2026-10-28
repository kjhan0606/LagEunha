"""Promote only a fully audited 54k-site t=.825 checkpoint toward t=10."""

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
    parser.add_argument('predecessor', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    receipt = json.loads((args.predecessor/'audit_t0825.json').read_text())
    if (not all(receipt['checks'].values()) or abs(receipt['time']-.825)>1e-12
            or receipt['n'] != 54312):
        raise ValueError('t=.825 accuracy audit did not pass')
    source = Path(receipt['endpoint'])
    if source.parent.resolve() != args.predecessor.resolve():
        raise ValueError('Receipt checkpoint is outside predecessor run')
    match = re.fullmatch(r'khout\.(\d+)\.dat', source.name)
    if not match or int(match.group(1)) <= 0:
        raise ValueError('Unexpected complete checkpoint name')
    with source.open('rb') as stream:
        magic, n, stride, stress = struct.unpack('<4i', stream.read(16))
        stream.seek(-52, 2)
        snapshot_time = struct.unpack('<d', stream.read(8))[0]
    if ((magic, n, stride, stress) != (-1, 54312, 576, 1) or
            source.stat().st_size != 16+n*stride+52 or
            abs(snapshot_time-.825)>1e-12):
        raise ValueError('Incomplete checkpoint')
    source_hash = sha256(source)
    if 'sha256' in receipt and source_hash != receipt['sha256'].get(source.name):
        raise ValueError('Checkpoint hash changed after audit')
    if args.destination.exists() and any(args.destination.glob('khout.*.dat')):
        raise ValueError('Destination already contains a checkpoint')
    args.destination.mkdir(parents=True, exist_ok=True)
    for name in (source.name, 'eunha2.exe', 'params.dat'):
        original = source if name == source.name else args.predecessor/name
        target = args.destination/name
        if target.exists():
            raise ValueError(f'Destination file already exists: {target}')
        shutil.copy2(original, target)
    inputs = {}
    for name in (source.name, 'eunha2.exe', 'params.dat'):
        inputs[name] = sha256(args.destination/name)
    plan = dict(predecessor_receipt=str(args.predecessor/'audit_t0825.json'),
                predecessor_checkpoint=str(source), predecessor_sha256=source_hash,
                time=.825, step=int(match.group(1)), target_time=10.,
                inputs_sha256=inputs)
    with (args.destination/'restart_plan.json').open('x') as stream:
        json.dump(plan, stream, indent=2)
        stream.write('\n')
    print(json.dumps(plan, indent=2))


if __name__ == '__main__':
    main()
