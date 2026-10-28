"""Stage exact W0 seed-bin pruning for matched Ada W0 and weighted replays."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time
from stage_kepler_q2_seed3_matched_t6501_20261027 import (
    DOCS, ROOT, PARENT, PARENT_RECEIPT, HELPERS, copy_checked,
)


BINARY = (ROOT / 'gpu_tess_seed_bin_prune_v9_build' /
          'eunha2_tess_binprune_v9.exe')
BINARY_SHA = '3f18131a31c235703965a094043747d9527fdaa3973d4f71a0c6e7618e465553'
BASE_RUNNER = DOCS / 'kepler_q2_tess_triple_ada_t6501_20261027.sh'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    parent = json.loads(PARENT_RECEIPT.read_text())
    checkpoint = Path(parent['endpoint'])
    if (parent['time'] != 6.5 or not all(parent['checks'].values())
            or checkpoint.parent != PARENT
            or digest(checkpoint) != parent['endpoint_sha256']
            or abs(snapshot_time(checkpoint, 55351)-6.5) > 1e-10):
        raise ValueError('Accepted t=6.5 parent changed')
    if digest(PARENT / 'params.dat') != '2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43':
        raise ValueError('Physical parameters changed')
    if digest(BINARY) != BINARY_SHA:
        raise ValueError('v9 CUDA binary changed')
    mask_receipt = json.loads((DOCS / 'KEPLER_Q2_BOUNDED_NATIVE_WDOT_PROBE_INPUT_T65_2026-10-27.json').read_text())
    mask = Path(mask_receipt['mask'])
    if mask.stat().st_size != 55351 or digest(mask) != mask_receipt['mask_sha256']:
        raise ValueError('Disk weight mask changed')
    text = BASE_RUNNER.read_text()
    old = 'for version in v3 v4 v6; do'
    if text.count(old) != 1:
        raise ValueError('Unexpected baseline runner loop')
    text = text.replace(old, 'for version in v9; do')
    dest.mkdir()
    common = dest / 'common'
    common.mkdir()
    shutil.copy2(mask, common / 'disk_active_mask.bin')
    (common / 'runner.sh').write_text(text)
    subprocess.run(['bash', '-n', str(common / 'runner.sh')], check=True)
    common_hashes = {p.name: digest(p) for p in
                     (common / 'disk_active_mask.bin', common / 'runner.sh')}
    with (common / 'inputs.sha256').open('x') as out:
        for name, value in sorted(common_hashes.items()):
            out.write(f'{value}  {name}\n')
    cases = {}
    for mode in ('w0', 'weighted'):
        cases[mode] = copy_checked(dest / f'{mode}_v9', {
            'eunha2.exe': BINARY,
            'params.dat': PARENT / 'params.dat',
            'khout.000001.dat': checkpoint,
            **{helper: DOCS / helper for helper in HELPERS},
        })
    with (dest / 'stage_manifest.json').open('x') as out:
        json.dump({'run': str(dest), 'host': 'lageunha', 'gpu': 'RTX 5000 Ada',
                   'time_start': 6.5, 'time_end': 6.501,
                   'parent_sha256': parent['endpoint_sha256'],
                   'binary_sha256': BINARY_SHA,
                   'runner_base_sha256': digest(BASE_RUNNER),
                   'common_inputs': common_hashes, 'cases': cases,
                   'scope': 'Native one-step W0/weighted gate for exact W0 seed-bin pruning; no long-run adoption.'},
                  out, indent=2)
        out.write('\n')
    print(json.dumps({'run': str(dest), 'binary_sha256': BINARY_SHA,
                      'cases': list(cases)}, indent=2))


if __name__ == '__main__':
    main()
