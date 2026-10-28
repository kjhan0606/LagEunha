"""Stage same-state v3/v4/v6 one-step CUDA comparisons on the idle Ada GPU."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time
from stage_kepler_q2_seed3_matched_t6501_20261027 import (
    DOCS, ROOT, PARENT, PARENT_RECEIPT, HELPERS, copy_checked,
)


BUILDS = {
    'v3': (ROOT / 'gpu_energy_neutral_native_build_v3_jacobi_choice' /
           'eunha2_energy_neutral_v3.exe',
           '7358c54996b820a7b45bd53bcf6f788f125754b3158845d2b4265afead570f38'),
    'v4': (ROOT / 'gpu_energy_neutral_native_build_v4_seed3' /
           'eunha2_energy_neutral_seed3_v4.exe',
           '945864b0c7065722d464420edaa92505665d94101caf6ef7cee666aff92aa989'),
    'v6': (ROOT / 'gpu_energy_neutral_native_build_v6_raw_limit1' /
           'eunha2_energy_neutral_raw_limit1_v6.exe',
           '91e55fb83584b8a0a246d44b7e2bfd20667712cd4b08dc22ba7643e309108519'),
}
RUNNER = 'kepler_q2_tess_triple_ada_t6501_20261027.sh'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(['bash', '-n', str(DOCS / RUNNER)], check=True)
    parent = json.loads(PARENT_RECEIPT.read_text())
    checkpoint = Path(parent['endpoint'])
    if (parent['time'] != 6.5 or not all(parent['checks'].values())
            or checkpoint.parent != PARENT
            or digest(checkpoint) != parent['endpoint_sha256']
            or abs(snapshot_time(checkpoint, 55351) - 6.5) > 1e-10):
        raise ValueError('Accepted jitter checkpoint changed')
    if digest(PARENT / 'params.dat') != '2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43':
        raise ValueError('Physical parameters changed')
    mask_receipt = json.loads((DOCS / 'KEPLER_Q2_BOUNDED_NATIVE_WDOT_PROBE_INPUT_T65_2026-10-27.json').read_text())
    mask = Path(mask_receipt['mask'])
    if mask.stat().st_size != 55351 or digest(mask) != mask_receipt['mask_sha256']:
        raise ValueError('Disk weight mask changed')
    for path, expected in BUILDS.values():
        if digest(path) != expected:
            raise ValueError(f'Native binary changed: {path}')
    dest.mkdir()
    common = copy_checked(dest / 'common', {
        'disk_active_mask.bin': mask,
        'runner.sh': DOCS / RUNNER,
    })
    cases = {}
    for mode in ('w0', 'weighted'):
        for label, (binary, _) in BUILDS.items():
            name = f'{mode}_{label}'
            cases[name] = copy_checked(dest / name, {
                'eunha2.exe': binary,
                'params.dat': PARENT / 'params.dat',
                'khout.000001.dat': checkpoint,
                **{helper: DOCS / helper for helper in HELPERS},
            })
    with (dest / 'stage_manifest.json').open('x') as out:
        json.dump({'run': str(dest), 'host': 'lageunha', 'gpu': 'RTX 5000 Ada',
                   'time_start': 6.5, 'time_end': 6.501,
                   'parent_receipt': str(PARENT_RECEIPT),
                   'parent_sha256': parent['endpoint_sha256'],
                   'binary_sha256': {k: v[1] for k, v in BUILDS.items()},
                   'common_inputs': common, 'cases': cases,
                   'scope': 'One-step same-state W0/weighted native comparison; not long-time validation.'},
                  out, indent=2)
        out.write('\n')
    print(json.dumps({'run': str(dest), 'cases': list(cases)}, indent=2))


if __name__ == '__main__':
    main()
