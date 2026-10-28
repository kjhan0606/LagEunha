"""Stage a same-checkpoint H200 tessellation speed and physics comparison."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT = ROOT / 'jitter005_to65_a10_v1'
PARENT_RECEIPT = DOCS / 'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json'
BUILDS = {
    'v3': (ROOT / 'gpu_energy_neutral_native_build_v3_jacobi_choice' /
           'eunha2_energy_neutral_v3.exe',
           '7358c54996b820a7b45bd53bcf6f788f125754b3158845d2b4265afead570f38'),
    'v4': (ROOT / 'gpu_energy_neutral_native_build_v4_seed3' /
           'eunha2_energy_neutral_seed3_v4.exe',
           '945864b0c7065722d464420edaa92505665d94101caf6ef7cee666aff92aa989'),
}
HELPERS = (
    'audit_kepler_bounded_short_replay_20261027.py',
    'audit_kepler_q2_weighted_short_replay_20261027.py',
    'audit_kepler_sparse_area_t0825_20261027.py',
    'audit_kepler_sparse_state_consistency_20261027.py',
    'laguerre_nonlinear_volume_probe_t6_20261026.py',
    'prepare_kepler_sparse_area_restart_20261027.py',
)
RUNNER = 'kepler_q2_seed3_matched_t6501_h200_20261027.slurm'


def copy_checked(directory, sources):
    directory.mkdir()
    hashes = {}
    for name, source in sources.items():
        shutil.copy2(source, directory / name)
        hashes[name] = digest(directory / name)
    with (directory / 'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    return hashes


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
        raise ValueError('The accepted jitter t=6.5 checkpoint changed')
    if digest(PARENT / 'params.dat') != '2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43':
        raise ValueError('Physical parameter file changed')
    mask_receipt = json.loads((DOCS / 'KEPLER_Q2_BOUNDED_NATIVE_WDOT_PROBE_INPUT_T65_2026-10-27.json').read_text())
    mask = Path(mask_receipt['mask'])
    if mask.stat().st_size != 55351 or digest(mask) != mask_receipt['mask_sha256']:
        raise ValueError('Weighted disk mask changed')
    for binary, expected in BUILDS.values():
        if digest(binary) != expected:
            raise ValueError(f'Built binary changed: {binary}')
    dest.mkdir()
    root_hashes = copy_checked(dest / 'common', {
        'disk_active_mask.bin': mask,
        'job.slurm': DOCS / RUNNER,
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
    manifest = {
        'parent_receipt': str(PARENT_RECEIPT),
        'parent_sha256': parent['endpoint_sha256'],
        'time_start': 6.5,
        'time_end': 6.501,
        'build_binary_sha256': {k: v[1] for k, v in BUILDS.items()},
        'common_inputs': root_hashes,
        'case_inputs': cases,
        'scope': 'Same jitter checkpoint and H200 GPU; v3 seed5 versus v4 seed3, W0 and weighted16, one RK4 step.',
    }
    with (dest / 'stage_manifest.json').open('x') as out:
        json.dump(manifest, out, indent=2)
        out.write('\n')
    print(json.dumps({'destination': str(dest), 'cases': list(cases),
                      'parent_sha256': parent['endpoint_sha256']}, indent=2))


if __name__ == '__main__':
    main()
