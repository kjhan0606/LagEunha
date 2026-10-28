"""Stage same Ada state for seed5/seed3 diagnostic tessellation work counts."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time
from stage_kepler_q2_seed3_matched_t6501_20261027 import (
    DOCS, ROOT, PARENT, PARENT_RECEIPT, HELPERS, copy_checked,
)


BUILD = ROOT / 'gpu_tess_work_diag_v8_build'
EXPECTED = {
    'seed5': '46a7cb7a45fb0bf4dcfe27cf5b8108b79ea7b4e570a002bd453e756adfa2cdba',
    'seed3': '57e9de574b4096cede88785bbe38dc057a50a0bf05f5b3a79f62d4bbd8bc534b',
}
RUNNER = DOCS / 'kepler_q2_tess_work_diag_v8_ada_t6501_20261027.sh'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(['bash', '-n', str(RUNNER)], check=True)
    parent = json.loads(PARENT_RECEIPT.read_text())
    checkpoint = Path(parent['endpoint'])
    if (checkpoint.parent != PARENT or parent['time'] != 6.5
            or not all(parent['checks'].values())
            or digest(checkpoint) != parent['endpoint_sha256']
            or abs(snapshot_time(checkpoint, 55351)-6.5) > 1e-10):
        raise ValueError('Accepted t=6.5 parent changed')
    if digest(PARENT / 'params.dat') != '2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43':
        raise ValueError('Physical parameters changed')
    for mode, expected in EXPECTED.items():
        binary = BUILD / f'eunha2_tess_work_{mode}.exe'
        if digest(binary) != expected:
            raise ValueError(f'Diagnostic binary changed: {mode}')
    dest.mkdir()
    shutil.copy2(RUNNER, dest / 'runner.sh')
    with (dest / 'inputs.sha256').open('x') as out:
        out.write(f'{digest(dest / "runner.sh")}  runner.sh\n')
    cases = {}
    for mode in EXPECTED:
        cases[mode] = copy_checked(dest / mode, {
            'eunha2.exe': BUILD / f'eunha2_tess_work_{mode}.exe',
            'params.dat': PARENT / 'params.dat',
            'khout.000001.dat': checkpoint,
            **{helper: DOCS / helper for helper in HELPERS},
        })
    with (dest / 'stage_manifest.json').open('x') as out:
        json.dump({'run': str(dest), 'parent_sha256': parent['endpoint_sha256'],
                   'time_start': 6.5, 'time_end': 6.501,
                   'binary_sha256': EXPECTED, 'case_inputs': cases,
                   'scope': 'Diagnostic one-step W0 replay; work counts are observational and binaries are not production.'},
                  out, indent=2)
        out.write('\n')
    print(json.dumps({'run': str(dest), 'cases': list(cases)}, indent=2))


if __name__ == '__main__':
    main()
