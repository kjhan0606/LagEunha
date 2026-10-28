"""Stage isolated W0 exact-bin-pruning CUDA build from pinned v3 source."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest
from stage_kepler_tess_work_diag_v8_20261027 import DOCS, PARENT, EXPECTED


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    runner = DOCS / 'build_kepler_tess_seed_bin_prune_v9_lageunha_20261027.sh'
    patcher = DOCS / 'patch_kepler_tess_seed_bin_prune_v9_20261027.py'
    subprocess.run(['bash', '-n', str(runner)], check=True)
    for name, expected in EXPECTED.items():
        if digest(PARENT / name) != expected:
            raise ValueError(f'Parent input changed: {name}')
    dest.mkdir()
    inputs = {}
    for name in EXPECTED:
        (dest / name).symlink_to((PARENT / name).resolve())
        inputs[name] = digest(dest / name)
    for name, source in [('patch_binprune.py', patcher), ('runner.sh', runner)]:
        shutil.copy2(source, dest / name)
        inputs[name] = digest(dest / name)
    with (dest / 'inputs.sha256').open('x') as out:
        for name in sorted(inputs):
            out.write(f'{inputs[name]}  {name}\n')
    with (dest / 'stage_manifest.json').open('x') as out:
        json.dump({'run': str(dest), 'parent': str(PARENT), 'inputs': inputs,
                   'scope': 'Experimental W0-certified 5x5 seed-bin pruning; build only.'},
                  out, indent=2)
        out.write('\n')
    print(json.dumps({'run': str(dest), 'inputs': inputs}, indent=2))


if __name__ == '__main__':
    main()
