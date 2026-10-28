"""Stage isolated, hash-guarded CUDA work-counter builds on lageunha."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
PARENT = ROOT / 'gpu_energy_neutral_native_build_v3_jacobi_choice'
EXPECTED = {
    'source.tar': '820d69d88ec1a5bd86dfc3d6a9ca15be36f8104f7089e406072a317d9914651f',
    'patch_guard.py': '67c74162615589b1c9611bc892ba5460a6821d99fa17d46904a6b0173298bd61',
    'patch_iterations.py': '78a407e4ae65f909d41e485e8cc965071cdcd0e8e272332bd89aa92351d10251',
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('destination', type=Path)
    args = ap.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    runner = DOCS / 'build_kepler_tess_work_diag_v8_lageunha_20261027.sh'
    patcher = DOCS / 'patch_kepler_tess_work_diag_v8_20261027.py'
    subprocess.run(['bash', '-n', str(runner)], check=True)
    for name, sha in EXPECTED.items():
        if digest(PARENT / name) != sha:
            raise ValueError(f'Changed parent input: {name}')
    dest.mkdir()
    inputs = {}
    for name in EXPECTED:
        (dest / name).symlink_to((PARENT / name).resolve())
        inputs[name] = digest(dest / name)
    for name, source in [('patch_diag.py', patcher), ('runner.sh', runner)]:
        shutil.copy2(source, dest / name)
        inputs[name] = digest(dest / name)
    with (dest / 'inputs.sha256').open('x') as out:
        for name in sorted(inputs):
            out.write(f'{inputs[name]}  {name}\n')
    with (dest / 'stage_manifest.json').open('x') as out:
        json.dump({'run': str(dest), 'parent': str(PARENT), 'inputs': inputs,
                   'purpose': 'Diagnostic-only CUDA seed and clipping work counts; not a production binary'},
                  out, indent=2)
        out.write('\n')
    print(json.dumps({'run': str(dest), 'inputs': inputs}, indent=2))


if __name__ == '__main__':
    main()
