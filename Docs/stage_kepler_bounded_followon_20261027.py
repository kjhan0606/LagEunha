"""Prestage a dependent bounded Kepler job without copying unaudited state."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
VARIANTS = {
    'aligned_to10': {
        'runner': 'kepler_q2_bounded_poly1024_to10_a10_20261027.slurm',
        'parent_job': 437317,
    },
    'jitter_to65': {
        'runner': 'kepler_q2_bounded_jitter005_to65_a10_20261027.slurm',
        'parent_job': 437320,
    },
}
HELPERS = (
    'prepare_kepler_bounded_followon_20261027.py',
    'prepare_kepler_sparse_area_restart_20261027.py',
)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('variant', choices=VARIANTS)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    config = VARIANTS[args.variant]
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(f'Refusing existing stage: {dest}')
    subprocess.run(['bash', '-n', str(ROOT / config['runner'])], check=True)
    dest.mkdir(parents=True)
    files = (*HELPERS, config['runner'])
    for name in files:
        shutil.copy2(ROOT / name, dest / name)
    manifest = {
        'variant': args.variant,
        'destination': str(dest),
        'parent_job': config['parent_job'],
        'dependency': f'afterok:{config["parent_job"]}',
        'runner': config['runner'],
        'prestaged_sha256': {name: digest(dest / name) for name in files},
        'physics_state': 'Copied only at job start after parent native audit passes',
    }
    with (dest / 'prestage_manifest.json').open('x') as stream:
        json.dump(manifest, stream, indent=2)
        stream.write('\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
