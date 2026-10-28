"""Try only the independent Kepler-acceleration cache on weighted v5."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS = Path(__file__).resolve().parent
V6 = ROOT/'gpu_energy_neutral_fast_capacity64_build_v6'
RUN = ROOT/'gpu_energy_neutral_accel_cache_capacity64_build_v7'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parent = json.loads((V6/'stage_manifest.json').read_text())
    template = (V6/'job.slurm').read_text()
    if digest(V6/'job.slurm') != parent['inputs']['job.slurm']:
        raise ValueError('v6 build runner changed')
    for old, new in {
        'LagQ2W64FastBuild': 'LagQ2W64CacheBuild',
        'lag_q2_fast_capacity64_': 'lag_q2_cache_capacity64_',
        'for patch_name in gap.patch accel.patch diagnostic.patch thermal.patch; do':
            'for patch_name in accel.patch; do',
        'eunha2_fast_capacity64_v6.exe': 'eunha2_cache_capacity64_v7.exe',
    }.items():
        expected_count = 2 if old == 'eunha2_fast_capacity64_v6.exe' else 1
        if template.count(old) != expected_count:
            raise ValueError(f'Unexpected v6 runner: {old!r}')
        template = template.replace(old, new)
    names = ('source.tar', 'patch_guard.py', 'patch_iterations.py', 'patch_zero.py',
             'patch_capacity.py', 'gfs_pressure_neutral_projection_20261027.h',
             'test_projection.cpp', 'accel.patch')
    for name in names:
        if digest(V6/name) != parent['inputs'][name]:
            raise ValueError(f'v6 input changed: {name}')
    if RUN.exists():
        raise FileExistsError(RUN)
    RUN.mkdir()
    (RUN/'job.slurm').write_text(template)
    subprocess.run(['bash', '-n', str(RUN/'job.slurm')], check=True)
    hashes = {'job.slurm': digest(RUN/'job.slurm')}
    for name in names:
        if name == 'source.tar':
            (RUN/name).symlink_to((V6/name).resolve())
        else:
            shutil.copy2(V6/name, RUN/name)
        hashes[name] = digest(RUN/name)
        if hashes[name] != parent['inputs'][name]:
            raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'parent': str(V6), 'inputs': hashes,
                   'scope': 'Experimental isolated acceleration-cache patch only; '
                            'native weighted paired validation is required.'}, out, indent=2)
        out.write('\n')
    print(RUN)


if __name__ == '__main__':
    main()
