"""Stage isolated v5 plus previously validated W0 CPU performance patches."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS = Path(__file__).resolve().parent
V5 = ROOT/'gpu_energy_neutral_capacity64_build_v5'
RUN = ROOT/'gpu_energy_neutral_fast_capacity64_build_v6'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parent = json.loads((V5/'stage_manifest.json').read_text())
    sources = {name: V5/name for name in ('source.tar', 'patch_guard.py',
               'patch_iterations.py', 'patch_zero.py', 'patch_capacity.py',
               'gfs_pressure_neutral_projection_20261027.h', 'test_projection.cpp')}
    sources.update({'gap.patch': DOCS/'KEPLER_RK4_GAP_PROFILE_2026-10-27.patch',
                    'accel.patch': DOCS/'KEPLER_ACCEL_CONFIG_CACHE_2026-10-27.patch',
                    'diagnostic.patch': DOCS/'KEPLER_RK4_STAGE_DIAG_TIMER_2026-10-27.patch',
                    'thermal.patch': DOCS/'KEPLER_RK4_DIAG_THERMAL_HOIST_2026-10-27.patch',
                    'job.slurm': DOCS/'build_kepler_weighted_fast_capacity64_v6_20261027.slurm'})
    for name in list(sources)[:7]:
        if digest(sources[name]) != parent['inputs'][name]['sha256']:
            raise ValueError(f'v5 build input changed: {name}')
    if RUN.exists():
        raise FileExistsError(RUN)
    subprocess.run(['bash', '-n', str(sources['job.slurm'])], check=True)
    RUN.mkdir()
    hashes = {}
    for name, source in sources.items():
        if name == 'source.tar':
            (RUN/name).symlink_to(source.resolve())
        else:
            shutil.copy2(source, RUN/name)
        hashes[name] = digest(RUN/name)
        if hashes[name] != digest(source):
            raise ValueError(f'Staged input changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'v5': str(V5), 'inputs': hashes, 'purpose': 'Experimental '
                   'weighted build using previously bitwise-tested W0 host optimizations; '
                   'must pass a weighted paired native comparison before use.'}, out, indent=2)
        out.write('\n')
    print(RUN)


if __name__ == '__main__':
    main()
