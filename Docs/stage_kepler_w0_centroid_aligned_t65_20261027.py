"""Stage a same-time aligned control for the W0 centroid diagnostic."""

import hashlib
from pathlib import Path
import shutil
import subprocess

DOCS = Path(__file__).resolve().parent
PARENT = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/w0_centroid_kinematics_v1')
RUN = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/w0_centroid_aligned_t65_v1')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    sources = {'job.slurm': DOCS/'kepler_w0_centroid_aligned_t65_20261027.slurm',
               'aligned_receipt.json': DOCS/'KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json'}
    for name in ('audit_kepler_w0_centroid_kinematics_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        if digest(PARENT/name) != digest(DOCS/name):
            raise ValueError(f'Kinematic audit input changed since previous run: {name}')
        sources[name] = PARENT/name
    if RUN.exists():
        raise FileExistsError(RUN)
    subprocess.run(['bash', '-n', str(sources['job.slurm'])], check=True)
    RUN.mkdir()
    hashes = {}
    for name, source in sources.items():
        shutil.copy2(source, RUN/name)
        hashes[name] = digest(RUN/name)
        if hashes[name] != digest(source):
            raise ValueError(f'Input changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    print(RUN)


if __name__ == '__main__':
    main()
