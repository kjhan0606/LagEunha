"""Stage CPU-only kinematic audit of two preserved W0 snapshots."""

import hashlib
from pathlib import Path
import shutil
import subprocess

DOCS = Path(__file__).resolve().parent
RUN = Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/w0_centroid_kinematics_v1')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    sources = {
        'job.slurm': DOCS/'kepler_w0_centroid_kinematics_20261027.slurm',
        'aligned_receipt.json': DOCS/'KEPLER_BOUNDED_FAST_RK4_2026-10-27/audit_t20.json',
        'jitter_receipt.json': DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json',
    }
    for name in ('audit_kepler_w0_centroid_kinematics_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        sources[name] = DOCS/name
    if RUN.exists():
        raise FileExistsError(RUN)
    subprocess.run(['bash', '-n', str(sources['job.slurm'])], check=True)
    RUN.mkdir()
    hashes = {}
    for name, source in sources.items():
        shutil.copy2(source, RUN/name)
        hashes[name] = digest(RUN/name)
        if hashes[name] != digest(source):
            raise ValueError(f'Input copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name, value in sorted(hashes.items()):
            out.write(f'{value}  {name}\n')
    print(RUN)


if __name__ == '__main__':
    main()
