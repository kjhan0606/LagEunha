"""Stage CPU face-slip lower-bound audits from accepted native states."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS=Path(__file__).resolve().parent
RUN=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/optimal_material_face_rate_t65_cpu_v1')


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    sources={
        'job.slurm':DOCS/'kepler_optimal_face_material_rate_t65_cpu_20261027.slurm',
        'aligned_receipt.json':DOCS/'KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json',
        'jitter_receipt.json':DOCS/'KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json',
    }
    for name in ('audit_kepler_optimal_face_material_rate_20261027.py',
                 'audit_kepler_bounded_pressure_force_t10_20261027.py',
                 'audit_kepler_sparse_area_t0825_20261027.py',
                 'laguerre_nonlinear_volume_probe_t6_20261026.py',
                 'prepare_kepler_sparse_area_restart_20261027.py'):
        sources[name]=DOCS/name
    if RUN.exists():
        raise FileExistsError(RUN)
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir()
    hashes={}
    for name,source in sources.items():
        shutil.copy2(source,RUN/name)
        hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):
            raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as f:
        for name,value in sorted(hashes.items()):f.write(f'{value}  {name}\n')
    manifest={'run':str(RUN),'inputs':{name:{'source':str(sources[name]),'sha256':hashes[name]}
                                     for name in sorted(hashes)},
              'scope':'CPU-only, read-only native checkpoint geometry and unconstrained '
                      'least-squares disk-face speed lower bound.'}
    with (RUN/'stage_manifest.json').open('x') as f:
        json.dump(manifest,f,indent=2);f.write('\n')
    print(json.dumps({'run':str(RUN),'input_count':len(hashes)},indent=2))


if __name__=='__main__':
    main()
