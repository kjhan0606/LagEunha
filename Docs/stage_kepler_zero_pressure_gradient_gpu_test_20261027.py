"""Stage old/fixed CUDA exact-uniform-pressure regression test."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS=Path(__file__).resolve().parent
RUN=ROOT/'zero_pressure_gradient_gpu_test_v1'
ARCHIVE=ROOT/'gpu_energy_neutral_native_build_v3_jacobi_choice/source.tar'
ARCHIVE_SHA='820d69d88ec1a5bd86dfc3d6a9ca15be36f8104f7089e406072a317d9914651f'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    if RUN.exists():raise FileExistsError(RUN)
    if digest(ARCHIVE)!=ARCHIVE_SHA:raise ValueError('v3 source archive changed')
    sources={'source.tar':ARCHIVE,
             'job.slurm':DOCS/'kepler_zero_pressure_gradient_gpu_test_20261027.slurm',
             'patch_zero.py':DOCS/'apply_kepler_zero_pressure_projection_20261027.py',
             'gfs_pressure_neutral_projection_20261027.h':DOCS/'gfs_pressure_neutral_projection_20261027.h',
             'test_kepler_zero_pressure_gradient_cuda_20261027.cu':DOCS/'test_kepler_zero_pressure_gradient_cuda_20261027.cu'}
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir()
    for name,source in sources.items():
        if name=='source.tar':(RUN/name).symlink_to(source.resolve())
        else:shutil.copy2(source,RUN/name)
        if digest(RUN/name)!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as f:
        for name in sorted(sources):f.write(f'{digest(RUN/name)}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as f:
        json.dump({'run':str(RUN),'inputs':{name:digest(RUN/name) for name in sources},
                   'scope':'Synthetic two-cell exact-uniform-pressure test of old and fixed '
                           'native CUDA weight-rate probe; not a fluid simulation.'},f,indent=2)
        f.write('\n')
    print(json.dumps({'run':str(RUN),'input_count':len(sources)}))


if __name__=='__main__':main()
