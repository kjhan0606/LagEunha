"""Stage an isolated CUDA build with the exact-zero projection repair."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS=Path(__file__).resolve().parent
V3=ROOT/'gpu_energy_neutral_native_build_v3_jacobi_choice'
RUN=ROOT/'gpu_energy_neutral_zero_pressure_build_v4'
EXPECTED={
    'source.tar':'820d69d88ec1a5bd86dfc3d6a9ca15be36f8104f7089e406072a317d9914651f',
    'patch_guard.py':'67c74162615589b1c9611bc892ba5460a6821d99fa17d46904a6b0173298bd61',
    'patch_iterations.py':'78a407e4ae65f909d41e485e8cc965071cdcd0e8e272332bd89aa92351d10251',
}


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def main():
    if RUN.exists():raise FileExistsError(RUN)
    sources={'source.tar':V3/'source.tar','patch_guard.py':V3/'patch_guard.py',
             'patch_iterations.py':V3/'patch_iterations.py',
             'patch_zero.py':DOCS/'apply_kepler_zero_pressure_projection_20261027.py',
             'gfs_pressure_neutral_projection_20261027.h':DOCS/'gfs_pressure_neutral_projection_20261027.h',
             'test_projection.cpp':DOCS/'test_gfs_pressure_neutral_projection_20261027.cpp',
             'job.slurm':DOCS/'build_kepler_zero_pressure_projection_20261027.slurm'}
    for name,value in EXPECTED.items():
        if digest(sources[name])!=value:raise ValueError(f'Changed v3 input: {name}')
    subprocess.run(['bash','-n',str(sources['job.slurm'])],check=True)
    RUN.mkdir()
    for name,source in sources.items():
        if name=='source.tar':(RUN/name).symlink_to(source.resolve())
        else:shutil.copy2(source,RUN/name)
    hashes={name:digest(RUN/name) for name in sources}
    if any(hashes[name]!=digest(source) for name,source in sources.items()):
        raise ValueError('Staged build input differs from source')
    with (RUN/'inputs.sha256').open('x') as f:
        for name,value in sorted(hashes.items()):f.write(f'{value}  {name}\n')
    manifest={'run':str(RUN),'inputs':{name:{'source':str(sources[name]),'sha256':hashes[name]}
                                     for name in sorted(hashes)},
              'scope':'Isolated v3 rebuild; exact-zero pressure gradient leaves '
                      'finite rates unchanged, nonzero arithmetic preserved.'}
    with (RUN/'stage_manifest.json').open('x') as f:
        json.dump(manifest,f,indent=2);f.write('\n')
    print(json.dumps({'run':str(RUN),'input_count':len(hashes)},indent=2))


if __name__=='__main__':
    main()
