"""Try the validated W0 wall-select GPU path on the isolated weighted build."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027')
DOCS=Path(__file__).resolve().parent
V7=ROOT/'gpu_energy_neutral_accel_cache_capacity64_build_v7'
RUN=ROOT/'gpu_energy_neutral_wallselect_cache_capacity64_build_v10'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    parent=json.loads((V7/'stage_manifest.json').read_text())
    template=(V7/'job.slurm').read_text()
    if digest(V7/'job.slurm')!=parent['inputs']['job.slurm']:raise ValueError('v7 runner changed')
    replacements={
        'LagQ2W64CacheBuild':'LagQ2W64WallBuild',
        'lag_q2_cache_capacity64_':'lag_q2_wall_capacity64_',
        'for patch_name in accel.patch; do':'for patch_name in accel.patch wallselect.patch; do',
        'eunha2_cache_capacity64_v7.exe':'eunha2_wallselect_cache_capacity64_v10.exe',
    }
    for old,new in replacements.items():
        count=2 if old=='eunha2_cache_capacity64_v7.exe' else 1
        if template.count(old)!=count:raise ValueError(f'Unexpected v7 anchor: {old!r}')
        template=template.replace(old,new)
    names=('source.tar','patch_guard.py','patch_iterations.py','patch_zero.py',
           'patch_capacity.py','gfs_pressure_neutral_projection_20261027.h',
           'test_projection.cpp','accel.patch')
    sources={name:V7/name for name in names}
    sources['wallselect.patch']=DOCS/'KEPLER_FIXEDPOLY_WALL_SELECT_PRECLIP_2026-10-27.patch'
    for name in names:
        if digest(sources[name])!=parent['inputs'][name]:raise ValueError(f'Changed v7 input: {name}')
    if RUN.exists():raise FileExistsError(RUN)
    RUN.mkdir();(RUN/'job.slurm').write_text(template)
    subprocess.run(['bash','-n',str(RUN/'job.slurm')],check=True)
    hashes={'job.slurm':digest(RUN/'job.slurm')}
    for name,source in sources.items():
        if name=='source.tar':(RUN/name).symlink_to(source.resolve())
        else:shutil.copy2(source,RUN/name)
        hashes[name]=digest(RUN/name)
        if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'parent':str(V7),'inputs':hashes,
                   'scope':'W0-tested host cache and wall-select changes applied '
                           'to weighted CUDA; native weighted paired gate required.'},out,indent=2)
        out.write('\n')
    print(RUN)


if __name__=='__main__':main()
