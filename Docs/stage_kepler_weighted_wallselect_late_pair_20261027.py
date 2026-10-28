"""Conditional v5/v10 comparison from the future accepted weighted t=6.5 state."""

import json
from pathlib import Path
import shutil
import subprocess

from stage_kepler_capacity64_recovery_20261027 import (
    BINARY_SHA, DOCS, INITIAL_SHA, MASK_SHA, ROOT, RUN as RECOVERY, digest,
)

PARENT=ROOT/'gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1'
SPEED=ROOT/'gpu_energy_neutral_wallselect_cache_pair_v1/speed_audit.json'
CAPACITY=ROOT/'gpu_energy_neutral_wallselect_capacity_pair_t006_v1/pair_audit.json'
BUILD=ROOT/'gpu_energy_neutral_wallselect_cache_capacity64_build_v10'
RUN=ROOT/'gpu_energy_neutral_wallselect_weighted_late_pair_v1'


def main():
    candidate=BUILD/'eunha2_wallselect_cache_capacity64_v10.exe'
    candidate_sha=(BUILD/'binary.sha256').read_text().split()[0]
    if digest(candidate)!=candidate_sha:raise ValueError('Candidate binary changed')
    manifest=json.loads((RECOVERY/'stage_manifest.json').read_text())
    if (manifest['binary_sha256']!=BINARY_SHA or manifest['initial_sha256']!=INITIAL_SHA
            or manifest['mask_sha256']!=MASK_SHA):raise ValueError('Recovery inputs changed')
    template=(RECOVERY/'job.slurm').read_text()
    if digest(RECOVERY/'job.slurm')!=manifest['inputs']['job.slurm']:
        raise ValueError('Recovery runner changed')
    header=template[:template.index('sha256sum -c inputs.sha256')]
    header=header.replace('LagQ2W64Recover005v5','LagQ2WallWeightedLate')
    env=template[template.index('export GFS_GPU_VOLW2='):template.index('nvidia-smi -L')]
    env=env.replace('HYDRO_TSTOP=.05 EUNHA_DUMP_DT=.05',
                    'HYDRO_TSTOP=6.525 EUNHA_DUMP_DT=6.525')
    job=header+env+'''root=$PWD
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 - <<'PY_RESTART'
import json
from pathlib import Path
import shutil
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

parent=Path('''+repr(str(PARENT))+''')
receipt=json.loads((parent/'audit_t65.json').read_text())
for path in (Path('''+repr(str(SPEED))+'''), Path('''+repr(str(CAPACITY))+''')):
    gate=json.loads(path.read_text())
    if not all(gate['checks'].values()): raise ValueError(f'Prior pair failed: {path}')
capacity=json.loads(Path('''+repr(str(CAPACITY))+''').read_text())
if (capacity['control_binary_sha256']!='''+repr(BINARY_SHA)+'''
        or capacity['candidate_binary_sha256']!='''+repr(candidate_sha)+'''):
    raise ValueError('Prior capacity pair used different binaries')
speed_manifest=json.loads((Path('''+repr(str(SPEED))+''').parent/'stage_manifest.json').read_text())
if speed_manifest['candidate_sha256']!='''+repr(candidate_sha)+''':
    raise ValueError('Prior speed pair used a different candidate')
if (len(receipt['checks'])!=21 or not all(receipt['checks'].values())
        or receipt['binary_sha256']!='''+repr(BINARY_SHA)+'''
        or receipt['initial_sha256']!='''+repr(INITIAL_SHA)+'''
        or abs(receipt['time']-6.5)>1e-12): raise ValueError('t=6.5 parent gate failed')
source=Path(receipt['endpoint'])
if (source.parent!=parent or digest(source)!=receipt['endpoint_sha256']
        or abs(snapshot_time(source,55351)-6.5)>1e-10): raise ValueError('Parent checkpoint changed')
for label in ('control','candidate'):
    target=Path(label)/'khout.000001.dat'
    if target.exists():raise FileExistsError(target)
    shutil.copy2(source,target)
    if digest(target)!=receipt['endpoint_sha256']:raise ValueError('Restart copy changed')
with Path('restart_origin.json').open('x') as out:
    json.dump({'receipt':str(parent/'audit_t65.json'),'source':str(source),
               'sha256':receipt['endpoint_sha256'],'time':6.5},out,indent=2)
    out.write('\\n')
PY_RESTART
nvidia-smi -L
for variant in control candidate; do
    cd "$root/$variant"
    sha256sum -c inputs.sha256
    export GFS_GPU_VOLW2_MASK="$PWD/disk_active_mask.bin"
    if [[ $variant == control ]]; then expected_binary='''+BINARY_SHA+'''; else expected_binary='''+candidate_sha+'''; fi
    TIMEFORMAT='[RUN_TIME] wall=%R'
    { time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1; } 2> runtime.log
    /home/kjhan/miniconda3/bin/python3 audit_kepler_q2_weighted_equilibrium_20261027.py \\
      initial.dat khout.000001.dat . --parent-time 6.5 --time 6.525 \\
      --expected-binary-sha256 "$expected_binary" --output audit.json > audit.log 2>&1
done
cd "$root"
/home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_fast_pair_20261027.py \\
  control/audit.json candidate/audit.json --control-binary-sha256 '''+BINARY_SHA+''' \\
  --candidate-binary-sha256 '''+candidate_sha+''' --parent-time 6.5 --time 6.525 \\
  --output pair_audit.json > pair_audit.log
'''
    if RUN.exists():raise FileExistsError(RUN)
    RUN.mkdir();(RUN/'job.slurm').write_text(job)
    subprocess.run(['bash','-n',str(RUN/'job.slurm')],check=True)
    root_names=('audit_kepler_weighted_fast_pair_20261027.py',
                'audit_kepler_sparse_area_t0825_20261027.py',
                'prepare_kepler_sparse_area_restart_20261027.py')
    root_hashes={'job.slurm':digest(RUN/'job.slurm')}
    for name in root_names:
        shutil.copy2(DOCS/name,RUN/name);root_hashes[name]=digest(RUN/name)
        if root_hashes[name]!=digest(DOCS/name):raise ValueError(f'Root input changed: {name}')
    with (RUN/'inputs.sha256').open('x') as out:
        for name,value in sorted(root_hashes.items()):out.write(f'{value}  {name}\n')
    names=[x for x in manifest['inputs'] if x not in ('job.slurm','parent_audit.json','eunha2.exe')]
    for label,exe in (('control',RECOVERY/'eunha2.exe'),('candidate',candidate)):
        directory=RUN/label;directory.mkdir();hashes={}
        sources={name:RECOVERY/name for name in names};sources['eunha2.exe']=exe
        for name,source in sources.items():
            shutil.copy2(source,directory/name);hashes[name]=digest(directory/name)
            if hashes[name]!=digest(source):raise ValueError(f'Copy changed: {label}/{name}')
        if hashes['initial.dat']!=INITIAL_SHA or hashes['disk_active_mask.bin']!=MASK_SHA:
            raise ValueError('Initial state or mask changed')
        with (directory/'inputs.sha256').open('x') as out:
            for name,value in sorted(hashes.items()):out.write(f'{value}  {name}\n')
    with (RUN/'stage_manifest.json').open('x') as out:
        json.dump({'parent':str(PARENT),'prior_speed_receipt':str(SPEED),
                   'prior_capacity_receipt':str(CAPACITY),'control_sha256':BINARY_SHA,
                   'candidate_sha256':candidate_sha,'time':6.5,'target':6.525},out,indent=2)
        out.write('\n')
    print(RUN)


if __name__=='__main__':main()
