"""Prepare three one-RHS probes with the unchanged bounded native binary.

Outside/on-wall probes are expected to stop at the old post-stage guard.
Their first K1 RHS is diagnostic only, never an accepted evolved state.
"""
from pathlib import Path
import hashlib
import json
import math
import re
import shutil
import struct
import subprocess
import sys

import numpy as np


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parent=Path('/gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029/control_capture')
    initial=parent/'khout.011580.dat'
    if (parent/'exit_status.txt').read_text().strip()!='98':
        raise ValueError('Expected original crossing guard terminal status')
    with initial.open('rb') as f:
        magic,n,stride,stress=struct.unpack('<4i',f.read(16))
        f.seek(-52,2);t=struct.unpack('<d',f.read(8))[0]
    if (magic,n,stride,stress)!=(-1,55351,576,1) or initial.stat().st_size!=16+n*stride+52:
        raise ValueError('Incomplete accepted predecessor')
    reader=Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
    layout=json.loads((reader/'snapshot_offsets.json').read_text())['fields']
    d=np.dtype({'names':['id'],'formats':['<u8'],'offsets':[layout['u4if']['offset']],'itemsize':stride})
    ids=np.memmap(initial,mode='r',dtype=d,offset=16,shape=(n,))['id'] & np.uint64((1<<56)-1)
    matches=np.flatnonzero(ids==55283)
    if len(matches)!=1:raise ValueError('Missing crossing generator')
    i=int(matches[0]);side=543;theta=(side+.5)*2*math.pi/1024
    nx,ny=math.cos(theta),math.sin(theta);apothem=2.75*math.cos(math.pi/1024)
    root=parent.parent/'rhs_extension_probes'
    root.mkdir()
    source_script=(parent/'run_ada.sh').read_text()
    exports='\n'.join(line for line in source_script.splitlines() if line.startswith(('export ','unset ')))
    exports=re.sub(r'HYDRO_TSTOP=48\.4',f'HYDRO_TSTOP={t+1.e-9:.17g}',exports)
    exports='\n'.join(line for line in exports.splitlines() if not line.startswith('export GFS_DIAG_CHECKPOINT_STEP='))
    cases=[]
    for name,delta in [('inside',-1.e-7),('on',0.),('outside',1.e-7)]:
        r=root/name;r.mkdir()
        (r/'eunha2.exe').symlink_to((parent/'eunha2.exe').resolve())
        shutil.copy2(parent/'params.dat',r/'params.dat')
        snapshot=r/'khout.011580.dat';shutil.copy2(initial,snapshot)
        offsets=[16+i*stride+layout[k]['offset'] for k in ['x','y']]
        with snapshot.open('r+b') as f:
            vals=[]
            for off in offsets:f.seek(off);vals.append(struct.unpack('<d',f.read(8))[0])
            old=nx*(vals[0]-3)+ny*(vals[1]-3)-apothem
            shifted=[vals[0]+(delta-old)*nx,vals[1]+(delta-old)*ny]
            for off,value in zip(offsets,shifted):f.seek(off);f.write(struct.pack('<d',value))
        script=f'''#!/bin/bash
set -euo pipefail
cd {r}
{exports}
unset GFS_RK4_FAILURE_SNAPSHOT GFS_DIAG_CHECKPOINT_STEP
export GFS_RK4_AUDIT_STEP_MIN=0 GFS_RK4_AUDIT_STEP_MAX=1000000
export GFS_RK4_AUDIT_STAGE=k1_rhs GFS_RK4_AUDIT_IDS=55283
set +e
/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 11580 > log 2>&1
rc=$?
set -e
printf '%s\\n' "$rc" > exit_status.txt
'''
        (r/'run.sh').write_text(script)
        subprocess.run(['bash','-n',str(r/'run.sh')],check=True)
        cases.append({'case':name,'signed_plane_offset':delta,'actual_plane_offset':nx*(shifted[0]-3)+ny*(shifted[1]-3)-apothem,
                      'snapshot_sha256':sha(snapshot),'original_xy':vals,'probe_xy':shifted})
    config={'initial':str(initial),'initial_sha256':sha(initial),'time':t,'site':55283,'side':543,
            'unchanged_binary_sha256':sha(parent/'eunha2.exe'),'cases':cases,
            'scope':'First native K1 RHS only at modified diagnostic coordinates; no production state, reflection, or accepted outside-domain evolution',
            'output_policy':'Three 32 MB probe inputs plus at most one 32 MB accepted inside endpoint and compact logs/stage TSV; no periodic dumps'}
    (root/'stage.json').write_text(json.dumps(config,indent=2)+'\n')
    runner=f'''#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
if [[ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]]; then echo 'GPU occupied' >&2;exit 31;fi
'''+''.join(f'bash {root/name}/run.sh\n' for name in ['inside','on','outside'])
    (root/'run_ada.sh').write_text(runner)
    print(root)


if __name__=='__main__':main()
