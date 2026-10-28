"""Stage first-RHS-only native probes and a flag-off parity replay."""
from pathlib import Path
import json
import math
import shutil
import struct
import subprocess
import sys
import numpy as np
from stage_kepler_wall_rhs_probe_20261029 import sha


def main():
    base=Path('/gpfs/kjhan/LagForce/Kepler_Q2_WallCrossing_Replay_20261029')
    old=base/'rhs_extension_probes'
    root=base/'native_band_probes';root.mkdir()
    candidate=Path('/gpfs/kjhan/LagForce/Kepler_SlipWall_EventsBuild_20261029/project/eunha2.exe')
    expected='3a9e98eaed65a903534a001f6c400c4ad815b7d873577ad52cd658ab1b580e45'
    if sha(candidate)!=expected:raise ValueError('Native probe candidate changed')
    source=base/'control_capture/khout.011580.dat'
    sys.path.insert(0,'/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
    from audit_snapshot import read_snapshot
    p,_=read_snapshot(source)
    ids=p['u4if'] & np.uint64((1<<56)-1)
    owner=int(np.flatnonzero(ids==55283)[0])
    xy=np.column_stack((p['x'],p['y']))
    distance=np.sum((xy-xy[owner])**2,axis=1);distance[owner]=np.inf
    other=int(np.argmin(distance))
    offsets=json.loads(Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/snapshot_offsets.json').read_text())['fields']
    theta=543.5*2*math.pi/1024.;normal=np.array([math.cos(theta),math.sin(theta)])
    apothem=2.75*math.cos(math.pi/1024.)
    tangent=np.array([-normal[1],normal[0]])
    tangential=float(np.dot(xy[owner]-3.,tangent))
    configs=[]
    choices=[('baseline',None),('inside',-1.e-7),('on',0.),('outside',1.e-7),
             ('maximum_band',1.e-6-1.e-12),('exceeds_band',1.1e-6),('empty_cell',1.e-7)]
    for name,offset in choices:
        r=root/name;r.mkdir();(r/'eunha2.exe').symlink_to(candidate)
        shutil.copy2(old/'inside/params.dat',r/'params.dat')
        snapshot=r/'khout.011580.dat';shutil.copy2(source,snapshot)
        updates={}
        if offset is not None:
            updates[owner]=3.+(apothem+offset)*normal+tangential*tangent
        if name=='empty_cell':
            updates[other]=3.+(apothem+5.e-8)*normal+tangential*tangent
        with snapshot.open('r+b') as f:
            for index,coord in updates.items():
                for axis,k in enumerate(['x','y']):
                    f.seek(16+index*576+offsets[k]['offset']);f.write(struct.pack('<d',coord[axis]))
        script=(old/'inside/run.sh').read_text().replace(str(old/'inside'),str(r))
        script=script.replace('set +e\n','export GFS_WALL_RHS_PROBE_ONLY=1 GFS_WALL_RHS_PROBE_BAND=1e-6\nset +e\n',1)
        (r/'run.sh').write_text(script)
        subprocess.run(['bash','-n',str(r/'run.sh')],check=True)
        configs.append({'case':name,'intended_offset':offset,'input_sha256':sha(snapshot),
                        'changed_site_ids':[int(ids[i]) for i in updates],
                        'expected_exit':1 if name in ['exceeds_band','empty_cell'] else 0})
    # Real accepted tiny step with both new settings absent, same input as old.
    r=root/'flag_off_parity';r.mkdir();(r/'eunha2.exe').symlink_to(candidate)
    for name in ['params.dat','khout.011580.dat']:shutil.copy2(old/'inside'/name,r/name)
    script=(old/'inside/run.sh').read_text().replace(str(old/'inside'),str(r))
    script=script.replace('set +e\n','unset GFS_WALL_RHS_PROBE_ONLY GFS_WALL_RHS_PROBE_BAND\nset +e\n',1)
    (r/'run.sh').write_text(script)
    config={'candidate':str(candidate),'candidate_sha256':expected,'cases':configs,
            'configured_band':1.e-6,'actual_band_test_margin':1.e-12,
            'source':str(source),'source_sha256':sha(source),
            'scope':'First RHS diagnostic snapshots only. No external state accepted or integrated. Separate flag-off inside replay must preserve old full native result.',
            'output_policy':'About 400 MB inputs/first-RHS snapshots plus 32 MB parity endpoint, compact TSV/logs; no periodic dumps'}
    (root/'stage.json').write_text(json.dumps(config,indent=2)+'\n')
    runner='''#!/bin/bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
if [[ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]]; then echo 'GPU occupied' >&2;exit 31;fi
'''+''.join(f'bash {root/name}/run.sh\n' for name in [x['case'] for x in configs]+['flag_off_parity'])
    (root/'run_ada.sh').write_text(runner)
    print(root)


if __name__=='__main__':main()
