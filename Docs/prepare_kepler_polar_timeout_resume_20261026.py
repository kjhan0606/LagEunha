"""Select the latest valid native polar checkpoint after a scheduler timeout.

This runs only inside an afternotok Slurm fallback whose parent state has
already been checked. It never accepts a failed RK4 stage diagnostic as a
restart state and never changes the parent run.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import struct
import sys

import numpy as np

ROOT = Path('/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020')
sys.path.insert(0, str(ROOT))
from audit_snapshot import read_snapshot


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def load_if_valid(path, mass_reference):
    match=re.fullmatch(r'khout\.(\d+)\.dat',path.name)
    if not match:
        return None
    try:
        with path.open('rb') as stream:
            magic,n,stride,_=struct.unpack('<4i',stream.read(16))
        if magic!=-1 or n!=384*384 or stride<=0 or path.stat().st_size!=16+n*stride+52:
            return None
        p,meta=read_snapshot(path)
        index=(p['u4if'].astype(np.uint64)&np.uint64((1<<56)-1)).astype(np.int64)
        if not np.array_equal(np.sort(index),np.arange(384*384)):
            return None
        p=p[np.argsort(index)]
        if not (np.array_equal(p['mass'],mass_reference) and np.all(p['w2']==0)
                and np.all(p['volume']>0) and np.all(p['ie']>0)
                and all(np.isfinite(p[k]).all() for k in ('x','y','vx','vy','mass','volume','ie','den','pressure'))):
            return None
        t=float(meta['t'])
        if not (.001-1e-12<=t<=.825+1e-12):
            return None
        return (t,int(match.group(1)),path)
    except (OSError,ValueError,struct.error):
        return None


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('preflight',type=Path)
    parser.add_argument('parent',type=Path)
    parser.add_argument('destination',type=Path)
    parser.add_argument('--parent-state',required=True)
    args=parser.parse_args()
    if args.parent_state not in ('TIMEOUT','PREEMPTED','NODE_FAIL','BOOT_FAIL'):
        raise ValueError(f'Unsafe parent state: {args.parent_state}')
    receipt=json.loads((args.preflight/'audit_t001.json').read_text())
    if not all(receipt['checks'].values()):
        raise ValueError('Preflight audit is not accepted')
    p0,meta0=read_snapshot(args.preflight/'khout.000000.dat')
    if meta0['t']!=0.:
        raise ValueError('Initial reference is not exact t=0')
    ids0=(p0['u4if'].astype(np.uint64)&np.uint64((1<<56)-1)).astype(np.int64)
    if not np.array_equal(np.sort(ids0),np.arange(384*384)):
        raise ValueError('Initial reference IDs changed')
    reference_mass=p0['mass'][np.argsort(ids0)]
    candidates=[]
    for path in args.parent.glob('khout.*.dat'):
        valid=load_if_valid(path,reference_mass)
        if valid is not None:
            candidates.append(valid)
    if not candidates:
        raise ValueError('No complete positive fixed-mass W=0 checkpoint')
    t,step,source=max(candidates)
    if args.destination.exists() and any(args.destination.glob('khout.*.dat')):
        raise ValueError('Destination already contains an evolution snapshot')
    args.destination.mkdir(parents=True,exist_ok=True)
    for original in (source,args.parent/'eunha2.exe',args.parent/'params.dat'):
        target=args.destination/original.name
        if target.exists():
            raise ValueError(f'Destination file exists: {target}')
        shutil.copy2(original,target)
    mode='audit_only' if abs(t-.825)<1e-12 else 'resume'
    input_hash={p.name:digest(args.destination/p.name) for p in
                (source,args.parent/'eunha2.exe',args.parent/'params.dat')}
    result=dict(mode=mode,parent_state=args.parent_state,source=str(source),
                source_sha256=digest(source),time=t,step=step,target_time=.825,
                input_sha256=input_hash)
    with (args.destination/'resume_plan.json').open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
