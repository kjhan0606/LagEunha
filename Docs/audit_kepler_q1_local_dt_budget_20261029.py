"""Frozen-state update-count budget for the completed early Q1 checkpoint."""
import hashlib
import json
import socket
import struct
from pathlib import Path
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global,ID_MASK

if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run on compute host')
docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
snapshot=Path('/gpfs/kjhan/LagForce/Kepler_Q1_OrbitalCap006_Ada_20261029/to0p4125/khout.005400.dat')
initial=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/ic_v1/khout.000001.dat')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
if sha(snapshot)!='df1eca2b0c341d0c4e61842c6e7ac8ec98ed88e24aea48ef63c1b19f6c57c924':raise ValueError('changed snapshot')
a,ma=ordered_global(initial);p,meta=ordered_global(snapshot,a['u4if']&ID_MASK)
if ma['t']!=0 or meta['t']!=.4125 or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('contract mismatch')
with snapshot.open('rb') as f:
    magic,n,stride,stress=struct.unpack('<4i',f.read(16))
if (magic,n,stride,stress)!=(-1,90352,576,1):raise ValueError('wrong layout')
dtype=np.dtype({'names':['u4if','x','y','dt'],'formats':['<u8','<f8','<f8','<f4'],'offsets':[8,24,32,128],'itemsize':576})
raw=np.memmap(snapshot,dtype=dtype,mode='r',offset=16,shape=(n,))
ids=raw['u4if']&ID_MASK;order=np.argsort(ids)
if not np.array_equal(ids[order],a['u4if']&ID_MASK):raise ValueError('changed IDs')
raw=raw[order];hydro=raw['dt'].astype(float)
if not np.all(np.isfinite(hydro)) or np.any(hydro<=0):raise ValueError('invalid saved CFL')
radius=np.hypot(raw['x']-3,raw['y']-3);r0=np.hypot(a['x']-3,a['y']-3)
orbital=.06*(radius*radius+.01175**2)**.75
local=np.minimum(hydro,orbital);base=float(local.min())
level=np.floor(np.log2(local/base)).astype(int);block=base*np.exp2(level)
span=.4125;global_steps=int(np.ceil(span/base));block_updates=np.ceil(span/block).astype(np.int64)
regions=[]
for lo,hi in ((0,.1),(.1,.5),(.5,.75),(.75,1),(1,2),(2,2.75)):
    mask=(r0>=lo)&(r0<hi)
    regions.append({'r0_min':lo,'r0_max':hi,'cells':int(mask.sum()),
                    'combined_dt_min':float(np.min(local[mask])),
                    'combined_dt_median':float(np.median(local[mask])),
                    'hydro_dt_min':float(np.min(hydro[mask])),
                    'block_cell_updates':int(np.sum(block_updates[mask],dtype=np.int64))})
all_updates=int(n*global_steps);block_total=int(np.sum(block_updates,dtype=np.int64))
result={'time':meta['t'],'snapshot_sha256':sha(snapshot),'initial_sha256':sha(initial),'cells':n,'eta':.06,
        'hydro_dt_min':float(hydro.min()),'orbital_dt_min':float(orbital.min()),'combined_dt_min':base,
        'span':span,'global_steps':global_steps,'global_cell_updates':all_updates,'block_cell_updates':block_total,
        'global_to_block_cell_update_ratio':all_updates/block_total,
        'levels':[{'level':int(k),'cells':int(np.sum(level==k))} for k in np.unique(level)],
        'regions':regions,'script_sha256':sha(Path(__file__)),
        'scope':'Frozen-state arithmetic using saved final-stage hydro dt and local orbital cap from one legacy-thermal checkpoint; assumes local independent power-of-two cell updates, does not model conservative face synchronization, geometry rebuild, evolving timesteps, gravity error, wall events, MPI, overhead or speed. Not an implemented individual time integrator or estimate of end-to-end runtime.'}
out=docs/'KEPLER_Q1_GRAVITY_CAP_2026-10-29/LOCAL_DT_BUDGET_04125.json'
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
