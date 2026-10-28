"""Instantaneous auxiliary work split on an unchanged saved physical state.

v = v_circular(x) + u is diagnostic only. B and z use the actual evolved
state; no force is altered and no reference velocity enters the simulation.
"""
import argparse,hashlib,json,struct,sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
from audit_canonical_pair_history import gravity_coefficient
from smooth_owner_memory import build_row
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('snapshot',type=Path);args=parser.parse_args();path=args.snapshot
a,t=read(path);ic,_=read(Path('/gpfs/kjhan/LagForce/Codex_canonical_initial_20261017/run/khout.000000.dat'))
data=Path(str(path)+'.rmem').read_bytes();header=struct.unpack('<6Q',data[:48])
assert header[:3]==(0x524d454d53544154,1,0x0102030405060708)
assert header[3]==len(a) and len(data)==56+24*len(a)
assert struct.unpack('<d',data[48:56])[0]==t
h=14695981039346656037
for byte in data[56:]:h=((h^byte)*1099511628211)&0xffffffffffffffff
assert h==header[5]
rec=np.frombuffer(data,dtype=[('id','<u8'),('x','<f8'),('y','<f8')],offset=56)
maskid=np.uint64(72057594037927935)
assert np.array_equal(rec['id'],a['id']&maskid) and np.array_equal(a['id']&maskid,ic['id']&maskid)
x=np.column_stack((a['x'],a['y']));xc=x-3;r=np.linalg.norm(xc,axis=1)
v=np.column_stack((a['vx'],a['vy']));z=np.column_stack((rec['x'],rec['y']))
assert np.isfinite(z).all()
vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-xc[:,1],xc[:,0]));u=v-vc
frozen=r>=2.5;tree=cKDTree(x%6,boxsize=6);gamma=float(np.float32(1.666667))
dz=np.zeros_like(v);dc=np.zeros_like(v);du=np.zeros_like(v);force=np.zeros_like(v);moment=0.
for i in np.flatnonzero(~frozen):
    ids,row,meta=build_row(x,i,tree);moment=max(moment,meta['moment_max'])
    b=np.sqrt(a['mass'][i])*np.sqrt(gamma*(gamma-1)*a['ie'][i]/a['mass'][i])/np.sqrt(a['volume'][i])*row
    dz[i]=b@v[ids];dc[i]=b@vc[ids];du[i]=b@u[ids]
    np.add.at(force,ids,-b[:,None]*z[i])
assert moment<1e-12
power=np.sum(z*dz,axis=1);pc=np.sum(z*dc,axis=1);pu=np.sum(z*du,axis=1)
closure=float(np.sum(force*v)+power.sum())
scale=float(np.sum(abs(force*v))+np.sum(abs(z*dz)))
assert abs(closure)<1e-12*max(scale,1e-300)
split=float(np.linalg.norm(dz-dc-du)/max(np.linalg.norm(dz),1e-300));assert split<1e-10
r0=np.linalg.norm(np.column_stack((ic['x']-3,ic['y']-3)),axis=1)
masks={'all_active':~frozen,'inner_background':(r0<.5)&(ic['den']<.5),
    'inner_shear_band':(r0>=.3)&(r0<.5),'inner_disk_edge':(r0<.6)&(ic['den']>.5),
    'disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
regions={}
for name,m in masks.items():
    free=m&~frozen
    regions[name]={'cells':int(m.sum()),'owner_auxiliary_power':float(power[m].sum()),
        'owner_power_from_circular_velocity':float(pc[m].sum()),
        'owner_power_from_departure_velocity':float(pu[m].sum()),
        'Bv_L2':float(np.linalg.norm(dz[m])),'Bvcircular_L2':float(np.linalg.norm(dc[m])),
        'Bdeparture_L2':float(np.linalg.norm(du[m])),
        'force_work_on_actual_velocity_free':float(np.sum(force[free]*v[free])),
        'force_work_on_circular_velocity_free':float(np.sum(force[free]*vc[free])),
        'force_work_on_departure_velocity_free':float(np.sum(force[free]*u[free]))}
out={'scope':__doc__,'snapshot':str(path),'snapshot_sha256':sha(path),
    'sidecar_sha256':hashlib.sha256(data).hexdigest(),'time':t,'regions':regions,
    'max_moment':moment,'velocity_split_relative_error':split,
    'global_power_closure_relative':abs(closure)/max(scale,1e-300),
    'frozen_boundary_external_power':-float(np.sum(force[frozen]*v[frozen])),
    'limitations':'Instantaneous auxiliary contribution only, not full departure-energy derivative or accumulated energy attribution. Owner regions and recipient-force regions differ because coupling crosses region boundaries. Circular-reference work is not automatically erroneous dissipation. SVD reconstruction is independently parity-checked elsewhere; no live state changed.'}
print(json.dumps(out,indent=2))
