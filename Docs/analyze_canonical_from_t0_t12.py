"""Existing near-t12 snapshots: from-t0 timestep cap versus adaptive HG run.
Not a temporal convergence order estimate; binaries also have audited changes.
"""
import json
import sys
sys.path.insert(0,"/gpfs/kjhan/LagForce")
import struct
from pathlib import Path
import numpy as np
from run_canonical_av_branches_20261018 import read,sha
from audit_canonical_pair_history import gravity_coefficient

base=Path('/gpfs/kjhan/LagForce')
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
x0=np.column_stack((ic['x']-3,ic['y']-3));r0=np.linalg.norm(x0,axis=1)
theta0=np.arctan2(x0[:,1],x0[:,0]);om0=np.sqrt(gravity_coefficient(r0))
masks={'inner_background':(r0<.5)&(ic['den']<.5),
       'inner_disk_edge':(r0<.6)&(ic['den']>.5),
       'disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
out={'scope':__doc__,'cases':{}}
for name,root in [('adaptive_HG',base/'Codex_HG_constrained_20261017/omp8_t20'),
                  ('from_t0_cap0025',base/'Codex_t0_maxdt0025_t20_20261017/run')]:
    candidates=[]
    for path in root.glob('khout.*.dat'):
        with path.open('rb') as fp:
            header=fp.read(16)
            if len(header)!=16:continue
            magic,n,size,stress=struct.unpack('<4i',header)
            if path.stat().st_size!=16+n*size+52:continue
            fp.seek(16+n*size);t=struct.unpack('<d',fp.read(8))[0]
            if abs(t-12)<.02:candidates.append((abs(t-12),path))
    assert candidates,name
    _,path=min(candidates);a,t=read(path)
    assert np.array_equal(a['mass'],ic['mass']) and np.all(a['w2']==0)
    assert all(np.isfinite(a[k]).all() for k in a.dtype.names)
    x=np.column_stack((a['x']-3,a['y']-3));r=np.linalg.norm(x,axis=1)
    v=np.column_stack((a['vx'],a['vy']))
    vc=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack((-x[:,1],x[:,0]))
    phase=np.angle(np.exp(1j*(np.arctan2(x[:,1],x[:,0])-theta0-om0*t)))
    row={'snapshot':str(path),'sha256':sha(path),'time':t,'mass_unchanged':True,'regions':{}}
    for label,mask in masks.items():
        rms=lambda z:float(np.sqrt(np.average(z[mask]**2,weights=ic['mass'][mask])))
        row['regions'][label]={'pressure_error_mass_RMS':rms(a['pressure']/ic['pressure']-1),
           'phase_error_mass_RMS_deg':rms(phase)*180/np.pi,
           'radial_displacement_mass_RMS':rms(r-r0),
           'circular_velocity_error_mass_RMS':rms(np.linalg.norm(v-vc,axis=1))}
    out['cases'][name]=row
target=base/'Codex_t0_maxdt0025_t20_20261017/t12_comparison.json'
with target.open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
