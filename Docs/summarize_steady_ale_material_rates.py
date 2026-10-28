"""Derived moving-cell rates, with mass chain rule and initial pressure rate."""
import json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base))
from run_canonical_av_branches_20261018 import read
h=np.load(base/'Codex_continuous_initial_rhs_quadratic_split_20261018/rhs_r0.55_f8.npz')
s=np.load(base/'Codex_steady_ale_initial_balance_20261018/reference_q32.npz')
assert np.array_equal(h['ids'],s['ids'])
ids=s['ids'];m=s['mass'];a=h['state'][:-4].reshape(-1,20)
r=np.column_stack([m,a[:,3:]]);d=s['numerical']-s['reference'][:,:18]
v=r[:,7:9]/m[:,None];c=r[:,1:3]/m[:,None]
U=r[:,17]-.5*(r[:,13]+r[:,16]);u=U/m
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
gamma=float(np.float32(1.666667));area=ic['volume'][ids]
gravity_thermal=s['gravity'][:,17]-.5*(s['gravity'][:,13]+s['gravity'][:,16])
fields=dict(bulk_velocity_rate=(d[:,7:9]-v*d[:,0,None])/m[:,None],
 centroid_velocity=(d[:,1:3]-c*d[:,0,None])/m[:,None],
 specific_thermal_rate=(d[:,17]-.5*(d[:,13]+d[:,16])-u*d[:,0])/m,
 pressure_rate=h['pressure_rates']+(gamma-1)*gravity_thermal/area)
radius=np.linalg.norm(a[:,:2],axis=1);active=s['active'];regions={};largest={}
for name,mask in [('active',active),('inner_r_lt_040',active&(radius<.4)),
 ('r040_048',active&(radius>=.4)&(radius<.48)),('contact_r048_055',active&(radius>=.48))]:
    regions[name]={'cells':int(mask.sum())}
    for key,val in fields.items():
        q=val[mask].reshape(int(mask.sum()),-1)
        regions[name][key]=dict(max_abs=float(abs(q).max()),mass_weighted_RMS=float(np.sqrt(np.average((q*q).sum(axis=1),weights=m[mask]))))
for key,val in fields.items():
    magnitude=np.linalg.norm(val.reshape(len(ids),-1),axis=1)
    top=np.flatnonzero(active)[np.argsort(magnitude[active])[-10:][::-1]]
    largest[key]=[dict(id=int(ids[i]),radius=float(radius[i]),norm=float(magnitude[i]),mass=float(m[i])) for i in top]
print(json.dumps(dict(regions=regions,largest=largest,
 scope='Same exact-projection initial state. Model minus exact moving-boundary derivatives include nonzero exact dm. Pressure reference is spatially constant Pdot=0; model pressure rate includes volume change and gravity work. Model density/velocity gravity source; numerical Psi evolution absent. No material-particle acceleration interpretation or long-term validation.'),indent=2))
