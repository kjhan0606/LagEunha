"""Prepare retained t2 remap chord states and independent long-double ledger."""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_moment_spin_backends_20261018'
a,b=[np.load(base/f'Codex_checkpoint_weight_remap_20261018/state_w{f}_q8.npz')['state'] for f in ('0.25','1')]
q=(a+b)/2;d=b-a
z=q.astype(np.longdouble);dz=d.astype(np.longdouble)
m=z[:,0];md=dz[:,0];c=z[:,1:3]/m[:,None];cd=(dz[:,1:3]-c*md[:,None])/m[:,None]
p=z[:,6:8];pd=dz[:,6:8]
I=z[:,3]+z[:,5]-m*np.sum(c*c,axis=1)
Id=dz[:,3]+dz[:,5]-md*np.sum(c*c,axis=1)-2*m*np.sum(c*cd,axis=1)
cross=lambda a,b:a[:,0]*b[:,1]-a[:,1]*b[:,0]
ell=z[:,9]-z[:,10]-cross(c,p);elld=dz[:,9]-dz[:,10]-cross(cd,p)-cross(c,pd)
omega=ell/I
ref=np.column_stack([I,ell,omega,.5*ell*omega,Id,elld,omega*elld-.5*omega**2*Id,dz[:,15]-.5*(dz[:,12]+dz[:,14])])
inputs=np.column_stack([q,d]);guard=[]
for case in ('zero_mass','negative_inertia','nan'):
    row=inputs[0].copy()
    if case=='zero_mass':row[0]=0.;code=1
    elif case=='negative_inertia':row[3]=row[5]=-1.;code=2
    else:row[7]=np.nan;code=1
    inputs=np.vstack([inputs,row]);guard.append(dict(case=case,code=code))
with (root/'input.bin').open('xb') as f:inputs.astype('<f8').tofile(f)
with (root/'reference.npy').open('xb') as f:np.save(f,np.asarray(ref,dtype=float))
(root/'cases.json').write_text(json.dumps(dict(valid=len(q),guards=guard)))
