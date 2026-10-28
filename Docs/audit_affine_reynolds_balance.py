"""Moving-cell Reynolds identity for stationary simple shear on closed t0 cells."""
import json
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
x=a[:,1:3]-3;n=len(a);active=np.linalg.norm(x,axis=1)<1.9
f=np.loadtxt(root/'raw_faces.dat.rank0000');f=f[active[f[:,0].astype(int)]]
i=f[:,0].astype(int);dx=f[:,4:6];va=f[:,6:8];vb=f[:,8:10]
d=np.linalg.norm(dx,axis=1);normal=dx/d[:,None]
line=vb-va;length=np.linalg.norm(line,axis=1);tangent=line/length[:,None]
cf=.5*(va+vb);cross=va[:,0]*vb[:,1]-va[:,1]*vb[:,0]
area=.5*abs(cross)
def fields(pos):
    xx,yy=pos.T;zero=np.zeros(len(pos));one=np.ones(len(pos))
    # v=(y,0), rho=1, specific thermal energy=.001.
    return np.column_stack((one,xx,yy,xx*xx,xx*yy,yy*yy,yy,zero,xx*yy,zero,yy*yy,zero,yy*yy,zero,zero,.001+.5*yy*yy,zero))
def accum(z):return np.column_stack([np.bincount(i,weights=z[:,k],minlength=n) for k in range(17)])
mom=np.zeros((n,17))
for u,v in [(1/6,1/6),(2/3,1/6),(1/6,2/3)]:mom+=accum((area/3)[:,None]*fields(x[i]+u*va+v*vb))
assert np.max(abs(mom[active,0]/a[active,6]-1))<1e-10
source=np.zeros_like(mom)
source[:,1:3]=mom[:,6:8];source[:,3]=2*mom[:,8];source[:,4]=mom[:,9]+mom[:,10];source[:,5]=2*mom[:,11]
source[:,8]=mom[:,12];source[:,9]=mom[:,13];source[:,10]=mom[:,13];source[:,11]=mom[:,14]
vi=np.column_stack((x[i,1],np.zeros(len(i))));vj=vi+np.column_stack((dx[:,1],np.zeros(len(i))))
mesh0=.5*np.sum((vi+vj)*normal,axis=1)+np.sum((vi-vj)*(cf-.5*dx),axis=1)/d
mesh1=np.sum((vi-vj)*tangent,axis=1)/d
transport=np.zeros_like(mom);boundary=np.zeros_like(mom)
nodes,weights=np.polynomial.legendre.leggauss(3)
for z0,w in zip(nodes,weights):
    z=.5*length*z0;pos=x[i]+cf+z[:,None]*tangent;h=fields(pos)
    fluidnormal=pos[:,1]*normal[:,0];mesh=mesh0+mesh1*z
    boundary+=accum((.5*length*w*mesh)[:,None]*h)
    transport+=accum((.5*length*w*(mesh-fluidnormal))[:,None]*h)
error=source+transport-boundary
scaled=abs(error[active])/a[active,6,None]
assert scaled.max()<1e-10,scaled.max()
out={'scope':__doc__,'cells':int(active.sum()),'max_volume_normalized_residual_per_component':scaled.max(axis=0).tolist(),
     'missing_source_error_max_per_component':(abs((transport-boundary)[active])/a[active,6,None]).max(axis=0).tolist(),
     'limits':'Closed per-cell faces, uniform-density manufactured shear, no gravity or evolved dynamics; shared-face global conservation not checked.'}
with (root/'affine_reynolds_balance.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
target=root/'reynolds_source_inputs.bin';mom[active].astype('<f8').tofile(target)
source[active].astype('<f8').tofile(root/'reynolds_source_expected.bin')
print(json.dumps(out,indent=2))
