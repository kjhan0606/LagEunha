"""Exact line-integrated affine transport on saved W=0 faces.
Uniform rho=1, simple shear v=(y,0); diagnostic only, no state conversion.
"""
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
f=np.loadtxt(root/'raw_faces.dat.common.rank0000');i,j=f[:,:2].astype(int).T
x=a[:,1:3]-3;r=np.linalg.norm(x,axis=1);n=len(a);active=r<1.9
dx=f[:,2:4];A=f[:,4];d=np.linalg.norm(dx,axis=1);e=dx/d[:,None]
tangent=np.column_stack((-e[:,1],e[:,0]));c=f[:,5:7];off=c-.5*dx
G=np.array([[0.,1.],[0.,0.]])
vi=x[i]@G.T;vj=vi+dx@G.T;vf=(x[i]+c)@G.T;gt=tangent@G.T
s0=.5*np.sum((vi+vj)*e,axis=1)+np.sum((vi-vj)*off,axis=1)/d
s1=np.sum((vi-vj)*tangent,axis=1)/d
B=np.bincount(i,weights=A*s0,minlength=n)
k=A/(2*d);L=coo_matrix((np.r_[k,-k],(np.r_[i,i],np.r_[i,j])),shape=(n,n)).tocsr()
wdot=np.zeros(n);wdot[active]=spsolve(L[active][:,active],-B[active])
delta0=s0+(wdot[i]-wdot[j])/(2*d)-np.sum(vf*e,axis=1)
delta1=s1-np.sum(gt*e,axis=1)
massface=A*delta0
p_mid=massface[:,None]*vf
p_correction=(A**3/12*delta1)[:,None]*gt
def accumulate(z):
    return np.column_stack([np.bincount(i,weights=z[:,k],minlength=n) for k in range(z.shape[1])])
mdot=np.bincount(i,weights=massface,minlength=n)
pdot=accumulate(p_mid+p_correction);pmid=accumulate(p_mid)
assert np.max(abs(mdot[active]/a[active,6]))<1e-7

target=Path('/gpfs/kjhan/LagForce/Codex_affine_flux_backend_20261018')
target.mkdir(exist_ok=False)
inputs=np.column_stack((A,delta0,delta1,x[i]+c,tangent,vf,gt))[active[i]]
rng=np.random.default_rng(20261018)
extra=rng.normal(size=(128,11));extra[:,0]=abs(extra[:,0]);extra[0,0]=0
inputs=np.vstack((inputs,extra))
inputs.astype('<f8').tofile(target/'inputs.bin')
expected=np.zeros((len(inputs),6))
for sign in [-1,1]:
    z=sign*inputs[:,0]/(2*np.sqrt(3))
    flux=.5*inputs[:,0]*(inputs[:,1]+inputs[:,2]*z)
    pos=inputs[:,3:5]+inputs[:,5:7]*z[:,None]
    vel=inputs[:,7:9]+inputs[:,9:11]*z[:,None]
    expected += flux[:,None]*np.column_stack((np.ones(len(inputs)),pos,vel,.5*np.sum(vel*vel,axis=1)))
expected.astype('<f8').tofile(target/'expected.bin')
(target/'manifest.json').write_text(json.dumps({'actual_geometry_faces':int(active[i].sum()),'manufactured_general_affine_cases':128,'scope':'Straight-face algebra only; no production integration or MPI.','reference':'Two-point Gauss, exact to degree3.'},indent=2)+'\n')
print(str(target))
