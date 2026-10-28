"""Frozen-geometry relaxed volume constraint; not a nonlinear simulation."""
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

root=Path('/gpfs/kjhan/LagForce/Codex_pair_hg_modes_t0_20261018')
a=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');a=a[np.argsort(a[:,0])]
f=np.loadtxt(root/'raw_faces.dat.common.rank0000');i,j=f[:,:2].astype(int).T
n=len(a);d=np.linalg.norm(f[:,2:4],axis=1);k=f[:,4]/(2*d)
L=coo_matrix((np.r_[k,-k],(np.r_[i,i],np.r_[i,j])),shape=(n,n)).tocsr()
r=np.linalg.norm(a[:,1:3]-3,axis=1);active=r<1.9
u=np.zeros(n);u[active]=spsolve(L[active][:,active],a[active,12])
res=float(np.linalg.norm((L@u-a[:,12])[active])/np.linalg.norm(a[active,12]))
assert res<1e-10,res
out={'scope':__doc__,'equation':'L Wdot = -Bv-L W/tau on r<1.9, outer W=0.',
     'solve_relative_residual':res,'cases':[],
     'limits':['Uniform incompressible density for mass interpretation, not actual density interface.',
       'Frozen L and Bv only: bounded solution here proves no nonlinear long-time bound.',
       'Full face flux still required; equal-weight original-solver reduction unproved.']}
touch=active[i]|active[j]
for tau in [.01,.1,1.]:
    W=-tau*u
    massrate=-L@W/tau # rho=1 manufactured reference
    err=float(np.max(abs((massrate-a[:,12])[active])))
    out['cases'].append({'tau':tau,'stationary_weight_difference_over_d2_max':float(np.max(abs((W[i]-W[j])[touch])/d[touch]**2)),
       'stationary_uniform_density_massrate_identity_max_error':err,
       'stationary_Wdot':0,
       'nonzero_massrate_L1_active':float(abs(massrate[active]).sum())})
with (root/'relaxed_weight_tangent.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
