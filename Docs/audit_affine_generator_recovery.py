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
qflux=accumulate(massface[:,None]*(x[i]+c)+(A**3/12*delta1)[:,None]*tangent)
recovered=(pdot-qflux@G.T)/a[:,6,None]
assert np.max(abs(recovered[active]))<1e-12
# Independent two-point Gauss integration of the quadratic momentum integrand.
pgauss=np.zeros_like(p_mid)
for sign in [-1,1]:
    z=sign*A/(2*np.sqrt(3))
    pgauss += (.5*A*(delta0+delta1*z))[:,None]*(vf+z[:,None]*gt)
quadrature_error=float(np.max(abs((pgauss-p_mid-p_correction)[active[i]])))
assert quadrature_error<1e-12
out={'gauss_vs_analytic_momentum_flux_max_error':quadrature_error,'scope':__doc__,'W':0,'relative_weight_rates_nonzero':True,
     'reduction_limit':'Equal weights at this instant do not invoke persistently-equal-weight reduction.',
     'regions':{}}
for name,mask in [('inner',r<.5),('edge',(r>=.5)&(r<.6)),('core',(r>=.6)&(r<1.9))]:
    rms=lambda z:float(np.sqrt(np.average(z[mask]**2,weights=a[mask,6])))
    out['regions'][name]={'fractional_net_mass_rate_RMS':rms(mdot/a[:,6]),
       'reconstructed_generator_acceleration_RMS':rms(np.linalg.norm(recovered,axis=1)),
       'momentum_transport_per_mass_RMS':rms(np.linalg.norm(pdot,axis=1)/a[:,6]),
       'midpoint_momentum_transport_error_RMS':rms(np.linalg.norm(pdot-pmid,axis=1)/a[:,6]),
       'absolute_face_mass_exchange':float(abs(massface[mask[i]]).sum())}
out['interpretation']='For manufactured steady simple shear G^2=0, centroid-aware v(X)=p/m-G(q/m-X) removes this transport-induced apparent generator acceleration. Frozen affine identity only; no evolved closure or Kepler cure.'
with (root/'affine_generator_recovery.json').open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
print(json.dumps(out,indent=2))
