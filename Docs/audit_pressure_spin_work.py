"""Pressure torque and inertia work on saved t0 polygons, manufactured fields.

Nonuniform affine density separates mass center from area center. Reuses
completed CPU/OpenMP/CUDA pressure-moment outputs, not a new fluid evolution.
"""
import json
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
root=base/'Codex_pressure_rate_backend_v2_20261018'
raw=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rhs.rank0000')
raw=raw[np.argsort(raw[:,0])];x=raw[:,1:3]-3;n=len(x)
assembly=np.load(root/'assembly.npz');active=assembly['active'];owner=assembly['owner']
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000')
faces=faces[active[faces[:,0].astype(int)]]
assert np.array_equal(owner,faces[:,0].astype(int))
a=faces[:,6:8];b=faces[:,8:10]
cross=a[:,0]*b[:,1]-a[:,1]*b[:,0]
orientation=np.sign(np.bincount(owner,weights=cross,minlength=n))
cross*=orientation[owner]
nodes,weights=np.polynomial.legendre.leggauss(4);nodes=(nodes+1)/2;weights=weights/2
m=np.zeros(n);q=np.zeros((n,2));Q=np.zeros((n,2,2));areaq=np.zeros((n,2))
def scatter(z):return np.bincount(owner,weights=z,minlength=n)
for u,wu in zip(nodes,weights):
 for v,wv in zip(nodes,weights):
    pos=x[owner]+u*a+(1-u)*v*b
    dv=cross*(1-u)*wu*wv
    rho=1+.08*pos[:,0]-.04*pos[:,1];assert np.all(rho>0)
    dm=dv*rho;m+=scatter(dm)
    for k in range(2):
        q[:,k]+=scatter(dm*pos[:,k]);areaq[:,k]+=scatter(dv*pos[:,k])
        for l in range(2):Q[:,k,l]+=scatter(dm*pos[:,k]*pos[:,l])
idx=np.flatnonzero(active);c=q[idx]/m[idx,None]
S=Q[idx]-m[idx,None,None]*c[:,:,None]*c[:,None,:]
assert np.linalg.eigvalsh(S).min()>0
G=assembly['G'];C=S@G.T;I=np.trace(S,axis1=1,axis2=2)
ell=C[:,0,1]-C[:,1,0];omega=ell/I
Idot=np.trace(G@S+S@G.T,axis1=1,axis2=2)
gp=np.array([.001,.002]);V=assembly['volume'][idx]
lever=areaq[idx]-V[:,None]*c
expected=-lever[:,0]*gp[1]+lever[:,1]*gp[0]
out={'scope':__doc__,'cells':len(idx),'density':'1+.08*x-.04*y, coordinates relative to (3,3)',
     'max_abs_expected_torque':float(abs(expected).max()),'backends':{}}
for name in ['cpu1','cpu4','cuda']:
    rates=np.fromfile(root/(name+'.bin'),dtype='<f8').reshape(-1,17)
    r=np.column_stack([scatter(rates[:,k]) for k in range(17)])[idx]
    # Pressure qdot=mdot=0: derivative of central angular momentum.
    tau=r[:,9]-r[:,10]-(c[:,0]*r[:,7]-c[:,1]*r[:,6])
    err=float(np.max(abs(tau-expected)/V));assert err<1e-10
    work=omega*tau-.5*omega**2*Idot
    fd=[]
    for h in [1e-3,5e-4,2.5e-4]:
        assert np.all(I-h*Idot>0)
        plus=.5*(ell+h*tau)**2/(I+h*Idot)
        minus=.5*(ell-h*tau)**2/(I-h*Idot)
        fd.append(float(np.max(abs((plus-minus)/(2*h)-work))))
    assert fd[-1]<fd[0]
    out['backends'][name]={'torque_max_volume_scaled_error':err,
       'inertia_rate_max_absolute':float(abs(Idot).max()),
       'torque_work_sum':float(np.sum(omega*tau)),
       'inertia_work_sum':float(np.sum(-.5*omega**2*Idot)),
       'rotation_energy_rate_sum':float(np.sum(work)),
       'energy_derivative_centered_difference_errors':fd}
out['limits']=['Manufactured smooth fields on actual polygons, not measured Kepler forces.',
 'Pressure moment backends reused; no new CUDA full evolution or MPI claim.',
 'Inertia rate is an instantaneous material-affine rate. Moving Voronoi/Laguerre cells require ALE flux terms.',
 'Finite-difference path checks the chain rule only, not Euler time evolution.',
 'Pressure torque is reversible exchange, not a specified dissipative relaxation law.']
print(json.dumps(out,indent=2))
