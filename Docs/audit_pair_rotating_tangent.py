"""Restricted pair tangent with central gravity and a frozen rotating frame."""
import sys,json
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from audit_canonical_pair_history import gravity_coefficient
base=Path('/gpfs/kjhan/LagForce')
data=json.loads((base/'Codex_pair_force_hessian_controls_20261018/result.json').read_text())
s=np.loadtxt(base/'Codex_qbp_t10_decomposition_20261018/raw_faces.dat.rhs.rank0000');s=s[np.argsort(s[:,0])]
i,j=data['pair'];x=s[[i,j],1:3]-3;m=s[[i,j],7];mu=m.prod()/m.sum();coef=np.array([-m[1],m[0]])/m.sum()
normal=(x[1]-x[0])/np.linalg.norm(x[1]-x[0]);B=np.array([normal,[-normal[1],normal[0]]])
center=(m[:,None]*x).sum(axis=0)/m.sum();omega=float(np.sqrt(gravity_coefficient(np.array([np.linalg.norm(center)])))[0])
def grav_stiff(eps):
 r=np.linalg.norm(x,axis=1);k=gravity_coefficient(r)
 rkprime=(gravity_coefficient(r*np.exp(eps))-gravity_coefficient(r*np.exp(-eps)))/(2*eps)
 hess=k[:,None,None]*np.eye(2)+rkprime[:,None,None]*x[:,:,None]*x[:,None,:]/r[:,None,None]**2
 return B@np.sum((m*coef**2)[:,None,None]*hess,axis=0)@B.T
G=grav_stiff(1e-5);G2=grav_stiff(5e-6);J=np.array([[0.,-1.],[1.,0.]])
out={'pair':[i,j],'reduced_mass':float(mu),'frame_omega':omega,
 'gravity_stiffness':G.tolist(),'gravity_derivative_check_max_abs':float(np.max(abs(G-G2))),
 'scope':'Instantaneous frozen coefficients; pair center-of-mass and all other generators constrained. Not a full variational trajectory, equilibrium or global eigenmode. Excludes AV/HG Jacobian and entropy dynamics.', 'models':{}}
for name,row in data['samples'][-1]['models'].items():
 H=np.array(row['matrix']);H=.5*(H+H.T);K=H+G-mu*omega**2*np.eye(2)
 A=np.block([[np.zeros((2,2)),np.eye(2)],[-K/mu,-2*omega*J]])
 eig=np.linalg.eigvals(A)
 out['models'][name]={'thermal_stiffness_per_mass_eigenvalues':np.linalg.eigvalsh(H/mu).tolist(),
 'effective_stiffness_per_mass_eigenvalues':np.linalg.eigvalsh(K/mu).tolist(),
 'tangent_eigenvalues_real_imag':[[float(v.real),float(v.imag)] for v in eig],
 'max_real_part':float(np.max(eig.real)),
 'eigen_residual':float(np.max(abs(np.linalg.det(A-eig[0]*np.eye(4)))))}
p=base/'Codex_pair_rotating_tangent_20261018';p.mkdir(exist_ok=False)
(p/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
