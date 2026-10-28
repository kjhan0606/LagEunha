"""Rigid-rotation projection of existing remapped moments; no evolution.

Uses mass-center inertia, not generator-centered inertia. These are conditional
lift/remap states, not newly measured physical spin in the production solver.
"""
import hashlib,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
base=Path('/gpfs/kjhan/LagForce')
root=base/'Codex_checkpoint_weight_remap_20261018'
ic,_=read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
J=np.array([[0.,-1.],[1.,0.]])
out={'scope':__doc__,'cases':{}}
for fraction in [.25,1.]:
    path=root/f'state_w{fraction:g}_q8.npz'
    with np.load(path) as data:q=data['state'];ids=data['ids']
    m=q[:,0];c=q[:,1:3]/m[:,None];v=q[:,6:8]/m[:,None]
    Q=np.stack((q[:,3],q[:,4],q[:,4],q[:,5]),axis=1).reshape(-1,2,2)
    Pi=np.stack((q[:,12],q[:,13],q[:,13],q[:,14]),axis=1).reshape(-1,2,2)
    S=Q-m[:,None,None]*c[:,:,None]*c[:,None,:]
    C=q[:,8:12].reshape(-1,2,2)-m[:,None,None]*c[:,:,None]*v[:,None,:]
    T=Pi-m[:,None,None]*v[:,:,None]*v[:,None,:]
    assert np.linalg.eigvalsh(S).min()>0
    G=np.linalg.solve(S,C).transpose(0,2,1)
    ell=C[:,0,1]-C[:,1,0];I=np.trace(S,axis1=1,axis2=2);omega=ell/I
    H=G-omega[:,None,None]*J
    rotation=.5*ell**2/I
    deformation=.5*np.trace(H@S@H.transpose(0,2,1),axis1=1,axis2=2)
    R=T-G@S@G.transpose(0,2,1)
    residual=.5*np.trace(R,axis1=1,axis2=2)
    internal=.5*np.trace(T,axis1=1,axis2=2)
    thermal=q[:,15]-.5*np.trace(Pi,axis1=1,axis2=2)
    remainder_cross=S@H.transpose(0,2,1)
    remaining_spin=remainder_cross[:,0,1]-remainder_cross[:,1,0]
    split_error=float(np.max(abs(rotation+deformation+residual-internal)))
    assert thermal.min()>0 and split_error<1e-12
    assert np.max(abs(remaining_spin))<1e-12
    position=np.column_stack((ic['x'][ids]-3,ic['y'][ids]-3));radius=np.linalg.norm(position,axis=1)
    row={'source':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'split_max_absolute_error':split_error,
         'remaining_affine_spin_max_absolute':float(np.max(abs(remaining_spin))),
         'minimum_residual_covariance_eigenvalue':float(np.linalg.eigvalsh(.5*(R+R.transpose(0,2,1))).min()),
         'regions':{}}
    for name,mask in [('inner',radius<.5),('edge',(radius>=.5)&(radius<.6))]:
        total=float(internal[mask].sum())
        row['regions'][name]={'cells':int(mask.sum()),'internal_kinetic_energy':total,
            'rotation':float(rotation[mask].sum()),'affine_nonrigid':float(deformation[mask].sum()),
            'residual':float(residual[mask].sum()),'rotation_share':float(rotation[mask].sum()/total),
            'rotation_over_thermal':float(rotation[mask].sum()/thermal[mask].sum()),
            'mass_RMS_difference_spin_vs_half_curl':float(np.sqrt(np.average(
                (omega[mask]-.5*(G[mask,1,0]-G[mask,0,1]))**2,weights=m[mask])))}
    out['cases'][str(fraction)]=row
out['limits']=['No torque evolution or relaxation inferred from a static projection.',
 'Spin omega=L/I differs from half-curl in anisotropic cells.',
 'Affine nonrigid component need not equal the symmetric-gradient component.',
 'Residual kinetic energy remains distinct from thermal energy.',
 'No production state, binary, force or running calculation changed.']
print(json.dumps(out,indent=2))
