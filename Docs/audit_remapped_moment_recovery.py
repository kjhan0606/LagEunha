"""Read existing actual-t2 remap states; no new remap or physical evolution."""
import hashlib,json
from pathlib import Path
import numpy as np
from run_canonical_av_branches_20261018 import read
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_checkpoint_weight_remap_20261018'
ic,_=read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
out={'scope':__doc__,'cases':{}}
for fraction in [.25,1.]:
    path=root/f'state_w{fraction:g}_q8.npz';data=np.load(path);q=data['state'];ids=data['ids']
    m=q[:,0];c=q[:,1:3]/m[:,None];v=q[:,6:8]/m[:,None]
    Q=np.array([[q[:,3],q[:,4]],[q[:,4],q[:,5]]]).transpose(2,0,1)
    Pi=np.array([[q[:,12],q[:,13]],[q[:,13],q[:,14]]]).transpose(2,0,1)
    D=q[:,8:12].reshape(-1,2,2)
    S=Q-m[:,None,None]*c[:,:,None]*c[:,None,:]
    C=D-m[:,None,None]*c[:,:,None]*v[:,None,:]
    T=Pi-m[:,None,None]*v[:,:,None]*v[:,None,:]
    assert np.linalg.eigvalsh(S).min()>0
    G=np.linalg.solve(S,C).transpose(0,2,1)
    resolved=G@S@G.transpose(0,2,1);R=T-resolved;R=.5*(R+R.transpose(0,2,1))
    eig=np.linalg.eigvalsh(R);U=q[:,15]-.5*(q[:,12]+q[:,14])
    assert U.min()>0 and eig.min()>-1e-12
    affine=.5*np.trace(resolved,axis1=1,axis2=2);residual=.5*np.trace(R,axis1=1,axis2=2)
    internal=.5*np.trace(T,axis1=1,axis2=2)
    oldv=np.column_stack((ic['vx'][ids],ic['vy'][ids]));x=np.column_stack((ic['x'][ids]-3,ic['y'][ids]-3))
    recovered=v+np.einsum('nij,nj->ni',G,x-c);r=np.linalg.norm(x,axis=1)
    row={'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'minimum_thermal_energy':float(U.min()),
         'minimum_residual_covariance_eigenvalue':float(eig.min()),
         'kinetic_split_max_abs_residual':float(abs(affine+residual-internal).max()),'regions':{}}
    for name,mask in [('inner',r<.5),('edge',(r>=.5)&(r<.6)),('outer_patch',(r>=.6)&(r<.7))]:
        row['regions'][name]={'cells':int(mask.sum()),'affine_K_over_thermal':float(affine[mask].sum()/U[mask].sum()),
            'residual_K_over_thermal':float(residual[mask].sum()/U[mask].sum()),
            'affine_share_of_internal_K':float(affine[mask].sum()/internal[mask].sum()),
            'recovered_generator_velocity_change_mass_RMS':float(np.sqrt(np.average(np.sum((recovered-oldv)[mask]**2,axis=1),weights=m[mask])))}
    out['cases'][str(fraction)]=row
out['limits']=['Conditional lifted/remapped patch only; not a direct unmodified point-state representation.',
               'Recovered gradient is supported by stored moments, not an independently imposed neighbor slope.',
               'Positive residual covariance is necessary, not a complete spatial realizability or Euler closure proof.',
               'Generator velocity change after remap is not established as an accuracy improvement.']
with (root/'recovery_audit.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
