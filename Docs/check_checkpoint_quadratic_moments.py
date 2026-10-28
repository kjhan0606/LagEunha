"""Projected curvature of actual t2 remap states; reintegrate stored moments."""
import hashlib,importlib.util,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'));sys.path.insert(0,str(base))
import moment_ale_evolve as model
import quadratic_velocity_reconstruction as quad
from run_canonical_av_branches_20261018 import read
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
geo=load('qmom_geometry',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
clip=load('qmom_clipper',base/'Codex_weighted_compression_geometry_20261017/audit.py')
ic,t=read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
x=np.column_stack([ic['x']-3,ic['y']-3]);w=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
os.environ.update(LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
alpha=np.sqrt(2)*np.array([[1.,0.],[-1.,0.],[0.,1.],[0.,-1.]])
rows=[]
for fraction in (.25,1.):
    path=base/f'Codex_checkpoint_covariance_factors_20261018/factors_w{fraction:g}.npz'
    data=np.load(path);ids=data['ids'];q=data['original_state'];mass=q[:,0];L=data['factors']
    polys=geo.polygons(x,fraction*w,ids,clip.clip);states=[]
    for poly,a in zip(polys,q):
        raw=np.r_[a[:3],a[3],a[4],a[4],a[5],a[6:12],a[12],a[13],a[13],a[14],a[15]]
        states.append(model.decode(poly,raw))
    affine=np.array([model.affine_velocity(s,point) for s,point in zip(states,x[ids])])
    diag=quad.prepare(states,polys,[],mass,L)
    integrated=[];change=[];residual_K=[]
    for i,s in enumerate(states):
        pts=s['density'].points;dm=s['density'].weighted_mass
        samples,_=quad.velocity_samples(s,L[i],mass[i],pts,alpha)
        mean=samples.mean(1);Pi=np.einsum('n,nki,nkj->ij',dm,samples,samples)/4
        p=dm@mean;D=pts.T@(dm[:,None]*mean)
        integrated.append(np.r_[p,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1]])
        change.append(model.velocity(s,x[ids[i]])-affine[i])
        meta=s.get('velocity_quadratic');M=meta['residual_map'] if meta else np.eye(2)
        residual_K.append(.5*np.sum((L[i]@M)**2))
    integrated=np.array(integrated);diff=integrated-q[:,6:15];change=np.array(change)
    assert np.max(abs(diff))<1e-12
    region=np.linalg.norm(x[ids],axis=1)<.5
    rows.append(dict(fraction=fraction,source_time=t,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        cells=len(ids),maximum_p_D_Pi_reintegration_difference=float(np.max(abs(diff))),
        global_momentum_difference=diff[:,:2].sum(0).tolist(),global_angular_momentum_difference=float(np.sum(diff[:,3]-diff[:,4])),
        global_kinetic_energy_difference=float(.5*np.sum(diff[:,6]+diff[:,8])),
        original_internal_residual_kinetic=float(.5*np.sum(L*L)),remaining_independent_residual_kinetic=float(sum(residual_K)),
        inner_generator_velocity_change_mass_RMS=float(np.sqrt(np.average(np.sum(change[region]**2,axis=1),weights=mass[region]))),
        maximum_generator_velocity_change=float(np.max(np.linalg.norm(change,axis=1))),diagnostic=diag))
print(json.dumps(dict(scope=__doc__,results=rows,
    limits='Uses reconstructed positive density, projected curvature, and original-quadrature covariance factors. No original checkpoint is rewritten. Changes in generator velocity are measured, not judged as accuracy improvements. No evolution, gravity closure, or zero-weight reduction proof.'),indent=2))
