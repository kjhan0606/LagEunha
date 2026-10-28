"""Recover positive residual factors from original t2 remap quadrature atoms.
Preserves original checkpoints. No eigenvalue clipping or thermal transfer.
"""
import ctypes,hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');source=base/'Codex_checkpoint_weight_remap_20261018'
output=base/'Codex_checkpoint_covariance_factors_20261018';output.mkdir(exist_ok=False)
sys.path.insert(0,str(base))
from run_canonical_av_branches_20261018 import read
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
mod=load('factor_remap',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
geo=load('factor_geometry',base/'Codex_weighted_compression_geometry_20261017/audit.py')
ic,t=read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
x=np.column_stack([ic['x']-3,ic['y']-3]);vel=np.column_stack([ic['vx'],ic['vy']])
lift=np.load(source/'lift_coefficients.npz');ids=lift['ids'];lam=lift['coefficient'];coeff=lam/lift['feature_scale'];means=lift['cell_feature_mean']
old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);volume=lift['volume'];rho=ic['mass'][ids]/volume;ud=ic['ie'][ids]/volume
lo=np.array([p.min(0) for p in old]);hi=np.array([p.max(0) for p in old]);phi0=mod.potential(x[ids])
weights=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
cpu=ctypes.CDLL(str(source/'cpu.so'));arr=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
cpu.wm_lift_batch.argtypes=[ctypes.c_size_t,arr,arr];cpu.wm_lift_batch.restype=ctypes.c_int
assert cpu.wm_set_threads(1)==1
reports=[]
for fraction in (.25,1.):
    started=time.monotonic();path=source/f'state_w{fraction:g}_q8.npz';saved=np.load(path)['state']
    polygons=mod.polygons(x,fraction*weights,ids,geo.clip)
    recovered=[];factors=[];centers=[];means_v=[];gradients=[];minimum_dm=float('inf');atoms=0
    reconstruction_error=0.;orthogonality=0.;coverage=[]
    for dest,poly in enumerate(polygons):
        candidate=np.flatnonzero(np.all(hi>=poly.min(0),axis=1)&np.all(lo<=poly.max(0),axis=1))
        blocks=[];points=[];velocities=[];total_area=0.
        for k in candidate:
            overlap=old[k].copy()
            for a,b in zip(poly,np.roll(poly,-1,axis=0)):
                e=b-a;n=np.array([e[1],-e[0]]);overlap=geo.clip(overlap,n,n@a)
                if len(overlap)<3:break
            if mod.area(overlap)<=1e-18:continue
            xy,dv=mod.quadrature(overlap,8);N=len(dv);total_area+=dv.sum()
            block=np.column_stack([dv,xy,np.tile(x[ids[k]],(N,1)),np.tile(vel[ids[k]],(N,1)),
                np.full(N,rho[k]),np.full(N,ud[k]),mod.potential(xy),np.full(N,phi0[k]),np.tile(coeff,(N,1)),np.full(N,means[k]@lam)])
            blocks.append(block);points.append(xy);velocities.append(np.tile(vel[ids[k]],(N,1)))
        data=np.ascontiguousarray(np.concatenate(blocks));xy=np.concatenate(points);vv=np.concatenate(velocities)
        one=np.empty((len(data),17));assert cpu.wm_lift_batch(len(data),data,one)==0
        dm=one[:,0];assert np.all(dm>0);minimum_dm=min(minimum_dm,float(dm.min()));atoms+=len(dm)
        q=one.sum(0);m=dm.sum();c=(dm[:,None]*xy).sum(0)/m;v=(dm[:,None]*vv).sum(0)/m
        X=xy-c;Y=vv-v;B=np.linalg.lstsq(np.sqrt(dm)[:,None]*X,np.sqrt(dm)[:,None]*Y,rcond=None)[0].T
        residual=Y-X@B.T;_,upper=np.linalg.qr(np.sqrt(dm)[:,None]*residual,mode='reduced');L=upper.T
        S=X.T@(dm[:,None]*X);Pi=np.outer(v,v)*m+B@S@B.T+L@L.T
        target=np.array([[q[12],q[13]],[q[13],q[14]]])
        reconstruction_error=max(reconstruction_error,float(np.max(abs(Pi-target))))
        orthogonality=max(orthogonality,float(np.max(abs(X.T@(dm[:,None]*residual)))),float(np.max(abs(dm@residual))))
        coverage.append(abs(total_area/mod.area(poly)-1));recovered.append(q);factors.append(L);centers.append(c);means_v.append(v);gradients.append(B)
    recovered=np.array(recovered);factors=np.array(factors)
    difference=recovered-saved
    assert max(coverage)<1e-8 and np.max(abs(difference))<1e-12 and reconstruction_error<1e-12
    row=dict(fraction=fraction,source_time=t,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),cells=len(ids),atoms=atoms,
        minimum_atom_mass=minimum_dm,max_coverage_error=max(coverage),max_reintegrated_moment_difference=float(np.max(abs(difference))),
        kinetic_tensor_reconstruction_error=reconstruction_error,residual_orthogonality_error=orthogonality,
        global_energy_difference=float(difference[:,15].sum()),global_angular_momentum_difference=float((difference[:,9]-difference[:,10]).sum()),
        factor_minimum_singular_value=float(np.linalg.svd(factors,compute_uv=False)[:,-1].min()),seconds=time.monotonic()-started)
    np.savez(output/f'factors_w{fraction:g}.npz',ids=ids,factors=factors,center=centers,mean_velocity=means_v,gradient=gradients,reintegrated_state=recovered,original_state=saved)
    reports.append(row);print(json.dumps(row),flush=True)
(output/'report.json').write_text(json.dumps(dict(scope=__doc__,results=reports,limits='Quadrature source recovery only, not evolution or automatic production state conversion. Reintegrated moments differ from stored sums by recorded rounding.'),indent=2))
