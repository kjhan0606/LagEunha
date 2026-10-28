"""Positive polynomial density fallback for actual remapped cell moments.

Each centroid triangle uses a Bernstein polynomial of degree n. Positive
coefficients imply positive density over the entire triangle. Component
integrals of 1,x,y,xx,xy,yy are analytic Dirichlet moments; no gravity
quadrature. A local LP maximizes a common density lower bound. Diagnostic
Python reconstruction only; not yet a production or CUDA density solver.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import warnings
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import eye,csr_matrix,hstack

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def components(poly,degree):
    center=poly.mean(axis=0);features=[];reference=[];metadata=[]
    count=(degree+1)*(degree+2)//2
    for aa,bb in zip(poly,np.roll(poly,-1,axis=0)):
        vertices=np.array([center,aa,bb]);d=aa-center;e=bb-center
        area=.5*(d[0]*e[1]-d[1]*e[0]);assert area>0
        for i in range(degree+1):
            for j in range(degree-i+1):
                powers=np.array([i,j,degree-i-j]);alpha=powers+1.;total=alpha.sum()
                mean=alpha@vertices/total
                second=vertices.T@(np.outer(alpha,alpha)+np.diag(alpha))@vertices/(total*(total+1))
                features.append([1.,*mean,second[0,0],second[0,1],second[1,1]])
                reference.append(area/count);metadata.append((vertices,powers))
    reference=np.array(reference);reference/=reference.sum()
    return np.array(features).T,reference,metadata

def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists();started=time.monotonic()
    base=Path('/gpfs/kjhan/LagForce');prior=base/'Codex_compact_reconstruction_20261018'
    previous=json.loads((prior/'result.json').read_text())
    modpath=base/'Codex_weight_remap_moments_centered_20261017/audit.py';geopath=base/'Codex_weighted_compression_geometry_20261017/audit.py';readerpath=base/'Codex_paired_response_20261017/run.py'
    mod=load('moments',modpath);geo=load('geo',geopath);reader=load('reader',readerpath)
    source=base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat';layout=base/'Codex_paired_response_20261017/layout.json'
    _,a,sort,_=reader.read(source,json.loads(layout.read_text()));a=a[sort];x=np.column_stack((a['x']-3,a['y']-3))
    targetpath=base/'Codex_evolving_weight_transport_20261017/weights_0.npz';target=np.load(targetpath)['weights']
    rows=[];paths=[Path(__file__),prior/'result.json',modpath,geopath,readerpath,source,layout,targetpath]
    for old in previous['rows']:
        fraction=old['weight_fraction'];ids=np.array(old['negative_density_ids'],dtype=int)
        statepath=base/f'Codex_checkpoint_weight_remap_20261018/state_w{fraction:g}_q8.npz';paths.append(statepath)
        data=np.load(statepath);index={int(i):k for k,i in enumerate(data['ids'])};polygons=mod.polygons(x,fraction*target,ids,geo.clip)
        for cell,poly in zip(ids,polygons):
            q=data['state'][index[int(cell)]];m=q[0];origin=poly.mean(axis=0);h=np.sqrt(mod.area(poly));p=(poly-origin)/h
            c=q[1:3]/m;Q=np.array([[q[3],q[4]],[q[4],q[5]]])/m
            centered=(Q-np.outer(origin,c)-np.outer(c,origin)+np.outer(origin,origin))/(h*h)
            b=np.r_[1.,(c-origin)/h,centered[0,0],centered[0,1],centered[1,1]]
            trials=[]
            for degree in (2,4,6):
                F,ref,metadata=components(p,degree);N=len(ref)
                eq=np.column_stack((F,np.zeros(6)));ub=hstack((-eye(N,format='csr'),csr_matrix(ref[:,None])),format='csr')
                objective=np.r_[np.zeros(N),-1.]
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore',message='Unrecognized options detected')
                    sol=linprog(objective,A_ub=ub,b_ub=np.zeros(N),A_eq=eq,b_eq=b,bounds=[(0,None)]*N+[(0,1)],method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9,'threads':1})
                trial={'degree':degree,'components':N,'success':bool(sol.success),'status':int(sol.status)};trials.append(trial)
                if not sol.success:continue
                beta=sol.x[:-1];res=float(np.max(abs(F@beta-b)));lower=float(np.min(beta/ref))
                trial.update({'moment_residual':res,'whole_cell_density_lower_bound_over_mean':lower,'lp_bound':float(sol.x[-1])})
                assert res<2e-8 and lower>=-1e-10
                # Independently integrate every polynomial basis of one triangle
                # with a sufficient Gauss rule; verifies analytic moment formula.
                import math
                check=0.
                for k,(vertices,powers) in enumerate(metadata[:(degree+1)*(degree+2)//2]):
                    xy,dv=mod.quadrature(vertices,degree+2)
                    bary=np.column_stack((xy,np.ones(len(xy))))@np.linalg.inv(np.column_stack((vertices,np.ones(3))))
                    bern=math.factorial(degree)/np.prod([math.factorial(int(z)) for z in powers])*np.prod(bary**powers,axis=1)
                    g=np.column_stack((np.ones(len(xy)),xy,xy[:,0]**2,xy[:,0]*xy[:,1],xy[:,1]**2))
                    observed=(dv*bern)@g/np.sum(dv*bern)
                    check=max(check,float(np.max(abs(observed-F[:,k]))))
                trial['analytic_basis_moment_independent_error']=check;assert check<1e-10
                np.savez(root/f'density_w{fraction:g}_cell{cell}.npz',beta=beta,reference=ref,origin=origin,h=h,degree=degree,polygon=p,mass=m)
                break
            rows.append({'weight_fraction':fraction,'id':int(cell),'trials':trials});print(json.dumps(rows[-1]),flush=True)
    result={'scope':'Only cells where the unconstrained quadratic density was negative, from actual t2 quarter/full finite weight remaps. No physical time evolution.',
            'representation':'Piecewise triangle Bernstein density with nonnegative coefficients. Match m,q,Q using analytic moments; maximize strictly positive lower bound. Density can jump at internal triangle boundaries.',
            'rows':rows,'seconds':time.monotonic()-started,
            'limits':['Local static feasibility, not a complete pressure/gravity/velocity covariance closure.', 'Does not match nonlinear potential moment.', 'LP is CPU/Python only; not claimed as CUDA reconstruction or performance-optimized production implementation.', 'Bounds use double precision with reported linear residual; not interval arithmetic.'],
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
