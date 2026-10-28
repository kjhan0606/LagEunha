"""Constructive positive density after actual full-moment material drift.
Uses prior quadratic fit and positive Bernstein fallback without resetting
any mass/q/Q. No gravity integration or Euler-evolution claim.
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

def main():
    root=Path(__file__).parent;assert not (root/'result.json').exists();base=Path('/gpfs/kjhan/LagForce');started=time.monotonic()
    modpath=base/'Codex_weight_remap_moments_centered_20261017/audit.py';geopath=base/'Codex_weighted_compression_geometry_20261017/audit.py';fitpath=base/'Codex_compact_reconstruction_20261018/audit.py';bernpath=base/'Codex_positive_density_fallback_20261018/audit.py'
    mod=load('moments',modpath);geo=load('geo',geopath);fit=load('fit',fitpath);bern=load('bern',bernpath)
    targetpath=base/'Codex_evolving_weight_transport_20261017/weights_0.npz';target=np.load(targetpath)['weights'];paths=[Path(__file__),modpath,geopath,fitpath,bernpath,targetpath];reports=[]
    for fraction in (.25,1.):
        for dt in (.005,.02):
            source=base/f'Codex_material_drift_support_20261018/drift_w{fraction:g}_dt{dt:g}.npz';paths.append(source)
            saved=np.load(source);ids=saved['ids'];q=saved['hydro_moments'];polys=mod.polygons(saved['sites'],fraction*target,ids,geo.clip)
            minimum=[];error=[];fallback=[];failures=[];representations=[]
            for i,p,z in zip(ids,polys,q):
                m=z[0];origin=p.mean(axis=0);V=mod.area(p);h=np.sqrt(V);poly=(p-origin)/h
                c=z[1:3]/m;Q=np.array([[z[3],z[4]],[z[4],z[5]]])/m
                C=(Q-np.outer(origin,c)-np.outer(c,origin)+np.outer(origin,origin))/(h*h)
                b=np.r_[1.,(c-origin)/h,C[0,0],C[0,1],C[1,1]]
                xy,dv=mod.quadrature(poly,3);F=fit.basis(xy);G=F.T@(dv[:,None]*F)/dv.sum();theta=np.linalg.solve(G,b)
                lo=fit.minimum_quadratic(poly,theta)
                if lo>1e-10:
                    minimum.append(lo);error.append(float(np.max(abs(G@theta-b))))
                    representations.append({'id':int(i),'kind':'quadratic','theta':theta.tolist(),'origin':origin.tolist(),'h':h,'mass':m});continue
                attempts=[]
                for degree in (2,4,6):
                    features,reference,_=bern.components(poly,degree);N=len(reference)
                    eq=np.column_stack((features,np.zeros(6)));ub=hstack((-eye(N,format='csr'),csr_matrix(reference[:,None])),format='csr')
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore',message='Unrecognized options detected')
                        sol=linprog(np.r_[np.zeros(N),-1.],A_ub=ub,b_ub=np.zeros(N),A_eq=eq,b_eq=b,bounds=[(0,None)]*N+[(0,1)],method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9,'threads':1})
                    attempts.append({'degree':degree,'status':int(sol.status)})
                    if not sol.success:continue
                    beta=sol.x[:-1];lo=float(np.min(beta/reference));res=float(np.max(abs(features@beta-b)))
                    assert lo>=-1e-10 and res<2e-8
                    minimum.append(lo);error.append(res);fallback.append({'id':int(i),'degree':degree,'lower_bound':lo,'residual':res,'attempts':attempts})
                    representations.append({'id':int(i),'kind':'bernstein','degree':degree,'beta':beta.tolist(),'reference':reference.tolist(),'origin':origin.tolist(),'h':h,'mass':m,'polygon':poly.tolist()});break
                else:failures.append({'id':int(i),'attempts':attempts})
            report={'weight_fraction':fraction,'dt':dt,'cells':len(ids),'positive_reconstructions':len(minimum),'minimum_density_factor':min(minimum),'moment_residual_max':max(error),'fallback':fallback,'failed':failures}
            reports.append(report);print(json.dumps({k:v for k,v in report.items() if k!='fallback'}),flush=True)
            (root/f'density_w{fraction:g}_dt{dt:g}.json').write_text(json.dumps(representations)+'\n')
    result={'scope':'Actual material-drifted central1590cell states in reconstructed moving power cells at dt .005,.02. Static end-state density reconstruction, no complete flow evolution.',
            'reports':reports,'seconds':time.monotonic()-started,
            'limits':['Matching m,q,Q and positive density does not establish face flux consistency or unique velocity reconstruction.', 'Residual velocity covariance, pressure/gravity work and zero-weight model reduction must be coupled before adoption.', 'Density LP is CPU/Python; only the preceding moment drift is shared C/OpenMP/CUDA.'],
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (root/'result.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
