"""Quantify uncapped failed face constraints by minimum uniform relaxation."""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from subcell_density_audit import triangles
from moment_aligned_density_audit import aligned_parts
from face_moment_feasibility_audit import build,face_matrix,OPTIONS


def minimax(p,F,target):
    k=p['k'];Aub=sparse.vstack([sparse.hstack([sparse.csr_matrix(F),-np.ones((len(target),1))]),
        sparse.hstack([sparse.csr_matrix(-F),-np.ones((len(target),1))])],format='csr');bub=np.r_[target,-target]
    Aeq=sparse.hstack([sparse.csr_matrix(p['A']),sparse.csr_matrix((len(p['b']),1))],format='csr')
    c=np.r_[np.zeros(k),1.];bounds=p['bounds'][:k]+[(0,None)]
    sol=linprog(c,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=p['b'],bounds=bounds,method='highs',options=OPTIONS)
    if not sol.success:return dict(gate=False,reason=sol.message)
    u=sol.x[:k];actual=float(np.max(abs(F@u-target)))
    moment=float(np.max(abs(p['A']@u-p['b'])))
    bound=max(0.,max(lo-x for (lo,hi),x in zip(bounds[:k],u)),max(x-hi for (lo,hi),x in zip(bounds[:k],u)))
    eqdual=sol.eqlin.marginals;ubdual=sol.ineqlin.marginals;lower=sol.lower.marginals;upper=sol.upper.marginals
    lo=np.array([a for a,b in bounds]);hi=np.array([np.inf if b is None else b for a,b in bounds]);finite=np.isfinite(hi)
    dual=float(p['b']@eqdual+bub@ubdual+lo@lower+hi[finite]@upper[finite])
    stationarity=float(np.max(abs(c-Aeq.T@eqdual-Aub.T@ubdual-lower-upper)))
    gap=abs(float(sol.fun)-dual);sign_error=float(max(0.,np.max(ubdual),-np.min(lower),np.max(upper)))
    return dict(gate=moment<1e-8 and bound<1e-8 and abs(actual-sol.fun)<1e-8 and gap<1e-8 and stationarity<1e-8 and sign_error<1e-8,
        minimum_normalized_residual=actual,dual_objective=dual,duality_gap=gap,dual_stationarity_error=stationarity,
        moment_error=moment,density_box_error_normalized=bound,solver_objective=float(sol.fun),dual_sign_error=sign_error)


def main():
    report=json.loads(Path('completion_report.json').read_text());previous=json.loads(Path('predecessor_completion.json').read_text());cases={}
    for name,result in report['cases'].items():
        cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config'];setup('shear_contact',cfg['n'],cfg['amplitude'])
        mass=np.asarray(cp['mass']);a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
        means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)]);low=float(means.min());high=float(means.max())
        centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in cells]);neighbors=[set() for _ in mass]
        faces=[[] for _ in mass]
        for edge_id,(i,j,va,vb,n,length) in enumerate(edges):
            faces[i].append((edge_id,va,vb))
            if j>=0:faces[j].append((edge_id,va,vb));neighbors[i].add(j);neighbors[j].add(i)
        neighbor_priors=[interface_prior(i,v,centers,means,neighbors)[1] for i,v in enumerate(cells)]
        moment_priors=previous['cases'][name]['variants']['moment_uniform']['inferred_priors'];rows=[]
        for original in result['rows']:
            mode='mean' if not original['mean_only']['feasible'] else 'joint' if not original['mean_and_first']['feasible'] else None
            if mode is None:continue
            i=original['cell'];v=cells[i];level=original['level'];tag=original['variant']
            metadata=(moment_priors if tag.startswith('moment') else neighbor_priors)[i]
            parts=aligned_parts(v,level,metadata) if tag.endswith('aligned') else triangles(v,level)
            p=build(v,np.r_[mass[i],a[i,3:9]],parts,metadata,means[i],low,high)
            T0=[];T1=[]
            if original['face_ids']!=[f[0] for f in faces[i]]:raise RuntimeError('face ordering mismatch')
            for edge_id,va,vb in faces[i]:
                w0,w1=face_matrix(parts,va,vb);T0.append(w0);T1.append(w1)
            target0=np.asarray(original['exact_face_means'])/high
            target1=np.asarray(original['exact_first_moments'])/high
            F=np.asarray(T0) if mode=='mean' else np.vstack([T0,T1])
            target=target0 if mode=='mean' else np.r_[target0,target1]
            solved=minimax(p,F,target)
            solved.update(cell=i,variant=tag,level=level,mode=mode,physical_density_scale=high)
            rows.append(solved)
        cases[name]=rows;print(name,'relaxation LPs',len(rows),flush=True)
    all_rows=[r for rows in cases.values() for r in rows]
    atomic('residual_report.json',dict(completed=True,all_gates=all(x['gate'] for x in all_rows),cases=cases,
        note='Minimum max absolute face-row residual normalized by H. Mean infeasibility implies joint infeasibility, avoiding duplicate LPs. Fixed m,q,Q and density box; no volume cap. Numerical primal/dual diagnostics, not exact-arithmetic certificates.'))
    raise SystemExit(0 if all(x['gate'] for x in all_rows) else 2)


if __name__=='__main__':main()
