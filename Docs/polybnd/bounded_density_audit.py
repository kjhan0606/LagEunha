"""Bounded moment reconstruction using numerical density extrema only."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from subcell_density_audit import fit,triangles,triangle_moments,normalized
from adaptive_coupled_mass import atomic


def minimum_relaxation(vertices,target,level,low,high):
    """Smallest eps in [0,1] for low*(1-eps)<=rho<=high*(1+eps).

    A finite-basis feasibility result, not a continuum existence theorem.
    """
    raw=np.array([triangle_moments(v) for v in triangles(vertices,level)])
    geo=raw.sum(0);center=geo[1:3]/geo[0];scale=np.sqrt(geo[0]);mass=target[0]
    A=np.column_stack([normalized(r,center,scale)/r[0] for r in raw])
    b=normalized(target,center,scale)/mass;k=len(raw);I=np.eye(k)
    lower=low*raw[:,0]/mass;upper=high*raw[:,0]/mass
    inequalities=np.vstack([np.column_stack([-I,-lower]),np.column_stack([I,-upper])])
    sol=linprog(np.r_[np.zeros(k),1.],A_ub=inequalities,b_ub=np.r_[-lower,upper],
        A_eq=np.column_stack([A,np.zeros(6)]),b_eq=b,bounds=[(0,None)]*k+[(0,1)],
        method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not sol.success:return dict(feasible=False,reason=sol.message)
    z=sol.x[:k];eps=float(sol.x[-1]);moment_error=float(np.max(abs(A@z-b)))
    bound_error=float(max(0.,np.max(lower*(1-eps)-z),np.max(z-upper*(1+eps))))
    return dict(feasible=True,minimum_relaxation=eps,moment_error=moment_error,
        fraction_bound_error=bound_error,gate=moment_error<1e-8 and bound_error<1e-8)


def audit(name):
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    if cfg['kind']!='shear_contact' or cp['t']!=.01:raise RuntimeError('configuration mismatch')
    rebuilt,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    a0=y0[:-4].reshape(-1,20);initial,_=model.geometry(a0[:,:2],a0[:,2],0.)
    refs=[];phases=np.where(a0[:,0]<.5,1.,cfg['amplitude'])
    for poly in initial:
        v=poly.copy();v[:,0]+=cp['t']*(.2+1.1*v[:,1]);refs.append(Polygon(v))
    a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
    centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in cells])
    means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)])
    low=float(means.min());high=float(means.max());neighbors=[set() for _ in mass]
    for i,j,*_ in edges:
        if j>=0:neighbors[i].add(j);neighbors[j].add(i)
    priors=[interface_prior(i,v,centers,means,neighbors) for i,v in enumerate(cells)]
    smooth=[model.Density(v,m,a[i,3:5],a[i,5:9].reshape(2,2)) for i,(v,m) in enumerate(zip(cells,mass))]
    relaxation=[minimum_relaxation(v,np.r_[mass[i],a[i,3:9]],2,low,high) for i,v in enumerate(cells)]
    rows=[]
    for margin in (0.,.001,.01,.05,.1,.5):
        bounds=(low*(1-margin),high*(1+margin))
        for blend in (.5,1.):
            records=[];error=0.
            for i,v in enumerate(cells):
                seed=lambda x:(1-blend)*smooth[i].value(x)+blend*priors[i][0](x)
                record,parts,rho=fit(v,np.r_[mass[i],a[i,3:9]],2,seed,density_bounds=bounds)
                record['cell']=i;records.append(record)
                if rho is None:continue
                for triangle,density in zip(parts,rho):
                    cell=Polygon(triangle)
                    for ref,exact in zip(refs,phases):error+=cell.intersection(ref).area*abs(density-exact)
            feasible=all(r['feasible'] for r in records)
            rows.append(dict(margin=margin,blend=blend,bounds=bounds,cells=records,all_feasible=feasible,
                all_gates=feasible and all(r['moment_gate'] and r['bound_gate'] for r in records),
                density_L1_relative=error/float(mass.sum()) if feasible else None))
    result=dict(case=name,status='completed',inferred_phase_densities=[low,high],minimum_relaxation=relaxation,
        comparisons=rows,note='Fixed numerical m,q,Q. Bounds inferred from current mean-density extrema; no exact phase bounds in fitting. Infeasible fits retained. Finite level-2 basis only.')
    atomic(name+'.json',result);return result


def main():
    if not json.loads(Path('predecessor_completion.json').read_text())['completed']:raise RuntimeError('predecessor incomplete')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(status='failed',error=repr(exc))
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',result);print(json.dumps(result),flush=True)
    raise SystemExit(0 if result['completed'] else 2)


if __name__=='__main__':main()
