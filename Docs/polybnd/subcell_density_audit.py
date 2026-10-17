"""Positive piecewise constant densities constrained by all six moments.
Linear programs use numerical moments only; exact density is for scoring.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from shear_interface_audit import moments
from adaptive_coupled_mass import atomic


def triangles(vertices,level):
    center=vertices.mean(0)
    result=[np.array([center,a,b]) for a,b in zip(vertices,np.roll(vertices,-1,axis=0))]
    for _ in range(level):
        refined=[]
        for a,b,c in result:
            ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
            refined.extend([np.array([a,ab,ca]),np.array([ab,b,bc]),np.array([ca,bc,c]),np.array([ab,bc,ca])])
        result=refined
    return result


def triangle_moments(v):
    a,b=v[1]-v[0],v[2]-v[0];area=abs(a[0]*b[1]-a[1]*b[0])/2
    s=v.sum(0);Q=area*(v.T@v+np.outer(s,s))/12
    return np.r_[area,area*s/3,Q.ravel()]


def normalized(raw,center,scale):
    m=raw[0];q=raw[1:3];Q=raw[3:7].reshape(2,2)
    S=(Q-np.outer(center,q)-np.outer(q,center)+m*np.outer(center,center))/scale**2
    return np.r_[m,(q-m*center)/scale,S[0,0],S[0,1],S[1,1]]


def fit(vertices,target,level,seed,density_bounds=None,part_averages=None):
    parts=triangles(vertices,level);raw=np.array([triangle_moments(v) for v in parts])
    geo=raw.sum(0);center=geo[1:3]/geo[0];scale=np.sqrt(geo[0]);mass=target[0]
    # Unknowns are subcell mass fractions, plus absolute deviations from seed.
    A=np.column_stack([normalized(r,center,scale)/r[0] for r in raw])
    b=normalized(target,center,scale)/mass
    values=np.asarray(seed(np.array([v.mean(0) for v in parts])) if part_averages is None else part_averages(parts))
    if values.shape!=(len(parts),) or not np.all(np.isfinite(values)) or np.any(values<0):
        raise ValueError('invalid prior density averages')
    prior=values*raw[:,0]/mass
    k=len(parts);I=np.eye(k)
    bounds=(0,None)
    if density_bounds is not None:
        lower,upper=density_bounds
        if not 0<=lower<=upper:raise ValueError('invalid density bounds')
        bounds=[(lower*r[0]/mass,upper*r[0]/mass) for r in raw]+[(0,None)]*k
    sol=linprog(np.r_[np.zeros(k),np.ones(k)],A_ub=np.block([[I,-I],[-I,-I]]),
        b_ub=np.r_[prior,-prior],A_eq=np.column_stack([A,np.zeros_like(A)]),b_eq=b,
        bounds=bounds,method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not sol.success:return dict(feasible=False,reason=sol.message),parts,None
    rho=sol.x[:k]*mass/raw[:,0];residual=float(np.max(abs(A@sol.x[:k]-b)))
    reconstructed=rho@raw
    bound_error=0. if density_bounds is None else max(0.,lower-float(rho.min()),float(rho.max())-upper)/max(1.,upper)
    return dict(feasible=True,max_scaled_moment_error=residual,moment_gate=residual<1e-8,
        prior_mass_relative_error=float(abs(prior.sum()-1)),
        bound_gate=bound_error<1e-8,max_scaled_bound_error=bound_error,
        mass_relative_error=float(abs(reconstructed[0]-mass)/mass),
        min_density=float(rho.min()),max_density=float(rho.max()),subcells=k),parts,rho


def audit(name):
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    if cfg['kind']!='shear_contact' or cp['t']!=.01:raise RuntimeError('configuration mismatch')
    mass,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);a0=y0[:-4].reshape(-1,20)
    saved_mass=np.asarray(cp['mass'])
    initialization_difference=float(np.max(abs(saved_mass-mass)))
    # Cross-host geometry may differ by a few ulps; the saved masses remain authoritative.
    np.testing.assert_allclose(saved_mass,mass,rtol=1e-13,atol=1e-15)
    mass=saved_mass
    initial,_=model.geometry(a0[:,:2],a0[:,2],0.);high=cfg['amplitude']
    phase=np.where(a0[:,0]<.5,1.,high);refs=[]
    for p in initial:
        v=p.copy();v[:,0]+=cp['t']*(.2+1.1*v[:,1]);refs.append(Polygon(v))
    a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,_=model.geometry(a[:,:2],a[:,2],cp['t'])
    rows=[];bounds=[]
    for i,p in enumerate(cells):
        mean=mass[i]/Polygon(p).area
        if mean<1 or mean>high:bounds.append(dict(cell=i,mean_density=mean,below_low=max(0.,1-mean),above_high=max(0.,mean-high)))
    for level in (1,2):
        records=[];error=0.
        for i,p in enumerate(cells):
            target=np.r_[mass[i],a[i,3:9]]
            smooth=model.Density(p,mass[i],a[i,3:5],a[i,5:9].reshape(2,2))
            record,parts,rho=fit(p,target,level,smooth.value);record['cell']=i;records.append(record)
            if rho is None:continue
            for vertices,density in zip(parts,rho):
                cell=Polygon(vertices)
                for ref,exact_density in zip(refs,phase):error+=cell.intersection(ref).area*abs(density-exact_density)
        feasible=all(r['feasible'] for r in records)
        rows.append(dict(level=level,cells=records,all_feasible=feasible,
            all_moments_preserved=feasible and all(r['moment_gate'] for r in records),
            density_L1_relative=error/float(mass.sum()) if feasible else None))
    result=dict(case=name,status='completed',phase_bound_violations=bounds,levels=rows,
        initial_mass_rebuild_max_difference=initialization_difference,
        note='Nonnegative subcell densities, no upper cap. Minimizes mass-weighted L1 distance to sampled smooth reconstruction. Feasibility is not accuracy; exact interface is not supplied to fit. No evolution modified.')
    atomic(name+'.json',result);return result


def main():
    if not json.loads(Path('predecessor_completion.json').read_text())['completed']:raise RuntimeError('predecessor incomplete')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(status='failed',error=repr(exc))
    r=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',r);print(json.dumps(r),flush=True)
    raise SystemExit(0 if r['completed'] else 2)


if __name__=='__main__':main()
