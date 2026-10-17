"""Exact spatial L1 objective to an inferred two-phase prior on fixed triangles."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from shapely.geometry import Polygon
from subcell_density_audit import triangles,triangle_moments,normalized
from area_prior_density_audit import interface_averages
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from density_sensitivity_audit import rotated_prior


def phase_fit(vertices,target,level,seed,density_bounds=None,part_averages=None):
    # The callable carries only the inferred prior, never the reference solution.
    parts=triangles(vertices,level);raw=np.array([triangle_moments(v) for v in parts])
    geo=raw.sum(0);center=geo[1:3]/geo[0];scale=np.sqrt(geo[0]);mass=target[0]
    A=np.column_stack([normalized(r,center,scale)/r[0] for r in raw]);b=normalized(target,center,scale)/mass
    low,high,fraction=part_averages.phase_data(parts)
    low=np.broadcast_to(low,len(parts));high=np.broadcast_to(high,len(parts));fraction=np.asarray(fraction)
    if np.any(low<0) or np.any(high<low) or np.any(fraction<0) or np.any(fraction>1):raise ValueError('invalid phases')
    lower=low*raw[:,0]/mass;upper=high*raw[:,0]/mass
    prior=(1-fraction)*lower+fraction*upper
    k=len(parts);I=np.eye(k)
    # Epigraph of (1-f)*|z-low*area/m| + f*|z-high*area/m|.
    Aub=np.vstack([np.column_stack([-I,-I]),np.column_stack([(1-2*fraction)[:,None]*I,-I]),np.column_stack([I,-I])])
    bub=np.r_[-prior,(1-fraction)*lower-fraction*upper,prior]
    if density_bounds is None:bounds=(0,None)
    else:
        lo,hi=density_bounds
        bounds=[(lo*r[0]/mass,hi*r[0]/mass) for r in raw]+[(0,None)]*k
    sol=linprog(np.r_[np.zeros(k),np.ones(k)],A_ub=Aub,b_ub=bub,
        A_eq=np.column_stack([A,np.zeros_like(A)]),b_eq=b,bounds=bounds,method='highs',
        options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    if not sol.success:return dict(feasible=False,reason=sol.message),parts,None
    z=sol.x[:k];rho=z*mass/raw[:,0];residual=float(np.max(abs(A@z-b)));reconstructed=rho@raw
    bound_error=0. if density_bounds is None else max(0.,lo-float(rho.min()),float(rho.max())-hi)/max(1.,hi)
    direct=float(np.sum((1-fraction)*abs(z-lower)+fraction*abs(z-upper)))
    objective_error=abs(direct-float(sol.fun))
    return dict(feasible=True,max_scaled_moment_error=residual,moment_gate=residual<1e-8,
        prior_mass_relative_error=float(abs(prior.sum()-1)),bound_gate=bound_error<1e-8,
        max_scaled_bound_error=bound_error,mass_relative_error=float(abs(reconstructed[0]-mass)/mass),
        min_density=float(rho.min()),max_density=float(rho.max()),subcells=k,
        inferred_prior_L1_relative=direct,objective_error=objective_error),parts,rho



def wrapped_fit(vertices,target,level,seed,density_bounds=None,part_averages=None):
    # Driver supplies a closure with prior data; keep exact field out of LP.
    def phase_data(parts):
        averages=part_averages(parts)
        prior_seed,metadata,low,high=part_averages.prior_data
        if metadata['flat']:
            value=float(prior_seed(np.array([vertices.mean(0)]))[0])
            return value,value,np.zeros(len(parts))
        return low,high,np.clip((averages-low)/(high-low),0.,1.)
    part_averages.phase_data=phase_data
    return phase_fit(vertices,target,level,seed,density_bounds,part_averages)


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
    base=[interface_prior(i,v,centers,means,neighbors) for i,v in enumerate(cells)]
    predecessor=json.loads(Path('predecessor_completion.json').read_text())['cases'][name]
    bounds=(low*.99,high*1.01);rows=[];relaxation={}
    for level in (1,2,3):
        relaxation[str(level)]=predecessor['minimum_relaxation'][str(level)]
        for angle in (-15.,0.,15.):
            records=[];error=0.
            for i,v in enumerate(cells):
                seed,metadata=rotated_prior(v,*base[i],low,high,angle)
                averages=lambda parts:interface_averages(parts,seed,metadata,low,high)
                averages.prior_data=(seed,metadata,low,high)
                record,parts,rho=wrapped_fit(v,np.r_[mass[i],a[i,3:9]],level,seed,density_bounds=bounds,part_averages=averages)
                record.update(cell=i,prior=metadata);records.append(record)
                if rho is None:continue
                for triangle,density in zip(parts,rho):
                    cell=Polygon(triangle)
                    for ref,exact in zip(refs,phases):error+=cell.intersection(ref).area*abs(density-exact)
            feasible=all(r['feasible'] for r in records)
            rows.append(dict(level=level,angle_degrees=angle,reused=False,cells=records,
                all_feasible=feasible,all_gates=feasible and all(r['moment_gate'] and r['bound_gate'] and r['prior_mass_relative_error']<1e-10 and r['objective_error']<1e-8 for r in records),
                area_average_L1_relative=next(r['density_L1_relative'] for r in predecessor['comparisons'] if r['level']==level and r['angle_degrees']==angle),
                density_L1_relative=error/float(mass.sum()) if feasible else None))
        atomic(name+'.partial.json',dict(comparisons=rows,minimum_relaxation=relaxation))
    result=dict(case=name,status='completed',bounds=bounds,minimum_relaxation=relaxation,comparisons=rows,
        note='Fixed numerical m,q,Q, 1% numerical mean-extrema bounds, full neighbor prior. Normal rotation preserves prior phase area fraction. Exact field used only for scoring. Objective is the exact spatial L1 to the inferred prior in the fixed triangle basis. All feasibility-relaxation results reused from predecessor; all reconstructions newly calculated.')
    atomic(name+'.json',result);return result


def main():
    if not json.loads(Path('predecessor_completion.json').read_text())['completed']:raise RuntimeError('predecessor incomplete')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(status='failed',error=repr(exc))
        print(name,cases[name]['status'],flush=True)
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',result)
    raise SystemExit(0 if result['completed'] else 2)


if __name__=='__main__':main()
