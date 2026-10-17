"""Neighbor-directed interface priors with unchanged mass/moment constraints."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from subcell_density_audit import fit
from discontinuous_density_audit import clip
from adaptive_coupled_mass import atomic


def phase_floor_threshold(low,high,fraction_floor=0.):
    return max(32*np.finfo(float).eps*max(1.,abs(low),abs(high)),fraction_floor*(high-low))


def phase_floor_blend(mean,low,high,fraction_floor=0.):
    threshold=phase_floor_threshold(low,high,fraction_floor)
    s=float(np.clip((min(mean-low,high-mean)/threshold-1)/3,0.,1.))
    return s*s*(3-2*s)


def interface_prior(i,vertices,centers,means,neighbors,smooth_phase_floor=False,phase_fraction_floor=0.):
    low=float(means.min());high=float(means.max());area=Polygon(vertices).area
    if high-low<1e-12:return lambda x:np.full(len(x),means[i]),dict(flat=True)
    # Roundoff in a uniform phase must not create a microscopic corner cut
    # whose face trace has a much larger error than its volume fraction.
    phase_roundoff=phase_floor_threshold(low,high,phase_fraction_floor if smooth_phase_floor else 0.)
    if min(means[i]-low,high-means[i])<=phase_roundoff:
        return lambda x:np.full(len(x),means[i]),dict(flat=True)
    ids=sorted(neighbors[i]);delta=centers[ids]-centers[i]
    distance=np.linalg.norm(delta,axis=1);weight=1/np.maximum(distance,1e-12)
    gradient=np.linalg.lstsq(delta*weight[:,None],(means[ids]-means[i])*weight,rcond=None)[0]
    norm=np.linalg.norm(gradient)
    if norm<1e-12:return lambda x:np.full(len(x),means[i]),dict(flat=True)
    normal=gradient/norm;fraction=(means[i]-low)/(high-low)
    if fraction<=0:return lambda x:np.full(len(x),low),dict(flat=True)
    if fraction>=1:return lambda x:np.full(len(x),high),dict(flat=True)
    lo,hi=float((vertices@normal).min()),float((vertices@normal).max())
    for _ in range(48):
        mid=(lo+hi)/2
        if clip(vertices,normal,mid).area>fraction*area:lo=mid
        else:hi=mid
    offset=(lo+hi)/2
    metadata=dict(flat=False,normal=normal.tolist(),offset=offset,fraction=float(fraction))
    blend=phase_floor_blend(means[i],low,high,phase_fraction_floor) if smooth_phase_floor else 1.
    if blend<1:
        metadata['blend']=blend
        return lambda x:(1-blend)*means[i]+blend*np.where(x@normal>=offset,high,low),metadata
    return lambda x:np.where(x@normal>=offset,high,low),metadata


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
    neighbors=[set() for _ in mass]
    for i,j,*_ in edges:
        if j>=0:neighbors[i].add(j);neighbors[j].add(i)
    priors=[interface_prior(i,v,centers,means,neighbors) for i,v in enumerate(cells)]
    smooth=[model.Density(v,m,a[i,3:5],a[i,5:9].reshape(2,2)) for i,(v,m) in enumerate(zip(cells,mass))]
    rows=[]
    for blend in (.25,.5,1.):
        records=[];error=0.
        for i,v in enumerate(cells):
            seed=lambda x:(1-blend)*smooth[i].value(x)+blend*priors[i][0](x)
            record,parts,rho=fit(v,np.r_[mass[i],a[i,3:9]],2,seed);record['cell']=i;records.append(record)
            if rho is None:continue
            for triangle,density in zip(parts,rho):
                cell=Polygon(triangle)
                for ref,exact in zip(refs,phases):error+=cell.intersection(ref).area*abs(density-exact)
        feasible=all(r['feasible'] for r in records)
        rows.append(dict(blend=blend,cells=records,all_feasible=feasible,
            all_moments_preserved=feasible and all(r['moment_gate'] for r in records),
            density_L1_relative=error/float(mass.sum()) if feasible else None))
    baseline=json.loads(Path('predecessor_completion.json').read_text())['cases'][name]['levels'][-1]
    result=dict(case=name,status='completed',inferred_phase_densities=[float(means.min()),float(means.max())],
        priors=[p[1] for p in priors],baseline_L1=baseline['density_L1_relative'],blends=rows,
        note='Neighbor gradients and current mean-density extrema only; no exact interface in optimization. Blend values preselected, no general optimality claim. Unchanged numerical masses and moments.')
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
