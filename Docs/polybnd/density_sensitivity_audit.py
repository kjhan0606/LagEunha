"""Finite-basis and neighbor-normal sensitivity of bounded reconstruction."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from discontinuous_density_audit import clip
from subcell_density_audit import fit
from bounded_density_audit import minimum_relaxation
from adaptive_coupled_mass import atomic


def rotated_prior(vertices,prior,metadata,low,high,angle):
    if metadata['flat'] or angle==0:return prior,metadata
    theta=np.deg2rad(angle);c,s=np.cos(theta),np.sin(theta)
    normal=np.array([[c,-s],[s,c]])@np.asarray(metadata['normal'])
    fraction=metadata['fraction'];area=Polygon(vertices).area
    lo,hi=float((vertices@normal).min()),float((vertices@normal).max())
    for _ in range(48):
        mid=(lo+hi)/2
        if clip(vertices,normal,mid).area>fraction*area:lo=mid
        else:hi=mid
    offset=(lo+hi)/2
    return lambda x:np.where(x@normal>=offset,high,low),dict(flat=False,normal=normal.tolist(),offset=offset,fraction=fraction)


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
        relaxation[str(level)]=predecessor['minimum_relaxation'] if level==2 else [
            minimum_relaxation(v,np.r_[mass[i],a[i,3:9]],level,low,high) for i,v in enumerate(cells)]
        for angle in (-15.,0.,15.):
            if level==2 and angle==0:
                row=dict(next(r for r in predecessor['comparisons'] if r['margin']==.01 and r['blend']==1.))
                row.update(level=level,angle_degrees=angle,reused=True);rows.append(row);continue
            records=[];error=0.
            for i,v in enumerate(cells):
                seed,metadata=rotated_prior(v,*base[i],low,high,angle)
                record,parts,rho=fit(v,np.r_[mass[i],a[i,3:9]],level,seed,density_bounds=bounds)
                record.update(cell=i,prior=metadata);records.append(record)
                if rho is None:continue
                for triangle,density in zip(parts,rho):
                    cell=Polygon(triangle)
                    for ref,exact in zip(refs,phases):error+=cell.intersection(ref).area*abs(density-exact)
            feasible=all(r['feasible'] for r in records)
            rows.append(dict(level=level,angle_degrees=angle,reused=False,cells=records,
                all_feasible=feasible,all_gates=feasible and all(r['moment_gate'] and r['bound_gate'] for r in records),
                density_L1_relative=error/float(mass.sum()) if feasible else None))
        atomic(name+'.partial.json',dict(comparisons=rows,minimum_relaxation=relaxation))
    result=dict(case=name,status='completed',bounds=bounds,minimum_relaxation=relaxation,comparisons=rows,
        note='Fixed numerical m,q,Q, 1% numerical mean-extrema bounds, full neighbor prior. Normal rotation preserves prior phase area fraction. Exact field used only for scoring. Level2/angle0 reused from predecessor.')
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
