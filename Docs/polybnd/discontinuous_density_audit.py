"""Mass-constrained two-density planar reconstruction, diagnostic only.
Uses the benchmark's initial phase densities, not its evolved interface.
No moments or geometry are overwritten in the simulation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from shear_interface_audit import moments
from adaptive_coupled_mass import atomic


def clip(vertices,normal,offset):
    out=[]
    for a,b in zip(vertices,np.roll(vertices,-1,axis=0)):
        fa=float(a@normal-offset);fb=float(b@normal-offset)
        if fa>=0:out.append(a)
        if (fa>=0)!=(fb>=0):out.append(a+(b-a)*fa/(fa-fb))
    return Polygon(out) if len(out)>=3 else Polygon()


def fit(poly,target,low,high):
    vertices=np.asarray(poly.exterior.coords[:-1]);area=poly.area;m=target[0]
    fraction=(m/area-low)/(high-low)
    if fraction < -1e-12 or fraction > 1+1e-12:
        return dict(feasible=False,fraction=float(fraction)),None
    fraction=float(np.clip(fraction,0,1));total=moments(poly)
    center=total[1:3]/area;scale=np.sqrt(area)
    def features(raw):
        q=raw[1:3];Q=raw[3:7].reshape(2,2)
        S=Q-np.outer(center,q)-np.outer(q,center)+raw[0]*np.outer(center,center)
        return np.r_[(q-raw[0]*center)/(m*scale),S.ravel()/(m*area)]
    goal=features(target)
    def candidate(theta):
        normal=np.array([np.cos(theta),np.sin(theta)])
        lo,hi=float((vertices@normal).min()),float((vertices@normal).max())
        if fraction==0:return Polygon()
        if fraction==1:return poly
        for _ in range(48):
            mid=(lo+hi)/2;cut=clip(vertices,normal,mid)
            if cut.area>fraction*area:lo=mid
            else:hi=mid
        return clip(vertices,normal,(lo+hi)/2)
    def objective(theta):
        raw=low*total+(high-low)*moments(candidate(theta))
        return float(np.sum((features(raw)-goal)**2))
    angles=np.linspace(0,2*np.pi,24,endpoint=False)
    scores=[objective(a) for a in angles];theta=angles[int(np.argmin(scores))]
    optimized=minimize_scalar(objective,bounds=(theta-np.pi/12,theta+np.pi/12),method='bounded',options={'xatol':1e-10})
    if optimized.success and optimized.fun<min(scores):theta=optimized.x
    cut=candidate(theta);raw=low*total+(high-low)*moments(cut)
    error=float(np.max(abs(features(raw)-goal)))
    return dict(feasible=True,fraction=fraction,mass_relative_error=float(abs(raw[0]-m)/m),
        max_scaled_moment_error=error,moment_gate=error<1e-8,angle=float(theta)),cut


def audit(name):
    c=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=c['config']
    if cfg['kind']!='shear_contact' or c['t']!=.01:raise RuntimeError('configuration mismatch')
    mass,y0=setup('shear_contact',cfg['n'],cfg['amplitude']);a0=y0[:-4].reshape(-1,20)
    np.testing.assert_allclose(c['mass'],mass,rtol=0,atol=1e-15)
    initial,_=model.geometry(a0[:,:2],a0[:,2],0.)
    low=1.;high=cfg['amplitude'];rho=np.where(a0[:,0]<.5,low,high)
    refs=[]
    for p in initial:
        v=p.copy();v[:,0]+=c['t']*(.2+1.1*v[:,1]);refs.append(Polygon(v))
    a=np.asarray(c['state'])[:-4].reshape(-1,20);cells,_=model.geometry(a[:,:2],a[:,2],c['t'])
    rows=[];l1=0.;valid=True
    for i,vertices in enumerate(cells):
        poly=Polygon(vertices);target=np.r_[mass[i],a[i,3:9]]
        row,cut=fit(poly,target,low,high);row['cell']=i;rows.append(row)
        if not row['feasible']:valid=False;continue
        for ref,density in zip(refs,rho):
            overlap=poly.intersection(ref);ha=overlap.intersection(cut).area
            l1+=abs(high-density)*ha+abs(low-density)*max(0.,overlap.area-ha)
    result=dict(case=name,status='completed',all_cells_feasible=valid,cells=rows,
        density_L1_relative=l1/float(mass.sum()) if valid else None,
        all_moments_preserved=valid and all(r.get('moment_gate',False) for r in rows),
        note='Known initial phase densities are supplied. Interface direction uses numerical moments only. Single-line local search is not a global feasibility proof; candidate is not used in evolution.')
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
