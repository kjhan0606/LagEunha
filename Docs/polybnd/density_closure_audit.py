"""Diagnostic density reconstruction using numerical versus exact moments.
Exact-moment fits are diagnostic only, never fed back into evolution.
"""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from shear_interface_audit import exact_on_cells
from adaptive_coupled_mass import atomic


def rule(poly,order):
    vertices=np.asarray(poly.exterior.coords[:-1]);center=vertices.mean(0)
    gx,gw=np.polynomial.legendre.leggauss(order);gx=(gx+1)/2;gw=gw/2
    points=[];weights=[]
    for a,b in zip(vertices,np.roll(vertices,-1,axis=0)):
        u=a-center;v=b-center;det=abs(u[0]*v[1]-u[1]*v[0])
        for s,w in zip(gx,gw):
            points.extend(center+s*u+((1-s)*gx)[:,None]*v)
            weights.extend(det*(1-s)*w*gw)
    return np.asarray(points),np.asarray(weights)


def audit(name):
    c=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=c['config']
    n=cfg['n'];ratio=cfg['amplitude'];t=c['t']
    if cfg['kind']!='shear_contact' or t!=.01:raise RuntimeError('configuration mismatch')
    mass,y0=setup('shear_contact',n,ratio);a0=y0[:-4].reshape(-1,20)
    np.testing.assert_allclose(c['mass'],mass,rtol=0,atol=1e-15)
    initial,_=model.geometry(a0[:,:2],a0[:,2],0.)
    rho=np.where(a0[:,0]<.5,1.,ratio)
    material=[]
    for p in initial:
        q=p.copy();q[:,0]+=t*(.2+1.1*q[:,1]);material.append(Polygon(q))
    a=np.asarray(c['state'])[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],t)
    exact,coverage,_=exact_on_cells(cells,initial,rho,t)
    if coverage>1e-11:raise RuntimeError('partition mismatch')
    fits=[];regions=[];mixed=0
    for i,poly in enumerate(cells):
        actual=model.Density(poly,mass[i],a[i,3:5],a[i,5:9].reshape(2,2))
        ideal=model.Density(poly,exact[i,0],exact[i,1:3],exact[i,3:7].reshape(2,2))
        fits.append((actual,ideal));pieces=[];labels=set()
        for ref,density in zip(material,rho):
            overlap=Polygon(poly).intersection(ref)
            if overlap.area>1e-25:
                if overlap.geom_type!='Polygon':raise RuntimeError('unexpected overlap type')
                pieces.append((overlap,density));labels.add(float(density))
        mixed+=len(labels)>1;regions.append(pieces)
    rows=[]
    for order in (8,16,32):
        errors=np.zeros(3);actual_mass=np.zeros(len(mass));ideal_mass=np.zeros(len(mass))
        for i,((actual,ideal),pieces) in enumerate(zip(fits,regions)):
            for poly,density in pieces:
                x,w=rule(poly,order);f=actual.value(x);g=ideal.value(x)
                errors+=np.array([w@abs(f-density),w@abs(g-density),w@abs(f-g)])
                actual_mass[i]+=w@f;ideal_mass[i]+=w@g
        rows.append(dict(order=order,actual_density_L1_relative=errors[0]/mass.sum(),
            exact_moment_fit_L1_relative=errors[1]/mass.sum(),
            actual_vs_exact_moment_fit_L1_relative=errors[2]/mass.sum(),
            actual_fit_max_mass_quadrature_error=float(np.max(abs(actual_mass-mass)/mass)),
            exact_fit_max_mass_quadrature_error=float(np.max(abs(ideal_mass-exact[:,0])/exact[:,0]))))
    keys=('actual_density_L1_relative','exact_moment_fit_L1_relative','actual_vs_exact_moment_fit_L1_relative')
    changes={k:abs(rows[-1][k]-rows[-2][k])/max(rows[-1][k],1e-30) for k in keys}
    result=dict(case=name,status='completed',mixed_cells=mixed,rows=rows,
        last_quadrature_relative_changes=changes,
        quadrature_stable=all(v<.01 for v in changes.values()),
        note='Exact-moment fitting is a reconstruction diagnostic, not a best-L1 approximation or an additive error decomposition. Actual fit uses unchanged fixed numerical masses.')
    atomic(name+'.json',result);return result


def main():
    if not json.loads(Path('predecessor_completion.json').read_text())['completed']:
        raise RuntimeError('temporal sensitivity incomplete')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(status='failed',error=repr(exc))
    result=dict(completed=all(c['status']=='completed' for c in cases.values()),cases=cases)
    atomic('completion_report.json',result);print(json.dumps(result),flush=True)
    raise SystemExit(0 if result['completed'] else 2)


if __name__=='__main__':main()
