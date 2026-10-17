"""Read-only exact-field moment audit on saved numerical cells.

Shear maps the initial piecewise constant density polygons with determinant 1.
Intersect those polygons with actual final cells before integrating m,q,Q.
Cell-average L1 is a lower bound on density L1, not an interface norm.
"""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from periodic2d_hllc import quadrature
from adaptive_coupled_mass import atomic


def moments(poly):
    if poly.is_empty or poly.area<1e-25:return np.zeros(7)
    if poly.geom_type!='Polygon':return sum((moments(p) for p in poly.geoms),np.zeros(7))
    if len(poly.interiors):raise RuntimeError('unexpected nonconvex overlap with holes')
    x,w=quadrature(np.asarray(poly.exterior.coords[:-1]))
    return np.r_[w.sum(),w@x,(x.T@(w[:,None]*x)).ravel()]


def exact_on_cells(cells,initial,rho,t):
    material=[]
    for p in initial:
        x=np.asarray(p).copy();x[:,0]+=t*(.2+1.1*x[:,1]);material.append(Polygon(x))
    values=[];coverage=[];shape=[]
    for i,poly in enumerate(cells):
        cell=Polygon(poly);row=np.zeros(7);area=0.
        for ref,density in zip(material,rho):
            overlap=cell.intersection(ref)
            area+=overlap.area;row+=density*moments(overlap)
        values.append(row);coverage.append(abs(area-cell.area))
        shape.append(cell.symmetric_difference(material[i]).area/material[i].area)
    return np.asarray(values),float(max(coverage)),float(max(shape))


def audit(case):
    source=Path('inputs')/case
    cp=json.loads((source/'checkpoint.json').read_text());cfg=cp['config']
    n=cfg['n'];ratio=cfg['amplitude'];t=cp['t']
    if cfg['kind']!='shear_contact' or t!=.01:raise RuntimeError('unexpected saved case')
    mass,y0=setup('shear_contact',n,ratio);a0=y0[:-4].reshape(-1,20)
    np.testing.assert_allclose(cp['mass'],mass,rtol=0,atol=1e-15)
    initial,_=model.geometry(a0[:,:2],a0[:,2],0.)
    rho=np.where(a0[:,0]<.5,1.,ratio)
    initial_exact,coverage0,_=exact_on_cells(initial,initial,rho,0.)
    raw0=np.column_stack([mass,a0[:,3:9]])
    if np.max(abs(initial_exact-raw0))>1e-11:raise RuntimeError('initial reference mismatch')
    a=np.asarray(cp['state'])[:-4].reshape(-1,20)
    cells,_=model.geometry(a[:,:2],a[:,2],t)
    exact,coverage,shape=exact_on_cells(cells,initial,rho,t)
    actual=np.column_stack([mass,a[:,3:9]]);error=actual-exact
    total=float(mass.sum());dm=error[:,0]
    result=dict(case=case,n=n,ratio=ratio,t=t,
        initial_reference_max_error=float(np.max(abs(initial_exact-raw0))),
        max_partition_area_error=max(coverage,coverage0),
        exact_global_mass_relative_error=float(abs(exact[:,0].sum()-total)/total),
        cell_average_density_L1_relative=float(np.sum(abs(dm))/total),
        max_cell_mass_relative_error=float(np.max(abs(dm)/exact[:,0])),
        q_L1_over_total_mass=float(np.sum(abs(error[:,1:3]))/total),
        Q_L1_over_total_mass=float(np.sum(abs(error[:,3:7]))/total),
        max_cell_material_symmetric_difference_fraction=shape,
        cell_mass_errors=dm.tolist(),
        reference_valid=bool(max(coverage,coverage0)<1e-11 and abs(exact[:,0].sum()-total)/total<1e-11),
        note='Exact sheared initial polygon field integrated on numerical cells. Density L1 is a cell-average lower bound. No spatial convergence claim: initial interfaces differ across n.')
    atomic(case+'.json',result)
    return result


def main():
    previous=json.loads(Path('predecessor_completion.json').read_text())
    if not previous.get('passed'):raise RuntimeError('predecessor failed')
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        try:cases[name]=audit(name)
        except Exception as exc:cases[name]=dict(reference_valid=False,error=repr(exc))
    report=dict(status='completed',reference_valid=all(c['reference_valid'] for c in cases.values()),
        cases=cases,note='Diagnostic report only; reference validity does not certify density accuracy.')
    atomic('completion_report.json',report);print(json.dumps(report),flush=True)
    raise SystemExit(0 if report['reference_valid'] else 2)


if __name__=='__main__':main()
