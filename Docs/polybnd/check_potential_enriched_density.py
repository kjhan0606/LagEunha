"""One-cell representability control for an extra potential-density moment.

Keeps geometry and original m,q,Q fixed. This tests only the missing tangent
direction in the density family, not a coupled timestep or gravity solver.
"""
import argparse,json,os
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from audit_density_potential_dual import basis
from audit_kepler_energy_step import phi,atomic
from kepler_disk_reference import polar_rule
from positive_exponential_tilt import fit


def check(data,steps=(.00025,.000125,.0000625)):
    poly=np.asarray(data['polygon']);raw=np.asarray(data['raw'])
    os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER']=str(data['density_quadrature_order'])
    density=base.decode(poly,raw)['density'];coefficient=np.asarray(data['coefficient']);sigma=data['sigma']
    x,w=polar_rule(poly,16,radii=(.25,.5));b=basis(x,density);psi=(phi(x)-b@coefficient)/sigma
    F=np.column_stack([b[:,1:],psi]);weights=w*density.value(x)/raw[0]
    native=basis(density.points,density)
    target_b=density.weighted_mass@native/raw[0]
    rows=[]
    for dt in steps:
        U=data['potential_mean']+dt*data['missing_rate']/raw[0]
        target=np.r_[target_b[1:],(U-coefficient@target_b)/sigma]
        try:theta,logZ,p,diag=fit(weights,F,target)
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            rows.append(dict(dt=dt,converged=False,error=str(exc)));continue
        checks=[]
        for order in (24,32):
            xx,ww=polar_rule(poly,order,radii=(.25,.5));bb=basis(xx,density)
            FF=np.column_stack([bb[:,1:],(phi(xx)-bb@coefficient)/sigma])
            rho=density.value(xx)*np.exp(FF@theta-logZ);dm=ww*rho
            moments=dm@bb/raw[0]
            checks.append(dict(order=order,mass_relative=float(dm.sum()/raw[0]-1),
                retained_moment_error=float(np.max(abs(moments-target_b))),
                potential_mean_error=float(dm@phi(xx)/raw[0]-U),minimum_sampled_density=float(rho.min())))
        rows.append(dict(dt=dt,converged=True,fit=diag,checks=checks,theta=theta))
    return dict(cell=data['cell'],rows=rows,full_validation=False,
        scope='fixed-cell extra potential moment representability; no coupled evolution')


def run(fixture):
    report=check(json.loads(Path(fixture).read_text()))
    atomic('report.json',report);print(json.dumps(report,default=lambda x:x.tolist()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');run(p.parse_args().fixture)
