"""Small-mesh derivative audit of bounded Voronoi cells across a site/wall crossing.

This tests the geometric energy extension, not an accepted outside-domain
fluid state or native GPU implementation. All cells must remain nonempty.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import audit_kepler_bounded_varpress_adjoint_20261027 as base


def analyze():
    sites = np.array([(x+.5,y+.5) for y in range(4) for x in range(4)])
    sites += np.random.default_rng(17).uniform(-.12,.12,size=sites.shape)
    theta = 4.5*2*np.pi/base.DOMAIN_SIDES
    normal = np.array([np.cos(theta),np.sin(theta)])
    apothem = base.RADIUS*np.cos(np.pi/base.DOMAIN_SIDES)
    sites[-1] = base.CENTER+(apothem-1.e-4)*normal
    entropy = np.linspace(.8,1.2,len(sites))*base.tessellate(sites)[1]**base.GAMMA
    velocity = np.column_stack((np.linspace(.1,.2,16),np.linspace(-.1,.03,16)))
    cases=[]
    for penetration in [-1.e-4,0.,1.e-4]:
        trial = sites.copy()
        trial[-1] = base.CENTER+(apothem+penetration)*normal
        polys,areas = base.tessellate(trial)
        pressure = entropy/areas**base.GAMMA
        force,rate,interior,wall,traction,closure,nwall,reciprocal=base.face_operator(trial,velocity,pressure,polys)
        errors=[]
        for h in [1.e-4,5.e-5,2.5e-5]:
            numerical=[]
            for axis in range(2):
                plus=trial.copy();minus=trial.copy()
                plus[-1,axis]+=h;minus[-1,axis]-=h
                numerical.append(-(base.energy(base.tessellate(plus)[1],entropy)-
                                   base.energy(base.tessellate(minus)[1],entropy))/(2*h))
            errors.append({'h':h,'gradient_error':float(np.linalg.norm(force[-1]-numerical))})
        power=float(np.sum(force*velocity)-np.dot(pressure,interior+wall))
        cases.append({'site_penetration':penetration,'minimum_cell_area':float(areas.min()),
                      'outside_site_cell_area':float(areas[-1]),'force_on_crossing_site':force[-1].tolist(),
                      'gradient_errors':errors,'adjoint_power_residual':power,
                      'partition_error':float(abs(sum(areas)-base.area(base.boundary()))),
                      'wall_traction_residual':float(np.linalg.norm(force.sum(axis=0)+traction)),
                      'wall_faces':nwall,'closure_residual':closure,'reciprocal_residual':reciprocal})
    checks={'all_cells_positive':all(c['minimum_cell_area']>0 for c in cases),
            'gradient_matches_inside_on_and_outside':all(c['gradient_errors'][-1]['gradient_error']<2.e-7 for c in cases),
            'static_wall_power_closes':all(abs(c['adjoint_power_residual'])<1.e-12 for c in cases),
            'area_partition_closes':all(c['partition_error']<1.e-12 for c in cases),
            'wall_traction_closes':all(c['wall_traction_residual']<1.e-12 for c in cases)}
    return {'cases':cases,'checks':checks,'passed':all(checks.values()),
            'scope':'Independent 16-cell fixed-polygon variational geometry with one generator just inside/on/outside the boundary. Outside site does not imply empty clipped cell. Supports a local mathematical RHS extension only while every cell stays positive; does not certify native GPU behavior, arbitrary excursions, event accuracy or long-time stability.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();r=analyze()
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r,indent=2))
    if not r['passed']:raise SystemExit(1)
