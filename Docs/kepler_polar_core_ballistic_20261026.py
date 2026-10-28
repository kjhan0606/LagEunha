"""Check central atmospheric Voronoi area noise of two q=2 layouts."""

from pathlib import Path
import json

import numpy as np

from kepler_hybrid_polar_ballistic_20261026 import read_glass, polygon_areas, EPS

SOURCES={'q2_center':Path('/tmp/kepler_hybrid_polar_q2_384.bin'),
         'q2_no_center':Path('/tmp/kepler_hybrid_polar_q2c_384.bin')}


def main():
    cases=[]
    for name,path in SOURCES.items():
        xy0=read_glass(path)
        dx0,dy0=xy0[:,0]-3.,xy0[:,1]-3.
        r0=np.hypot(dx0,dy0)
        selected=np.flatnonzero(r0<.5)
        area0=polygon_areas(xy0,selected)
        mass=(.01+(r0[selected]/.5)**3)*area0
        rows=[]
        for t in (.825,2.221441469079183,120.,600.):
            theta=t*(r0*r0+EPS*EPS)**(-.75)
            xy=np.column_stack((3.+dx0*np.cos(theta)-dy0*np.sin(theta),
                                3.+dx0*np.sin(theta)+dy0*np.cos(theta)))
            area=polygon_areas(xy,selected)
            err=abs(area0/area-1.)
            rows.append(dict(time=t,central_density_change_mass_weighted=float(np.sum(mass*err)/np.sum(mass)),
                             central_density_change_p95=float(np.quantile(err,.95)),
                             central_density_change_max=float(np.max(err)),
                             min_area_ratio=float(np.min(area/area0))))
        cases.append(dict(layout=name,central_cells=len(selected),
                          initial_area_min=float(np.min(area0)),
                          initial_area_max=float(np.max(area0)),results=rows))
    result=dict(cases=cases,limitation='Exact-orbit Voronoi geometry only; numerical forces and RK4 positivity are not tested')
    path=Path(__file__).with_name('KEPLER_POLAR_CORE_BALLISTIC_2026-10-26.json')
    with path.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
