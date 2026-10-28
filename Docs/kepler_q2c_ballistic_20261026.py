"""Exact-orbit volume control for the centerless q=2 polar disk variant."""

from pathlib import Path
import json

import numpy as np

from kepler_hybrid_polar_ballistic_20261026 import read_glass, polygon_areas, EPS, TIMES

SOURCE = Path('/tmp/kepler_hybrid_polar_q2c_384.bin')


def main():
    xy0 = read_glass(SOURCE)
    dx0,dy0=xy0[:,0]-3.,xy0[:,1]-3.
    radius=np.hypot(dx0,dy0)
    selected=np.flatnonzero((radius>=.5)&(radius<2.))
    a0=polygon_areas(xy0,selected)
    mass=1.01*a0
    rows=[]
    for t in TIMES:
        angle=t*(radius*radius+EPS*EPS)**(-.75)
        xy=np.column_stack((3.+dx0*np.cos(angle)-dy0*np.sin(angle),
                            3.+dx0*np.sin(angle)+dy0*np.cos(angle)))
        at=polygon_areas(xy,selected)
        error=abs(a0/at-1.)
        bins=[]
        for lo,hi in ((.5,.75),(.75,1.),(1.,2.),(.5,2.)):
            take=(radius[selected]>=lo)&(radius[selected]<hi)
            bins.append(dict(rlo=lo,rhi=hi,cells=int(np.sum(take)),
                             density_change_mass_weighted=float(np.sum(mass[take]*error[take])/np.sum(mass[take])),
                             density_change_p95=float(np.quantile(error[take],.95)),
                             density_change_max=float(np.max(error[take]))))
        rows.append(dict(time=t,bins=bins))
    result=dict(layout='polar_q2_no_center',source=str(SOURCE),cells=len(xy0),
                disk_cells=len(selected),initial_disk_area_min=float(np.min(a0)),
                initial_disk_area_max=float(np.max(a0)),results=rows,
                limitation='Exact-orbit unweighted Voronoi geometry only; no hydro or native timestep measurement')
    output=Path(__file__).with_name('KEPLER_Q2C_BALLISTIC_2026-10-26.json')
    with output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
