"""Diagnose centroid aliasing of neighbor-interface priors; no LP rerun."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_moment_batch import setup
import moment_ale_evolve as model
from neighbor_density_audit import interface_prior
from density_sensitivity_audit import rotated_prior
from subcell_density_audit import triangles,triangle_moments
from adaptive_coupled_mass import atomic


def main():
    cases={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
        rebuilt,_=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
        np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
        a=np.asarray(cp['state'])[:-4].reshape(-1,20);cells,edges=model.geometry(a[:,:2],a[:,2],cp['t'])
        centers=np.array([np.asarray(Polygon(v).centroid.coords)[0] for v in cells])
        means=np.array([m/Polygon(v).area for m,v in zip(mass,cells)])
        low=float(means.min());high=float(means.max());neighbors=[set() for _ in mass]
        for i,j,*_ in edges:
            if j>=0:neighbors[i].add(j);neighbors[j].add(i)
        base=[interface_prior(i,v,centers,means,neighbors) for i,v in enumerate(cells)]
        rows=[]
        for level in (1,2,3):
            for angle in (-15.,0.,15.):
                records=[]
                for i,v in enumerate(cells):
                    seed,meta=rotated_prior(v,*base[i],low,high,angle)
                    parts=triangles(v,level);points=np.array([p.mean(0) for p in parts]);areas=np.array([triangle_moments(p)[0] for p in parts])
                    values=seed(points);original=base[i][0](points)
                    records.append(dict(cell=i,centroid_changed=int(np.count_nonzero(values!=original)),subcells=len(parts),
                        sampled_prior_mass_relative_error=float(abs(values@areas-mass[i])/mass[i]),flat=meta['flat']))
                rows.append(dict(level=level,angle_degrees=angle,cells=records))
        cases[name]=rows
    atomic('sampling_report.json',dict(completed=True,cases=cases,note='Analytic prior mass equals stored cell mass by area-fraction construction; this compares centroid-sampled mass before LP correction.'))


if __name__=='__main__':main()
