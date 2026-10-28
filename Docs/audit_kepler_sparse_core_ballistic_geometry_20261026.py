"""Exact-orbit Voronoi geometry of the innermost Kepler cells.

Only the local r<0.5 sites are tessellated; all selected core cell vertices
must stay inside r<0.18, so exterior sites cannot set their bisectors.
This is a geometry-only control and does not evolve pressure or gravity.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi, cKDTree

from audit_kepler_sparse_core_force_20261026 import glass, area, CENTER, EPS


def orbit(local, t):
    offset = local-CENTER
    r2 = np.sum(offset*offset,axis=1)
    phase = t*(r2+EPS*EPS)**(-.75)
    c,s = np.cos(phase),np.sin(phase)
    return CENTER + np.column_stack((offset[:,0]*c-offset[:,1]*s,
                                     offset[:,0]*s+offset[:,1]*c))


def core_areas(local, core_local, t):
    positions = orbit(local,t)
    vor = Voronoi(positions, qhull_options="Qbb Qc Qz")
    areas, vertex_radii = [], []
    for i in core_local:
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            raise ValueError(f"Unbounded core cell at t={t}, i={i}")
        vertices = vor.vertices[region]
        areas.append(area(vertices))
        vertex_radii.append(np.max(np.linalg.norm(vertices-CENTER,axis=1)))
    a = np.asarray(areas)
    if not np.all(np.isfinite(a)&(a>0)) or max(vertex_radii)>=.18:
        raise ValueError(f"Invalid local core geometry at t={t}")
    return a,float(max(vertex_radii))


def evaluate(name, path, positions, tolerance, times):
    xy = glass(path)
    distance, core = cKDTree(xy).query(positions)
    if np.max(distance)>tolerance or len(set(core.tolist()))!=len(positions):
        raise ValueError(f"Missing angular counterpart core sites in {path}")
    initial_radius = np.linalg.norm(xy-CENTER,axis=1)
    selected = np.flatnonzero(initial_radius<.5)
    local = xy[selected]
    lookup = {int(g):int(i) for i,g in enumerate(selected)}
    core_local = np.asarray([lookup[int(i)] for i in core],dtype=int)
    a0, max_vertex0 = core_areas(local,core_local,0.)
    rho = .01+(initial_radius[core]/.5)**3
    mass = rho*a0
    rows=[]
    for t in times:
        at,max_vertex = core_areas(local,core_local,float(t))
        change = a0/at-1.
        rows.append({
            "time":float(t),
            "density_change_signed_mass_weighted":float(np.sum(mass*change)/np.sum(mass)),
            "density_change_abs_mass_weighted":float(np.sum(mass*abs(change))/np.sum(mass)),
            "density_change_abs_max":float(np.max(abs(change))),
            "area_min":float(np.min(at)),
            "maximum_vertex_radius":max_vertex,
        })
    return {
        "name":name,"glass":str(path),"sites":len(xy),
        "local_sites_within_r0p5":len(local),
        "core_sites":len(core),
        "core_site_radius_minmax":[float(np.min(initial_radius[core])),
                                   float(np.max(initial_radius[core]))],
        "max_distance_to_profile_sampled_ring":float(np.max(distance)),
        "initial_area_minmax":[float(np.min(a0)),float(np.max(a0))],
        "initial_max_vertex_radius":max_vertex0,
        "series":rows,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("reference",type=Path)
    parser.add_argument("narrow",type=Path)
    parser.add_argument("wide",type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    narrow=glass(args.narrow)
    rr=np.linalg.norm(narrow-CENTER,axis=1)
    core_positions=narrow[abs(rr-np.min(rr))<1e-9]
    if len(core_positions)!=8:
        raise ValueError("Expected eight profile-sampled innermost sites")
    times=np.concatenate((np.linspace(0.,1.2,61),[2.221441469079183,6.5,10.,120.,600.]))
    result={
        "softening":EPS,
        "times":times.tolist(),
        "cases":[
            evaluate("full_q2_nearby",args.reference,core_positions,.02,times),
            evaluate("narrow",args.narrow,core_positions,1e-12,times),
            evaluate("wide",args.wide,core_positions,1e-8,times),
        ],
        "scope":"Exact softened circular-orbit local Voronoi geometry of innermost sites. Density uses one-time site-density times initial area as fixed mass. No hydrodynamic forces, heating, or native runtime are represented.",
    }
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2)
        stream.write("\n")
    print(json.dumps({"cases":[{
        "name":c["name"],
        "core_site_radius":c["core_site_radius_minmax"],
        "worst_density_abs":max(row["density_change_abs_max"] for row in c["series"]),
        "worst_mass_weighted_density_abs":max(row["density_change_abs_mass_weighted"] for row in c["series"]),
        "at_t600":c["series"][-1]["density_change_abs_mass_weighted"],
    } for c in result["cases"]]},indent=2))


if __name__=="__main__":
    main()
