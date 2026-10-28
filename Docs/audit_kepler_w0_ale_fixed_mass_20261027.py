"""Incompressible shear counterexample for W0 Voronoi vs exact ALE mass."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import Voronoi


def polygon_area(vor,index):
    region=vor.regions[vor.point_region[index]]
    if not region or -1 in region:
        raise ValueError("Selected cell is unbounded")
    p=vor.vertices[region]
    return .5*abs(float(np.sum(p[:,0]*np.roll(p[:,1],-1)
                              -p[:,1]*np.roll(p[:,0],-1))))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    rng=np.random.default_rng(17)
    grid=np.array([(x,y) for y in range(-3,4) for x in range(-3,4)],dtype=float)
    sites=grid+rng.uniform(-.12,.12,grid.shape)
    velocity=np.column_stack((np.zeros(len(sites)),sites[:,0]))
    cell=24
    epsilon=1e-5
    vor=Voronoi(sites)
    area=polygon_area(vor,cell)
    derivative_fd=(polygon_area(Voronoi(sites+epsilon*velocity),cell)
                   -polygon_area(Voronoi(sites-epsilon*velocity),cell))/(2*epsilon)
    geometric_rate=0.
    fluid_flux=0.
    closure=np.zeros(2)
    faces=[]
    for pair,verts in zip(vor.ridge_points,vor.ridge_vertices):
        if cell not in pair:
            continue
        if -1 in verts or len(verts)!=2:
            raise ValueError("Selected cell has a non-finite face")
        j=int(pair[1] if pair[0]==cell else pair[0])
        c=.5*(vor.vertices[verts[0]]+vor.vertices[verts[1]])
        length=float(np.linalg.norm(vor.vertices[verts[1]]-vor.vertices[verts[0]]))
        delta=sites[j]-sites[cell]
        d=float(np.linalg.norm(delta))
        normal=delta/d
        w_n=float((np.dot(sites[j],velocity[j])-np.dot(sites[cell],velocity[cell])
                   -np.dot(velocity[j]-velocity[cell],c))/d)
        fluid=np.array([0.,c[0]])
        fluid_n=float(np.dot(fluid,normal))
        geometric_rate+=length*w_n
        fluid_flux+=length*fluid_n
        closure+=length*normal
        faces.append({"neighbor":j,"length":length,"face_normal_speed":w_n,
                      "fluid_normal_speed":fluid_n,"relative_normal_speed":fluid_n-w_n})
    rho=1.
    mass_rate_ale=-rho*(fluid_flux-geometric_rate)
    checks={"finite_closed_cell":len(faces)>=3 and np.linalg.norm(closure)<1e-12,
            "incompressible_fluid_flux":abs(fluid_flux)<1e-12,
            "face_speed_matches_area_derivative":abs(geometric_rate-derivative_fd)<1e-8,
            "ale_mass_rate_matches_rho_area_rate":abs(mass_rate_ale-rho*derivative_fd)<1e-8,
            "nonzero_mass_rate_at_W0":abs(mass_rate_ale)>1e-4}
    checks={name:bool(value) for name,value in checks.items()}
    result={"field":"v=(0,x), rho=1, divergence(v)=0, W=0, exact geometric Voronoi faces",
            "cell":cell,"area":area,"finite_difference_epsilon":epsilon,
            "area_rate_finite_difference":derivative_fd,"area_rate_from_faces":geometric_rate,
            "integrated_fluid_normal_flux":fluid_flux,"ale_cell_mass_rate":mass_rate_ale,
            "relative_mass_rate":mass_rate_ale/(rho*area),
            "face_vector_closure_norm":float(np.linalg.norm(closure)),"faces":faces,"checks":checks,
            "scope":"Kinematic counterexample for exact geometric ALE flux, not a native "
                    "Kepler run or a proof that fixed-mass schemes cannot converge."}
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps({"area":area,"area_rate":derivative_fd,"ale_mass_rate":mass_rate_ale,
                      "failed":[k for k,v in checks.items() if not v]},indent=2))
    if not all(checks.values()):
        raise SystemExit("W0 ALE counterexample checks failed")


if __name__=='__main__':
    main()
