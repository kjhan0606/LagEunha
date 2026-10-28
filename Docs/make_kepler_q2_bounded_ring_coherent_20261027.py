"""Make a bounded q2 control with one 5% angular phase offset per ring."""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.spatial import cKDTree

from make_kepler_bounded_polygon_ic_20261027 import (
    BOUNDARY_RADIUS, BOUNDARY_SIDES, BOX, CENTER, GUARD_RADII,
    bounded_areas,
)


SOURCE_SHA256="37fcc5f1600bca9083b90aa2a24f98f973acff517266fe623d07e43fb92558c2"
N=55351
ETA=.05
SEED=20261027


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("source",type=Path)
    ap.add_argument("output_prefix",type=Path)
    args=ap.parse_args()
    glass_path=Path(str(args.output_prefix)+".glass.bin")
    area_path=Path(str(args.output_prefix)+".areas.bin")
    manifest_path=Path(str(args.output_prefix)+".json")
    if any(p.exists() for p in (glass_path,area_path,manifest_path)):
        raise FileExistsError("Refusing to overwrite coherent-ring IC")
    raw=args.source.read_bytes()
    if sha(raw)!=SOURCE_SHA256 or len(raw)!=20+16*N:
        raise ValueError("Source is not accepted bounded q2 glass")
    if struct.unpack_from("<idd",raw)!=(N,BOX,BOX):
        raise ValueError("Unexpected source header")
    x=np.frombuffer(raw,dtype="<f8",count=N,offset=20).copy()
    y=np.frombuffer(raw,dtype="<f8",count=N,offset=20+8*N).copy()
    x0,y0=x.copy(),y.copy()
    rx,ry=x-CENTER[0],y-CENTER[1]
    radius=np.hypot(rx,ry)
    polar=(radius>=.48)&(radius<2.02)
    disk=(radius>=.5)&(radius<2.)
    rings,inverse,population=np.unique(np.round(radius[polar],8),
                                       return_inverse=True,return_counts=True)
    if len(rings)!=70 or population.min()<200 or int(disk.sum())!=49038:
        raise ValueError("Unexpected ring structure")
    ring_offset=np.random.default_rng(SEED).uniform(-ETA,ETA,size=len(rings))
    angular_shift=ring_offset[inverse]*2*np.pi/population[inverse]
    theta=np.arctan2(ry[polar],rx[polar])+angular_shift
    x[polar]=CENTER[0]+radius[polar]*np.cos(theta)
    y[polar]=CENTER[1]+radius[polar]*np.sin(theta)
    if not (np.array_equal(x[~polar],x0[~polar]) and
            np.array_equal(y[~polar],y0[~polar])):
        raise ValueError("Nonpolar sites changed")
    radial_error=float(np.max(np.abs(np.hypot(x-CENTER[0],y-CENTER[1])-radius)))
    if radial_error>2e-15 or radius.max()>=BOUNDARY_RADIUS*np.cos(np.pi/BOUNDARY_SIDES):
        raise ValueError("Radius or wall changed")
    xy=np.column_stack((x,y))
    nearest=cKDTree(xy).query(xy,k=2,workers=1)[0][:,1]
    if nearest.min()<=.005:
        raise ValueError("Coherent phase created a near collision")
    area4,_,clipped4=bounded_areas(xy,GUARD_RADII[0])
    area5,_,clipped5=bounded_areas(xy,GUARD_RADII[1])
    area_error=float(np.max(np.abs(area4/area5-1.)))
    polygon_area=.5*BOUNDARY_SIDES*BOUNDARY_RADIUS**2*np.sin(2*np.pi/BOUNDARY_SIDES)
    area_sum=float(np.sum(area5,dtype=np.longdouble))
    if area_error>1e-8 or abs(area_sum-polygon_area)>1e-7:
        raise ValueError("Independent bounded area check failed")
    glass=struct.pack("<idd",N,BOX,BOX)+x.astype("<f8").tobytes()+y.astype("<f8").tobytes()
    sidecar=struct.pack("<idd",N,BOX,BOX)+area5.astype("<f8").tobytes()
    report={
        "source":str(args.source),"source_sha256":sha(raw),
        "glass":str(glass_path),"glass_sha256":sha(glass),
        "areas":str(area_path),"areas_sha256":sha(sidecar),
        "sites":N,"disk_sites":int(disk.sum()),
        "perturbed_sites":int(polar.sum()),"perturbed_rings":len(rings),
        "phase_mode":"one shared angle offset per ring",
        "phase_fraction_of_ring_spacing":ETA,"seed":SEED,
        "max_radius_error":radial_error,"minimum_pair_distance":float(nearest.min()),
        "boundary_cells_per_guard":[clipped4,clipped5],
        "area_sum":area_sum,"polygon_area":polygon_area,
        "max_guard_area_relative_difference":area_error,
        "scope":"Same bounded q2 exterior/pressure/orbits as aligned and site-jitter controls; W0 one-time fixed masses use the perturbed power areas."
    }
    with glass_path.open("xb") as stream: stream.write(glass)
    with area_path.open("xb") as stream: stream.write(sidecar)
    with manifest_path.open("x") as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
