"""Export accepted W0 FaceCSR and independent CPU projected-rate reference."""

import argparse
import json
from pathlib import Path
import struct

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from export_kepler_q2_weight_rate_csr_t65_20261027 import (
    finite_power_edges, write_array,
)
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
MAGIC = 0x4B574550


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accepted_receipt", type=Path)
    parser.add_argument("reference_rates", type=Path)
    parser.add_argument("case", choices=("aligned", "jitter"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if receipt["time"] != 6.5 or not all(receipt["checks"].values()):
        raise ValueError("Expected accepted t=6.5 endpoint")
    if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
        raise ValueError("Endpoint hash changed")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3.,p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3.,p["y"].astype(float)-3.))
    r0 = np.linalg.norm(xy0, axis=1)
    active = ((r0 >= .5) & (r0 < 2.)).astype(np.uint8)
    lookup = np.full(N,-1,dtype=int)
    lookup[active.astype(bool)] = np.arange(active.sum())
    with np.load(args.reference_rates) as saved:
        if not np.array_equal(saved["active"],active):
            raise ValueError("Active set changed")
        reference = saved[args.case].copy()
    area, meta, triangles, centers = power_areas(
        xy,np.zeros(N),np.flatnonzero(active),lookup,return_geometry=True)
    if meta["hidden_active"] or np.any(area <= 0):
        raise ValueError("Invalid W0 geometry")
    i,j,c1,c2 = finite_power_edges(triangles,centers,active.astype(bool))
    if len(i) != 148557:
        raise ValueError(f"Unexpected unique face count {len(i)}")
    length = np.linalg.norm(c2-c1,axis=1)
    d = np.linalg.norm(xy[j]-xy[i],axis=1)
    g = length/(2*d)
    pressure = p["pressure"].astype("<f4")
    gradient = np.zeros(N)
    pair = g*(pressure[i].astype(float)-pressure[j].astype(float))
    np.add.at(gradient,i,pair)
    np.add.at(gradient,j,-pair)
    gradient[~active.astype(bool)] = 0.
    numerator = np.dot(gradient,reference)
    denominator = np.dot(gradient,gradient)
    if not denominator > 0:
        raise ValueError("Zero pressure gradient")
    alpha = numerator/denominator
    projected = reference-alpha*gradient
    projected[~active.astype(bool)] = 0.
    owner = np.r_[i,j]
    neighbor = np.r_[j,i]
    corner1 = np.vstack((c1-xy[i],c2-xy[j]))
    corner2 = np.vstack((c2-xy[i],c1-xy[j]))
    order = np.lexsort((neighbor,owner))
    owner,neighbor = owner[order],neighbor[order]
    corner1,corner2 = corner1[order],corner2[order]
    offset = np.r_[0,np.cumsum(np.bincount(owner,minlength=N))]
    with args.output.open("xb") as stream:
        stream.write(struct.pack("<4i",MAGIC,N,len(owner),16))
        write_array(stream,offset,"<i4")
        write_array(stream,neighbor,"<i4")
        write_array(stream,np.zeros(len(owner)),"<i4")
        for array in (corner1[:,0],corner1[:,1],corner2[:,0],corner2[:,1],
                      p["x"],p["y"],p["vx"],p["vy"]):
            write_array(stream,array,"<f8")
        write_array(stream,pressure,"<f4")
        write_array(stream,active,"u1")
        for array in (reference,gradient,projected):
            write_array(stream,array,"<f8")
    result = {
        "case":args.case,"binary":str(args.output.resolve()),
        "sha256":digest(args.output),
        "accepted_endpoint_sha256":receipt["endpoint_sha256"],
        "reference_rates_sha256":digest(args.reference_rates),
        "sites":N,"active":int(active.sum()),
        "undirected_faces_touching_active":len(i),
        "directed_faces":len(owner),
        "pressure_work_before":float(numerator),
        "pressure_work_after":float(np.dot(gradient,projected)),
        "projection_alpha":float(alpha),
        "scope":"CPU-exported actual W0 t=6.5 power geometry and pressure-gradient reference for a separate CUDA projection probe; no hydro state change.",
    }
    with args.manifest.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
