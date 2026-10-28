"""Midpoint mass-transport and off-center torque proxy for native weights."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weighted_face_slip_t65_20261027 import finite_faces
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
DT = .001
EPS = .01175


def cross(a, b):
    return a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]


def case_metrics(xy, velocity, r0, active, lookup, weight, rate,
                 density, mass, disk_angular_momentum):
    _, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"]:
        raise ValueError("Hidden disk power cell")
    i,j,c,length = finite_faces(triangles, centers, lookup>=0)
    dr = xy[j]-xy[i]
    d = np.linalg.norm(dr, axis=1)
    n = dr/d[:, None]
    u_fixed = (np.einsum("ij,ij->i", xy[j], velocity[j])
               - np.einsum("ij,ij->i", xy[i], velocity[i])
               - np.einsum("ij,ij->i", velocity[j]-velocity[i], c))/d
    u_face = u_fixed+(rate[i]-rate[j])/(2*d)
    rface2 = np.einsum("ij,ij->i", c, c)
    fluid = np.column_stack((-c[:, 1], c[:, 0]))
    fluid *= (rface2+EPS*EPS)[:, None]**(-.75)
    slip = np.einsum("ij,ij->i", fluid, n)-u_face
    edge_rate = .5*(density[i]+density[j])*length*slip
    net = np.zeros(N)
    np.add.at(net, i, edge_rate)
    np.add.at(net, j, -edge_rate)
    gross = float(np.sum(np.abs(edge_rate), dtype=np.longdouble))
    disk_mass = float(np.sum(mass[active], dtype=np.longdouble))
    if abs(float(np.sum(net, dtype=np.longdouble))) > 1e-12*max(gross, 1e-30):
        raise ValueError("Face mass transport failed pairwise cancellation")
    fraction = DT*np.abs(net[active])/mass[active]
    momentum = edge_rate[:, None]*fluid
    orbital_torque = cross(xy[j]-xy[i], momentum)
    spin_i = -cross(c-xy[i], momentum)
    spin_j = cross(c-xy[j], momentum)
    torque_residual = orbital_torque+spin_i+spin_j
    return {
        "disk_disk_faces": len(i),
        "length_weighted_rms_face_slip": float(np.sqrt(np.average(slip*slip, weights=length))),
        "gross_edge_mass_rate_over_disk_mass": gross/disk_mass,
        "one_step_gross_edge_mass_over_disk_mass": DT*gross/disk_mass,
        "one_step_net_cell_mass_fraction_p50": float(np.quantile(fraction, .5)),
        "one_step_net_cell_mass_fraction_p99": float(np.quantile(fraction, .99)),
        "one_step_net_cell_mass_fraction_max": float(np.max(fraction)),
        "one_step_absolute_face_orbital_torque_over_disk_L":
            DT*float(np.sum(np.abs(orbital_torque), dtype=np.longdouble))/abs(disk_angular_momentum),
        "max_pair_orbital_plus_spin_torque_residual": float(np.max(np.abs(torque_residual))),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=("aligned", "jitter"))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    base = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
    run = base/"gpu_energy_neutral_weighted_t6501_guard_a10_v2"/args.case
    receipt = json.loads((run/"audit_t6501.json").read_text())
    if not all(receipt["checks"].values()):
        raise ValueError("Unaccepted native weighted step")
    start_path,end_path = Path(receipt["restart"]),Path(receipt["endpoint"])
    if (digest(start_path)!=receipt["restart_sha256"]
            or digest(end_path)!=receipt["endpoint_sha256"]):
        raise ValueError("Changed native snapshots")
    docs = Path(__file__).resolve().parent
    initial_name = ("KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json" if args.case=="aligned"
                    else "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json")
    initial_receipt = json.loads((docs/initial_name).read_text())
    if (initial_receipt["endpoint_sha256"]!=receipt["restart_sha256"]
            or digest(Path(initial_receipt["initial"]))!=initial_receipt["initial_sha256"]):
        raise ValueError("Changed accepted lineage")
    initial, _ = ordered(Path(initial_receipt["initial"]), N)
    start, _ = ordered(start_path, N)
    end, _ = ordered(end_path, N)
    xy0 = np.column_stack((initial["x"].astype(float)-3., initial["y"].astype(float)-3.))
    xy_start = np.column_stack((start["x"].astype(float)-3.,start["y"].astype(float)-3.))
    xy_end = np.column_stack((end["x"].astype(float)-3.,end["y"].astype(float)-3.))
    xy = .5*(xy_start+xy_end)
    velocity = .5*(np.column_stack((start["vx"],start["vy"]))
                   +np.column_stack((end["vx"],end["vy"])))
    density = .5*(start["den"].astype(float)+end["den"].astype(float))
    mass = start["mass"].astype(float)
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0>=.5)&(r0<2.))
    lookup = np.full(N,-1,dtype=int)
    lookup[active] = np.arange(len(active))
    weight0 = start["w2"].astype(float)
    weight1 = end["w2"].astype(float)
    weight = .5*(weight0+weight1)
    rate = (weight1-weight0)/DT
    disk_l = float(np.sum(mass[active]*cross(xy[active],velocity[active]),
                          dtype=np.longdouble))
    zero = np.zeros(N)
    w0 = case_metrics(xy,velocity,r0,active,lookup,zero,zero,density,mass,disk_l)
    weighted = case_metrics(xy,velocity,r0,active,lookup,weight,rate,density,mass,disk_l)
    result = {
        "case": args.case, "start_sha256": receipt["restart_sha256"],
        "end_sha256": receipt["endpoint_sha256"], "dt": DT,
        "W0_on_same_midpoint_sites": w0,
        "applied_weight_rate_midpoint": weighted,
        "gross_mass_rate_ratio_weighted_over_W0":
            weighted["gross_edge_mass_rate_over_disk_mass"]/
            w0["gross_edge_mass_rate_over_disk_mass"],
        "net_cell_p99_ratio_weighted_over_W0":
            weighted["one_step_net_cell_mass_fraction_p99"]/
            w0["one_step_net_cell_mass_fraction_p99"],
        "scope": "Analytic face-centroid Kepler velocity, arithmetic face density, "
                 "disk-disk faces only. Midpoint endpoint-difference rate, not native "
                 "Riemann flux or exact RK4 stage audit. Spin terms are a pairwise "
                 "bookkeeping proxy, not implemented cell rotation.",
    }
    with args.output.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    print(json.dumps({"case": args.case,
                      "gross_rate_ratio": result["gross_mass_rate_ratio_weighted_over_W0"],
                      "net_p99_ratio": result["net_cell_p99_ratio_weighted_over_W0"],
                      "W0_p99": w0["one_step_net_cell_mass_fraction_p99"],
                      "weighted_p99": weighted["one_step_net_cell_mass_fraction_p99"]}, indent=2))


if __name__ == "__main__":
    main()
