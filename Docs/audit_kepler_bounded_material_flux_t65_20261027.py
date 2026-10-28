"""Frozen-state proxy for conservative mass transport through Laguerre faces."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weighted_face_slip_t65_20261027 import finite_faces
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
EPS = 0.01175
DT = 0.00193983


def cross(a, b):
    return a[:,0]*b[:,1]-a[:,1]*b[:,0]


def measure(i, j, c, xy, face_velocity, length, density, mass,
            disk_angular_momentum, slip, active, label):
    face_density = .5*(density[i]+density[j])
    edge_rate = face_density*length*slip
    net = np.zeros(N)
    np.add.at(net, i, edge_rate)
    np.add.at(net, j, -edge_rate)
    per_step = DT*np.abs(net[active])/mass[active]
    gross_rate = np.sum(np.abs(edge_rate), dtype=np.longdouble)
    disk_mass = np.sum(mass[active], dtype=np.longdouble)
    if abs(float(np.sum(net, dtype=np.longdouble))) > 1e-12*float(gross_rate):
        raise ValueError(f"{label}: interior-face mass flux is not pair-conservative")
    momentum_rate = edge_rate[:,None]*face_velocity
    orbital_torque = cross(xy[j]-xy[i], momentum_rate)
    spin_i = -cross(c-xy[i], momentum_rate)
    spin_j = cross(c-xy[j], momentum_rate)
    angular_residual = orbital_torque+spin_i+spin_j
    if np.max(np.abs(angular_residual)) > 1e-12*max(np.max(np.abs(orbital_torque)),1e-30):
        raise ValueError(f"{label}: off-center spin torque does not cancel orbital torque")
    absolute_orbital_torque = np.sum(np.abs(orbital_torque),dtype=np.longdouble)
    signed_orbital_torque = np.sum(orbital_torque,dtype=np.longdouble)
    return {"disk_disk_faces": len(i),
            "absolute_edge_mass_transfer_per_unit_time": float(gross_rate),
            "absolute_edge_mass_transfer_rate_over_disk_mass": float(gross_rate/disk_mass),
            "one_step_gross_mass_transfer_over_disk_mass": float(DT*gross_rate/disk_mass),
            "one_step_absolute_net_cell_mass_fraction_median": float(np.median(per_step)),
            "one_step_absolute_net_cell_mass_fraction_p95": float(np.quantile(per_step,.95)),
            "one_step_absolute_net_cell_mass_fraction_p99": float(np.quantile(per_step,.99)),
            "one_step_absolute_net_cell_mass_fraction_max": float(np.max(per_step)),
            "net_all_disk_disk_mass_rate": float(np.sum(net, dtype=np.longdouble)),
            "signed_orbital_angular_momentum_rate_from_advective_flux": float(signed_orbital_torque),
            "absolute_face_orbital_torque_rate_sum": float(absolute_orbital_torque),
            "one_step_signed_orbital_angular_momentum_fraction": float(DT*signed_orbital_torque/disk_angular_momentum),
            "one_step_absolute_face_torque_over_disk_angular_momentum": float(DT*absolute_orbital_torque/abs(disk_angular_momentum)),
            "max_pair_orbital_plus_spin_torque_residual": float(np.max(np.abs(angular_residual)))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accepted_receipt", type=Path)
    ap.add_argument("static_weights", type=Path)
    ap.add_argument("volume_rates", type=Path)
    ap.add_argument("material_rates", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads(args.accepted_receipt.read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError("Expected accepted native t=6.5 state")
    for key, hash_key in (("initial", "initial_sha256"),
                          ("endpoint", "endpoint_sha256")):
        if digest(Path(receipt[key])) != receipt[hash_key]:
            raise ValueError(f"Changed native {key} snapshot")
    p0, _ = ordered(Path(receipt["initial"]), N)
    p, _ = ordered(Path(receipt["endpoint"]), N)
    xy0 = np.column_stack((p0["x"].astype(float)-3., p0["y"].astype(float)-3.))
    xy = np.column_stack((p["x"].astype(float)-3., p["y"].astype(float)-3.))
    velocity = np.column_stack((p["vx"].astype(float),p["vy"].astype(float)))
    r0 = np.linalg.norm(xy0, axis=1)
    active = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    arrays = []
    for path, key in ((args.static_weights,"weight"),
                      (args.volume_rates,"weight_rate"),
                      (args.material_rates,"material_rate")):
        with np.load(path) as saved:
            if not np.array_equal(saved["active"],active):
                raise ValueError(f"Wrong active IDs: {path}")
            arrays.append(saved[key].copy())
    weight, volume_rate, material_rate = [np.zeros(N) for _ in range(3)]
    weight[active],volume_rate[active],material_rate[active] = arrays
    area, meta, triangles, centers = power_areas(
        xy,weight,active,lookup,return_geometry=True)
    if (meta["hidden_active"] or np.any(area <= 0)
            or np.max(np.abs(area/p0["volume"][active]-1)) > 1e-7):
        raise ValueError("Changed static target-area solution")
    i,j,c,length = finite_faces(triangles,centers,lookup >= 0)
    d = np.linalg.norm(xy[j]-xy[i],axis=1)
    n = (xy[j]-xy[i])/d[:,None]
    u_fixed = (np.einsum("ij,ij->i",xy[j],velocity[j])
               -np.einsum("ij,ij->i",xy[i],velocity[i])
               -np.einsum("ij,ij->i",velocity[j]-velocity[i],c))/d
    vc = np.column_stack((-c[:,1],c[:,0]))
    vc *= (np.einsum("ij,ij->i",c,c)+EPS*EPS)[:,None]**(-.75)
    material_n = np.einsum("ij,ij->i",vc,n)
    mass = p["mass"].astype(float)
    density = p["den"].astype(float)
    density[active] = mass[active]/area
    disk_angular_momentum = np.sum(mass[active]*cross(xy[active],velocity[active]),
                                   dtype=np.longdouble)
    results = {}
    for label,rate in (("fixed_weight",np.zeros(N)),
                       ("volume_rate",volume_rate),
                       ("material_least_squares",material_rate)):
        u_face = u_fixed+(rate[i]-rate[j])/(2*d)
        results[label] = measure(i,j,c,xy,vc,length,density,mass,
                                 disk_angular_momentum,material_n-u_face,
                                 active,label)
    result={"accepted_receipt":str(args.accepted_receipt.resolve()),
            "endpoint_sha256":receipt["endpoint_sha256"],
            "static_weights_sha256":digest(args.static_weights),
            "volume_rates_sha256":digest(args.volume_rates),
            "material_rates_sha256":digest(args.material_rates),
            "time":6.5,"one_step_dt":DT,
            "face_density":"Arithmetic mean of cell density; disk density reconstructed from fixed mass/solved static power area.",
            "material_velocity":"Analytic softened Kepler velocity at the face centroid.",
            "advected_momentum":"Mass flux times analytic face-centroid Kepler velocity; pressure and Riemann stress are excluded. Orbital torque uses site lever arms, and equal/opposite spin torque uses face-centroid lever arms.",
            "scope":"Disk-disk face transport proxy only. It excludes disk/core and disk/exterior faces, native reconstructed face states, and any actual hydro evolution.",
            **results}
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps({k:result[k]["absolute_edge_mass_transfer_rate_over_disk_mass"]
                      for k in results}))


if __name__ == "__main__":
    main()
