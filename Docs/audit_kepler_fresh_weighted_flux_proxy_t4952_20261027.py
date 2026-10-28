"""Estimate face material-mass residual in a same-state W0/power comparison."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_weighted_face_slip_t65_20261027 import EPS
from audit_kepler_sparse_area_t0825_20261027 import ordered
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
BINS = ((.5, .52), (.52, .55), (.55, .65), (.5, .65),
        (.65, .75), (.75, 1.), (1., 2.), (.5, 2.))


def touching_finite_faces(triangles, centers, active_flag):
    edges = np.vstack((triangles[:, [0, 1]], triangles[:, [1, 2]],
                       triangles[:, [2, 0]]))
    triangle_ids = np.tile(np.arange(len(triangles)), 3)
    edges.sort(axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, triangle_ids = edges[order], triangle_ids[order]
    starts = np.r_[0, 1+np.flatnonzero(np.any(edges[1:] != edges[:-1], axis=1))]
    counts = np.diff(np.r_[starts, len(edges)])
    touching = active_flag[edges[starts, 0]] | active_flag[edges[starts, 1]]
    if np.any(counts[touching] != 2):
        raise ValueError("Non-finite or nonmanifold face touches active disk")
    starts = starts[touching]
    i, j = edges[starts, 0], edges[starts, 1]
    a, b = centers[triangle_ids[starts]], centers[triangle_ids[starts+1]]
    midpoint = .5*(a+b)
    length = np.linalg.norm(a-b, axis=1)
    if np.any(~np.isfinite(length) | (length <= 0)):
        raise ValueError("Invalid finite face length")
    return i, j, midpoint, length


def summarize_flux(i, j, length, signed_flux, mass, active, r0, dt):
    cell_rate = np.zeros(N)
    np.add.at(cell_rate, i, signed_flux)
    np.add.at(cell_rate, j, -signed_flux)
    gross = float(np.sum(np.abs(signed_flux), dtype=np.longdouble))
    if abs(float(np.sum(cell_rate, dtype=np.longdouble))) > gross*1e-12:
        raise ValueError("Pairwise mass flux fails closure")
    rows = {}
    for lo, hi in BINS:
        take = (r0 >= lo) & (r0 < hi)
        fraction = dt*np.abs(cell_rate[take])/mass[take]
        rows[f"{lo:g}:{hi:g}"] = {
            "cells": int(take.sum()),
            "median_one_step_abs_net_cell_mass_fraction": float(np.median(fraction)),
            "p99_one_step_abs_net_cell_mass_fraction": float(np.quantile(fraction, .99)),
            "max_one_step_abs_net_cell_mass_fraction": float(np.max(fraction)),
        }
    disk_mass = float(np.sum(mass[active], dtype=np.longdouble))
    return {"faces": len(i), "gross_mass_flux_rate": gross,
            "one_step_gross_mass_flux_over_disk_mass": dt*gross/disk_mass,
            "net_pairwise_mass_flux_rate": float(np.sum(cell_rate, dtype=np.longdouble)),
            "bins": rows}


def face_metrics(label, xy, velocity, weight, weight_rate, density, mass,
                 active, lookup, r0, dt):
    _, meta, triangles, centers = power_areas(
        xy, weight, active, lookup, return_geometry=True)
    if meta["hidden_active"]:
        raise ValueError(f"{label}: hidden active cell")
    active_flag = lookup >= 0
    i, j, midpoint, length = touching_finite_faces(triangles, centers, active_flag)
    both_disk = active_flag[i] & active_flag[j]
    boundary = ~both_disk
    if int(np.count_nonzero(both_disk)) != 145671:
        raise ValueError(f"{label}: unexpected disk-disk face count")
    separation = xy[j] - xy[i]
    d = np.linalg.norm(separation, axis=1)
    normal = separation / d[:, None]
    face_speed = (np.einsum("ij,ij->i", xy[j], velocity[j])
                  - np.einsum("ij,ij->i", xy[i], velocity[i])
                  - np.einsum("ij,ij->i", velocity[j]-velocity[i], midpoint)) / d
    face_speed += (weight_rate[i]-weight_rate[j])/(2*d)
    kepler_velocity = np.column_stack((-midpoint[:, 1], midpoint[:, 0]))
    kepler_velocity *= (np.einsum("ij,ij->i", midpoint, midpoint)
                        + EPS*EPS)[:, None]**(-.75)
    material_speed = np.einsum("ij,ij->i", kepler_velocity, normal)
    rho_face = .5*(density[i]+density[j])
    signed_flux = rho_face*length*(material_speed-face_speed)
    disk = summarize_flux(i[both_disk], j[both_disk], length[both_disk],
                          signed_flux[both_disk], mass, active, r0, dt)
    all_faces = summarize_flux(i, j, length, signed_flux, mass, active, r0, dt)
    other = np.where(active_flag[i], j, i)
    return {"disk_disk": disk, "all_faces_touching_disk": all_faces,
            "disk_core_faces": int(np.count_nonzero(boundary & (r0[other] < .5))),
            "disk_exterior_faces": int(np.count_nonzero(boundary & (r0[other] >= 2.))),
            "max_abs_face_speed_slip": float(np.max(np.abs(material_speed-face_speed))),
            "max_abs_disk_disk_face_speed_slip": float(np.max(
                np.abs(material_speed[both_disk]-face_speed[both_disk])))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text())
    if len(receipt["checks"]) != 21 or not all(receipt["checks"].values()):
        raise ValueError("Native one-step gate is not accepted")
    for key, hash_key in (("initial", "initial_sha256"),
                          ("parent", "parent_sha256"),
                          ("endpoint", "endpoint_sha256")):
        if digest(Path(receipt[key])) != receipt[hash_key]:
            raise ValueError(f"{key} snapshot digest differs")
    initial, _ = ordered(Path(receipt["initial"]), N)
    before, mb = ordered(Path(receipt["parent"]), N)
    after, ma = ordered(Path(receipt["endpoint"]), N)
    dt = float(ma["t"] - mb["t"])
    if abs(dt-.001) > 1e-10 or not np.array_equal(before["mass"], after["mass"]):
        raise ValueError("Expected one 0.001 step with fixed cell masses")
    active_mask = np.fromfile(args.receipt.parent/"disk_active_mask.bin", dtype=np.uint8)
    if len(active_mask) != N or int(np.count_nonzero(active_mask)) != 49038:
        raise ValueError("Unexpected active disk mask")
    active = np.flatnonzero(active_mask)
    lookup = np.full(N, -1, dtype=int)
    lookup[active] = np.arange(len(active))
    r0 = np.hypot(initial["x"].astype(float)-3., initial["y"].astype(float)-3.)
    xy = .5*np.column_stack((before["x"].astype(float)+after["x"].astype(float)-6.,
                              before["y"].astype(float)+after["y"].astype(float)-6.))
    velocity = .5*np.column_stack((before["vx"].astype(float)+after["vx"].astype(float),
                                    before["vy"].astype(float)+after["vy"].astype(float)))
    density = .5*(before["den"].astype(float)+after["den"].astype(float))
    mass = before["mass"].astype(float)
    weight = .5*(before["w2"].astype(float)+after["w2"].astype(float))
    rate = (after["w2"].astype(float)-before["w2"].astype(float))/dt
    zero = np.zeros(N)
    w0 = face_metrics("W0", xy, velocity, zero, zero, density, mass,
                      active, lookup, r0, dt)
    moving = face_metrics("moving_power", xy, velocity, weight, rate,
                          density, mass, active, lookup, r0, dt)
    result = {"receipt": str(args.receipt), "receipt_sha256": digest(args.receipt),
              "parent_sha256": receipt["parent_sha256"],
              "endpoint_sha256": receipt["endpoint_sha256"],
              "start_time": float(mb["t"]), "end_time": float(ma["t"]),
              "dt": dt, "W0_same_state": w0, "moving_power": moving,
              "scope": ("Midpoint independent disk-disk power faces, arithmetic "
                        "cell density at each face and analytic softened Kepler "
                        "velocity. Proxy only; no native Riemann states or "
                        "disk/core and disk/exterior faces. W0 uses the same "
                        "weighted-run generator state, not a separate W0 trajectory.")}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for label in ("W0_same_state", "moving_power"):
        row = result[label]["all_faces_touching_disk"]
        print(label, "gross_fraction", row["one_step_gross_mass_flux_over_disk_mass"],
              "disk_p99_net_cell_fraction",
              row["bins"]["0.5:2"]["p99_one_step_abs_net_cell_mass_fraction"],
              "first_ring_p99", row["bins"]["0.5:0.52"]["p99_one_step_abs_net_cell_mass_fraction"])


if __name__ == "__main__":
    main()
