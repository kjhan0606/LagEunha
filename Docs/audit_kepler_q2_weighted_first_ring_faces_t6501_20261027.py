"""Audit selected native weighted first-ring face pressure forces and torque."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_faceforce_pair_t6501_20261027 import (
    gravity_only, rows, weighted_summary)
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


CATEGORIES = ("same_ring", "cross_ring", "core_interface")
FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")


def cross(a, b):
    return a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]


def weighted_rms(v, mass):
    w = mass/mass.sum()
    return float(np.sqrt(np.sum(w*v*v)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("face_run", type=Path)
    parser.add_argument("baseline_receipt", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run = args.face_run
    manifest = json.loads((run/"stage_manifest.json").read_text())
    audit = json.loads((run/"audit_t6501.json").read_text())
    baseline = json.loads(args.baseline_receipt.read_text())
    if (len(audit["checks"]) != 21 or not all(audit["checks"].values())
            or audit["parent_time"] != 6.5 or audit["time"] != 6.501
            or audit["parent_sha256"] != manifest["parent_endpoint_sha256"]
            or manifest["source_endpoint_sha256"] != baseline["endpoint_sha256"]):
        raise ValueError("Accepted one-step provenance or physics gate failed")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(audit[key])) != audit[sha_key]:
            raise ValueError(f"Changed {key} snapshot")
    for name, sha in manifest["inputs"].items():
        if digest(run/name) != sha:
            raise ValueError(f"Changed staged input {name}")
    n = audit["n"]
    p0, _ = ordered(Path(audit["initial"]), n)
    start, _ = ordered(Path(audit["parent"]), n)
    end, _ = ordered(Path(audit["endpoint"]), n)
    no_dump, _ = ordered(Path(baseline["endpoint"]), n)
    equal = {field: bool(np.array_equal(end[field], no_dump[field])) for field in FIELDS}
    if not all(equal.values()):
        raise ValueError("Face-dump instrumentation changed a physics field")
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    selected = manifest["selected_ids"]
    if (len(selected) != 64 or len(set(selected)) != 64
            or not all(.5 <= r0[i] < .52 for i in selected)):
        raise ValueError("Unexpected selected first-ring IDs")
    input_path = run/"face_input_r0_call1.tsv"
    face_path = run/"face_dump_r0_call1.tsv"
    inputs = {int(row["index"]): row for row in rows(input_path)}
    faces = rows(face_path)
    forces = {category: {i: np.zeros(2) for i in selected}
              for category in CATEGORIES}
    counts = {category: 0 for category in CATEGORIES}
    selected_set = set(selected)
    face_count = {i: 0 for i in selected}
    max_state_gap = 0.
    max_pressure_gap = 0.
    max_weight_gap = 0.
    max_pstar_gap = 0.
    for row in faces:
        i, j = int(row["i"]), int(row["j"])
        owner = int(row["id_i"])
        if owner not in selected_set or j < 0 or int(row["wall"]) != 0:
            raise ValueError("Unexpected selected owner or ghost/wall face")
        if i not in inputs or j not in inputs:
            raise ValueError("Face input missing owner or neighbor")
        a, b = inputs[i], inputs[j]
        if int(a["id"]) != owner or int(b["id"]) != int(row["id_j"]):
            raise ValueError("Face IDs differ from captured inputs")
        for key in ("x", "y", "vx", "vy"):
            max_state_gap = max(max_state_gap,
                                abs(float(a[key])-float(start[key][owner])))
        max_pressure_gap = max(max_pressure_gap,
                               abs(float(a["P"])-float(start["pressure"][owner])))
        max_weight_gap = max(max_weight_gap,
                             abs(float(a["W"])-float(start["w2"][owner])))
        max_pstar_gap = max(max_pstar_gap, abs(float(row["pstar"])-float(a["P"])))
        category = ("core_interface" if r0[int(row["id_j"])] < .5-1e-8
                    else "same_ring" if abs(r0[owner]-r0[int(row["id_j"])]) < 1e-9
                    else "cross_ring")
        xi = np.array((float(a["x"]), float(a["y"])))
        xj = np.array((float(b["x"]), float(b["y"])))
        x1 = np.array((float(row["x1"]), float(row["y1"])))
        x2 = np.array((float(row["x2"]), float(row["y2"])))
        ds = np.array((float(row["dSx"]), float(row["dSy"])))
        if np.max(np.abs(ds-np.array((x2[1]-x1[1], -(x2[0]-x1[0]))))) > 1e-13:
            raise ValueError("Face area vector differs from corners")
        distance = np.linalg.norm(xj-xi)
        if distance <= 0:
            raise ValueError("Zero generator separation")
        coeff = .5*ds+np.linalg.norm(ds)*(.5*(x1+x2-xi-xj))/distance
        forces[category][owner] += (float(a["P"])-float(b["P"]))*coeff
        counts[category] += 1
        face_count[owner] += 1
    if (min(face_count.values()) < 3 or max_state_gap != 0
            or max_pressure_gap > 2e-13 or max_weight_gap != 0
            or max_pstar_gap != 0):
        raise ValueError("First-call face inputs differ from native restart")

    ids = np.asarray(selected, dtype=int)
    mass = start["mass"][ids].astype(float)
    centre_r = np.column_stack((start["x"][ids].astype(float)-3.,
                                start["y"][ids].astype(float)-3.))
    start_v = np.column_stack((start["vx"][ids].astype(float),
                               start["vy"][ids].astype(float)))
    end_v = np.column_stack((end["vx"][ids].astype(float),
                             end["vy"][ids].astype(float)))
    reference_v = np.empty_like(end_v)
    for k, i in enumerate(ids):
        state = (centre_r[k, 0], centre_r[k, 1],
                 start_v[k, 0], start_v[k, 1])
        reference_v[k] = gravity_only(state, .001)[2:]
    observed_acc = (end_v-reference_v)/.001
    face_acc = {category: np.array([forces[category][int(i)]/mass[k]
                                    for k, i in enumerate(ids)])
                for category in CATEGORIES}
    total_acc = sum(face_acc.values())
    observed_torque = cross(centre_r, observed_acc)
    total_torque = cross(centre_r, total_acc)
    obs_norm2 = weighted_rms(observed_torque, mass)**2
    categories = {}
    for category in CATEGORIES:
        torque = cross(centre_r, face_acc[category])
        categories[category] = {
            "directed_faces": counts[category],
            "pressure_acceleration_rms": float(np.sqrt(np.sum(
                (mass/mass.sum())*np.sum(face_acc[category]**2, axis=1)))),
            "orbital_torque_rms": weighted_rms(torque, mass),
            "orbital_torque_projection_on_observed_fraction": float(
                np.sum((mass/mass.sum())*torque*observed_torque)/obs_norm2),
        }
    result = {"face_run": str(run), "baseline_receipt": str(args.baseline_receipt),
              "start_sha256": audit["parent_sha256"],
              "endpoint_sha256": audit["endpoint_sha256"],
              "face_dump_sha256": digest(face_path),
              "face_input_sha256": digest(input_path),
              "selected_cells": len(ids), "directed_face_rows": len(faces),
              "physics_fields_bitwise_equal_to_no_dump": equal,
              "max_first_call_state_gap": max_state_gap,
              "max_first_call_pressure_gap": max_pressure_gap,
              "max_first_call_weight_gap": max_weight_gap,
              "max_pstar_minus_owner_pressure": max_pstar_gap,
              "pressure_acceleration": weighted_summary(total_acc, observed_acc, mass),
              "observed_orbital_torque_rms": weighted_rms(observed_torque, mass),
              "reconstructed_face_orbital_torque_rms": weighted_rms(total_torque, mass),
              "torque_relative_l2_difference": (weighted_rms(
                  observed_torque-total_torque, mass)/weighted_rms(observed_torque, mass)),
              "categories": categories,
              "scope": ("64 selected first-ring cells, one accepted weighted GPU "
                        "RK4 first-call face dump. Local fixed-weight geometric "
                        "pressure-force reconstruction may omit nonlocal weight "
                        "projection; agreement must be checked against the "
                        "gravity-subtracted native one-step acceleration. Not "
                        "an accumulated causal attribution or a full-disk force audit.")}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print("relative_force_l2", result["pressure_acceleration"]["nongrav_minus_face_relative_l2"])
    print("relative_torque_l2", result["torque_relative_l2_difference"])
    for category in CATEGORIES:
        print(category, categories[category])


if __name__ == "__main__":
    main()
