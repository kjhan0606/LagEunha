"""Reconstruct native W0 face pressure forces and a one-step orbit defect."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


N = 55351
EXPECTED_BINARY = "3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def gravity_only(position_velocity, dt):
    def rhs(_, state):
        x, y, vx, vy = state
        inv = (x*x+y*y+.01175**2)**(-1.5)
        return (vx, vy, -x*inv, -y*inv)
    solution = solve_ivp(rhs, (0., dt), position_velocity,
                         method="DOP853", rtol=1e-12, atol=1e-14)
    if not solution.success:
        raise ValueError("Gravity-only reference integration failed")
    return solution.y[:, -1]


def weighted_summary(face, observed, mass):
    w = mass/np.sum(mass)
    fmag = np.linalg.norm(face, axis=1)
    omag = np.linalg.norm(observed, axis=1)
    diff = observed-face
    f2 = float(np.sum(w*fmag*fmag))
    o2 = float(np.sum(w*omag*omag))
    d2 = float(np.sum(w*np.sum(diff*diff, axis=1)))
    cross = float(np.sum(w*np.sum(face*observed, axis=1)))
    return {
        "cells": len(mass),
        "face_pressure_acceleration_mean_abs": float(np.sum(w*fmag)),
        "face_pressure_acceleration_rms": float(np.sqrt(f2)),
        "nongrav_acceleration_mean_abs": float(np.sum(w*omag)),
        "nongrav_acceleration_rms": float(np.sqrt(o2)),
        "nongrav_minus_face_acceleration_rms": float(np.sqrt(d2)),
        "nongrav_minus_face_relative_l2": float(np.sqrt(d2/o2)) if o2 else None,
        "vector_cosine": float(cross/np.sqrt(f2*o2)) if f2 and o2 else None,
        "nongrav_on_face_slope": float(cross/f2) if f2 else None,
    }


def evaluate_case(root, label, selected, highest):
    directory = root/label
    audit = json.loads((directory/"audit_t6501.json").read_text())
    if len(audit["checks"]) != 11 or not all(audit["checks"].values()):
        raise ValueError(f"{label}: one-step native physics gate failed")
    if digest(directory/"eunha2.exe") != EXPECTED_BINARY:
        raise ValueError(f"{label}: native binary changed")
    for name, sha in (("restart", "restart_sha256"), ("endpoint", "endpoint_sha256")):
        if digest(Path(audit[name])) != audit[sha]:
            raise ValueError(f"{label}: {name} checkpoint changed")
    dt = float(audit["endpoint_time"]-audit["start_time"])
    if abs(audit["start_time"]-6.5)>1e-10 or abs(dt-.001)>1e-10:
        raise ValueError(f"{label}: wrong force-probe time")
    start, _ = ordered(Path(audit["restart"]), N)
    end, _ = ordered(Path(audit["endpoint"]), N)

    inputs_path = directory/"face_input_r0_call1.tsv"
    faces_path = directory/"face_dump_r0_call1.tsv"
    inputs = {int(row["index"]): row for row in rows(inputs_path)}
    faces = rows(faces_path)
    face_force = {cell: np.zeros(2, dtype=float) for cell in selected}
    face_count = {cell: 0 for cell in selected}
    max_pressure_discrepancy = 0.
    max_state_discrepancy = 0.
    for row in faces:
        i, j = int(row["i"]), int(row["j"])
        owner = int(row["id_i"])
        if owner not in face_force or j < 0 or int(row["wall"]) != 0:
            raise ValueError(f"{label}: unexpected owner, ghost, or wall face")
        if i not in inputs or j not in inputs:
            raise ValueError(f"{label}: face input missing owner or neighbor")
        a, b = inputs[i], inputs[j]
        if int(a["id"]) != owner or int(b["id"]) != int(row["id_j"]):
            raise ValueError(f"{label}: face ID differs from captured input")
        for key, field in (("x", "x"), ("y", "y"), ("vx", "vx"), ("vy", "vy")):
            max_state_discrepancy = max(max_state_discrepancy,
                                        abs(float(a[key])-float(start[field][owner])))
        pi, pj = float(a["P"]), float(b["P"])
        max_pressure_discrepancy = max(max_pressure_discrepancy,
                                       abs(float(row["pstar"])-pi))
        if float(a["W"]) != 0 or float(b["W"]) != 0:
            raise ValueError(f"{label}: nonzero Laguerre weight")
        xi = np.array((float(a["x"]), float(a["y"])))
        xj = np.array((float(b["x"]), float(b["y"])))
        x1 = np.array((float(row["x1"]), float(row["y1"])))
        x2 = np.array((float(row["x2"]), float(row["y2"])))
        ds = np.array((float(row["dSx"]), float(row["dSy"])))
        if np.max(np.abs(ds-np.array((x2[1]-x1[1], -(x2[0]-x1[0]))))) > 1e-13:
            raise ValueError(f"{label}: face area vector differs from corners")
        distance = np.linalg.norm(xj-xi)
        if distance <= 0:
            raise ValueError(f"{label}: zero generator distance")
        offset = .5*(x1+x2-xi-xj)
        b_i = .5*ds + np.linalg.norm(ds)*offset/(distance+1e-30)
        face_force[owner] += (pi-pj)*b_i
        face_count[owner] += 1
    if set(face_force) != set(selected) or min(face_count.values()) < 3:
        raise ValueError(f"{label}: missing selected-cell faces")
    if max_pressure_discrepancy != 0 or max_state_discrepancy != 0:
        raise ValueError(f"{label}: first face input differs from restart")

    face_acc = []
    nongrav_acc = []
    mass = []
    cell_rows = []
    for cell in selected:
        m = float(start["mass"][cell])
        acceleration = face_force[cell]/m
        state0 = (float(start["x"][cell])-3., float(start["y"][cell])-3.,
                  float(start["vx"][cell]), float(start["vy"][cell]))
        gravity_end = gravity_only(state0, dt)
        observed = (np.array((float(end["vx"][cell]), float(end["vy"][cell])))
                    - gravity_end[2:])/dt
        face_acc.append(acceleration)
        nongrav_acc.append(observed)
        mass.append(m)
        cell_rows.append({"id": cell, "group": "high_pressure_change" if cell in highest else "random_inner",
                          "face_count": face_count[cell], "face_pressure_acceleration": acceleration.tolist(),
                          "nongrav_acceleration": observed.tolist(), "mass": m})
    face_acc = np.asarray(face_acc)
    nongrav_acc = np.asarray(nongrav_acc)
    mass = np.asarray(mass)
    high_mask = np.array([cell in highest for cell in selected])
    return {
        "restart_sha256": audit["restart_sha256"],
        "endpoint_sha256": audit["endpoint_sha256"],
        "face_dump_sha256": digest(faces_path),
        "face_input_sha256": digest(inputs_path),
        "dt": dt,
        "face_rows": len(faces), "face_input_rows": len(inputs),
        "max_first_call_pressure_discrepancy": max_pressure_discrepancy,
        "max_first_call_state_discrepancy": max_state_discrepancy,
        "all_selected": weighted_summary(face_acc, nongrav_acc, mass),
        "high_pressure_change": weighted_summary(face_acc[high_mask], nongrav_acc[high_mask], mass[high_mask]),
        "random_inner": weighted_summary(face_acc[~high_mask], nongrav_acc[~high_mask], mass[~high_mask]),
        "cells": cell_rows,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("staged_pair", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = args.staged_pair.resolve()
    manifest = json.loads((root/"stage_manifest.json").read_text())
    if digest(root/"selected_ids.txt") != manifest["selected_ids_sha256"]:
        raise ValueError("Selected ID list changed")
    selected = list(manifest["selected_ids"])
    highest = set(manifest["highest_pressure_change_ids"])
    if len(selected) != 64 or len(highest) != 32:
        raise ValueError("Wrong selected-face sample")
    cases = {label: evaluate_case(root, label, selected, highest)
             for label in ("aligned", "jitter")}
    ratio = (cases["jitter"]["all_selected"]["face_pressure_acceleration_rms"]
             / cases["aligned"]["all_selected"]["face_pressure_acceleration_rms"])
    result = {
        "staged_pair": str(root),
        "selected_ids_sha256": manifest["selected_ids_sha256"],
        "binary_sha256": EXPECTED_BINARY,
        "cases": cases,
        "jitter_over_aligned_face_acceleration_rms": ratio,
        "scope": ("64 preselected inner-disk cells, one accepted-state W0 face-force call and one "
                  "6.5→6.501 native step. Gravity-only ODE reference isolates effective nongrav velocity change; "
                  "it does not prove sustained or whole-disk behavior."),
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for label in ("aligned", "jitter"):
        print(label, json.dumps(cases[label]["all_selected"], indent=2))
    print("jitter_over_aligned_face_acceleration_rms", ratio)


if __name__ == "__main__":
    main()
