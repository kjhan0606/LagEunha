"""Reconstruct the accepted GPU W=0 face pressure force for matched sites."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


EPS = 0.01175


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def records(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty native capture: {path}")
    return rows


def reconstruct(run, ids):
    inputs = {int(row["index"]): row for row in records(run / "face_input_r0_call1.tsv")}
    faces = records(run / "face_dump_r0_call1.tsv")
    meta = json.loads((run / "face_input_meta_r0_call1.json").read_text())
    if (meta["call"] != 1 or meta["phase1"] != 1 or meta["av_mode"] != 5
            or meta["use_muscl"] != 1 or meta["nu_phys"] != 0):
        raise ValueError("Native capture is not the accepted smooth W=0 branch")
    forces = {i: np.zeros(2) for i in ids}
    sums = {i: np.zeros(2) for i in ids}
    counts = {i: 0 for i in ids}
    max_work_error = max_pressure_error = max_geometry_error = 0.0
    for row in faces:
        owner = int(row["id_i"])
        if (owner not in ids or int(row["wall"]) != 0
                or int(row["j"]) < 0 or int(row["call"]) != 1):
            raise ValueError("Unexpected owner, wall, ghost, or force call")
        i, j = int(row["i"]), int(row["j"])
        a, b = inputs[i], inputs[j]
        if int(a["id"]) != owner or int(b["id"]) != int(row["id_j"]):
            raise ValueError("Face/input owner or neighbor differs")
        xi = np.array([float(a["x"]), float(a["y"])])
        xj = np.array([float(b["x"]), float(b["y"])])
        vi = np.array([float(a["vx"]), float(a["vy"])])
        ds = np.array([float(row["dSx"]), float(row["dSy"])])
        x1 = np.array([float(row["x1"]), float(row["y1"])])
        x2 = np.array([float(row["x2"]), float(row["y2"])])
        length = np.linalg.norm(ds)
        separation = np.linalg.norm(xj - xi)
        if length <= 0 or separation <= 0:
            raise ValueError("Degenerate native face")
        face_center = 0.5 * (x1 + x2)
        b_i = 0.5 * ds + length * (face_center - 0.5 * (xi + xj)) / separation
        pi, pj = float(a["P"]), float(b["P"])
        f = (pi - pj) * b_i
        predicted_work = float(vi @ f) - pi * float(row["udotds"])
        max_work_error = max(max_work_error, abs(predicted_work - float(row["work"])))
        max_pressure_error = max(max_pressure_error, abs(pi - float(row["pstar"])))
        max_geometry_error = max(max_geometry_error,
                                 np.linalg.norm(ds - np.array([x2[1]-x1[1], x1[0]-x2[0]])))
        forces[owner] += f
        sums[owner] += ds
        counts[owner] += 1
    by_id = {}
    for id_ in sorted(ids):
        owner_rows = [row for row in inputs.values() if int(row["id"]) == id_]
        if len(owner_rows) != 1 or counts[id_] < 3:
            raise ValueError(f"Missing owner or face rows for {id_}")
        a = owner_rows[0]
        x = np.array([float(a["x"])-3., float(a["y"])-3.])
        r = np.linalg.norm(x)
        m = float(a["mass"])
        if r <= 0 or m <= 0:
            raise ValueError("Invalid native cell radius or mass")
        er = x / r
        et = np.array([-er[1], er[0]])
        acc = forces[id_] / m
        gravity = r / (r*r + EPS*EPS)**1.5
        by_id[id_] = {
            "pressure": float(a["P"]), "radius": float(r), "mass": m,
            "faces": counts[id_], "pressure_acceleration_x": float(acc[0]),
            "pressure_acceleration_y": float(acc[1]),
            "pressure_acceleration_radial": float(acc @ er),
            "pressure_acceleration_tangential": float(acc @ et),
            "pressure_acceleration_norm": float(np.linalg.norm(acc)),
            "pressure_acceleration_over_gravity": float(np.linalg.norm(acc) / gravity),
            "face_vector_closure": float(np.linalg.norm(sums[id_])),
        }
    norms = np.array([v["pressure_acceleration_norm"] for v in by_id.values()])
    checks = {"all_64_selected": len(by_id) == 64,
              "force_finite": bool(np.all(np.isfinite(norms))),
              "face_vector_closed": max(v["face_vector_closure"] for v in by_id.values()) < 1e-10,
              "native_pressure_matches_input": max_pressure_error < 1e-12,
              "native_work_matches_variational_formula": max_work_error < 1e-14,
              "native_area_vectors_match_endpoints": max_geometry_error < 1e-12}
    checks = {key: bool(value) for key, value in checks.items()}
    return {"face_rows": len(faces), "input_rows": len(inputs),
            "face_dump_sha256": digest(run / "face_dump_r0_call1.tsv"),
            "input_dump_sha256": digest(run / "face_input_r0_call1.tsv"),
            "max_face_work_error": max_work_error,
            "max_pressure_input_error": max_pressure_error,
            "max_area_vector_error": max_geometry_error,
            "median_pressure_acceleration": float(np.median(norms)),
            "p90_pressure_acceleration": float(np.quantile(norms, .9)),
            "max_pressure_acceleration": float(np.max(norms)),
            "checks": checks, "cells": by_id}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.run / "staging.json").read_text())
    selected = json.loads((args.run / "selected_ids.json").read_text())
    ids = {row["id"] for row in selected}
    if len(ids) != 64:
        raise ValueError("Wrong selected ID count")
    captures = {}
    for name in ("aligned", "jitter"):
        run = args.run / name
        short = json.loads((run / "audit_short.json").read_text())
        if not all(short["checks"].values()):
            raise ValueError(f"Failed short replay gate: {name}")
        if (digest(run / "khout.000001.dat") != manifest["sources"][name]["accepted_endpoint_sha256"]
                or digest(run / "eunha2.exe") != manifest["sources"][name]["binary_sha256"]):
            raise ValueError(f"Changed accepted input: {name}")
        captures[name] = reconstruct(run, ids)
    aligned = captures["aligned"]["cells"]
    jitter = captures["jitter"]["cells"]
    paired = []
    for row in selected:
        id_ = row["id"]
        a = aligned[id_]
        j = jitter[id_]
        paired.append({**row,
                       "aligned_pressure_acceleration_norm": a["pressure_acceleration_norm"],
                       "jitter_pressure_acceleration_norm": j["pressure_acceleration_norm"],
                       "jitter_over_aligned_acceleration":
                       j["pressure_acceleration_norm"] / a["pressure_acceleration_norm"]
                       if a["pressure_acceleration_norm"] > 0 else None,
                       "acceleration_vector_difference": float(np.hypot(
                           j["pressure_acceleration_x"]-a["pressure_acceleration_x"],
                           j["pressure_acceleration_y"]-a["pressure_acceleration_y"]))})
    differences = np.array([row["acceleration_vector_difference"] for row in paired])
    result = {"run": str(args.run.resolve()), "time": 6.5,
              "auditor_sha256": digest(Path(__file__)),
              "selected_ids_sha256": digest(args.run / "selected_ids.json"),
              "captures": captures, "paired_cells": paired,
              "median_vector_acceleration_difference": float(np.median(differences)),
              "p90_vector_acceleration_difference": float(np.quantile(differences, .9)),
              "checks": {name: all(capture["checks"].values())
                         for name, capture in captures.items()},
              "scope": "One GPU RK4 force call at matched accepted t=6.5 states; selected 64 cells only. Native face inputs reconstruct the variational pressure force. This is not a time-integrated force or causal attribution."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print("aligned median |a_P|", captures["aligned"]["median_pressure_acceleration"])
    print("jitter median |a_P|", captures["jitter"]["median_pressure_acceleration"])
    print("paired median |delta a_P|", result["median_vector_acceleration_difference"])
    if not all(result["checks"].values()):
        raise SystemExit("Native force input audit failed")


if __name__ == "__main__":
    main()
