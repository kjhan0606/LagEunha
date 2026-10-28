"""Rebuild selected W0 face forces using exact-orbit adiabatic pressures."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_faceforce_pair_t6501_20261027 import rows
from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
GAMMA = float(np.float32(1.666667))
PARENTS = {
    "aligned": DOCS/"KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
    "jitter": DOCS/"KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json",
}


def compare(predicted, native, mass):
    w = mass/np.sum(mass)
    pred2 = float(np.sum(w*np.sum(predicted*predicted, axis=1)))
    nat2 = float(np.sum(w*np.sum(native*native, axis=1)))
    diff2 = float(np.sum(w*np.sum((native-predicted)**2, axis=1)))
    cross = float(np.sum(w*np.sum(predicted*native, axis=1)))
    return {
        "cells": len(mass),
        "ballistic_pressure_force_acceleration_rms": float(np.sqrt(pred2)),
        "native_pressure_force_acceleration_rms": float(np.sqrt(nat2)),
        "native_minus_ballistic_acceleration_rms": float(np.sqrt(diff2)),
        "native_minus_ballistic_relative_l2": float(np.sqrt(diff2/nat2)) if nat2 else None,
        "vector_cosine": float(cross/np.sqrt(pred2*nat2)) if pred2 and nat2 else None,
        "native_on_ballistic_slope": float(cross/pred2) if pred2 else None,
    }


def evaluate(label, root, parent, original, selected, high):
    accepted = json.loads(parent.read_text())
    if not all(accepted["checks"].values()) or abs(accepted["time"]-6.5)>1e-10:
        raise ValueError(f"{label}: accepted parent failed")
    if digest(Path(accepted["initial"])) != accepted["initial_sha256"]:
        raise ValueError(f"{label}: accepted initial checkpoint changed")
    p0, _ = ordered(Path(accepted["initial"]), 55351)
    directory = root/label
    face_path = directory/"face_dump_r0_call1.tsv"
    input_path = directory/"face_input_r0_call1.tsv"
    if digest(face_path) != original["face_dump_sha256"] or digest(input_path) != original["face_input_sha256"]:
        raise ValueError(f"{label}: selected-face capture changed")
    face_rows = rows(face_path)
    inputs = {int(q["index"]):q for q in rows(input_path)}
    ids = sorted({int(q["id"]) for q in inputs.values()})
    if len(ids) != len(inputs):
        raise ValueError(f"{label}: repeated face-input ID")
    _, area, safety = disk_geometry(p0, 6.5, np.asarray(ids, dtype=int))
    if safety <= 1:
        raise ValueError(f"{label}: exact-orbit area unsafe")
    pball = {cell:float(p0["pressure"][cell]) *
             (float(p0["volume"][cell])/a)**GAMMA
             for cell,a in zip(ids,area)}
    predicted = {cell:np.zeros(2) for cell in selected}
    for q in face_rows:
        i,j,owner = int(q["i"]),int(q["j"]),int(q["id_i"])
        if owner not in predicted or i not in inputs or j not in inputs:
            raise ValueError(f"{label}: selected-face inconsistency")
        owner_id,neighbor_id=int(inputs[i]["id"]),int(inputs[j]["id"])
        if owner != owner_id or int(q["id_j"]) != neighbor_id:
            raise ValueError(f"{label}: selected-face ID mismatch")
        xi = np.array((float(inputs[i]["x"]),float(inputs[i]["y"])))
        xj = np.array((float(inputs[j]["x"]),float(inputs[j]["y"])))
        x1 = np.array((float(q["x1"]),float(q["y1"])))
        x2 = np.array((float(q["x2"]),float(q["y2"])))
        ds = np.array((float(q["dSx"]),float(q["dSy"])))
        distance=np.linalg.norm(xj-xi)
        offset=.5*(x1+x2-xi-xj)
        b_i=.5*ds+np.linalg.norm(ds)*offset/(distance+1e-30)
        predicted[owner]+=(pball[owner_id]-pball[neighbor_id])*b_i/float(inputs[i]["mass"])
    original_cells={int(q["id"]):q for q in original["cells"]}
    pred=np.asarray([predicted[cell] for cell in selected])
    native=np.asarray([original_cells[cell]["face_pressure_acceleration"] for cell in selected])
    mass=np.asarray([float(original_cells[cell]["mass"]) for cell in selected])
    high_mask=np.asarray([cell in high for cell in selected])
    pressure_difference=np.asarray([float(inputs[i]["P"])-pball[int(inputs[i]["id"])] for i in inputs])
    pressure_native=np.asarray([float(inputs[i]["P"]) for i in inputs])
    return {
        "initial_sha256":accepted["initial_sha256"],
        "face_dump_sha256":original["face_dump_sha256"],
        "face_input_sha256":original["face_input_sha256"],
        "ballistic_area_image_safety_ratio":safety,
        "captured_input_cells":len(ids),
        "captured_input_pressure_relative_residual_rms":float(np.sqrt(np.mean((pressure_difference/pressure_native)**2))),
        "all_selected":compare(pred,native,mass),
        "high_pressure_change":compare(pred[high_mask],native[high_mask],mass[high_mask]),
        "random_inner":compare(pred[~high_mask],native[~high_mask],mass[~high_mask]),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("force_pair_receipt",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    original=json.loads(args.force_pair_receipt.read_text())
    root=Path(original["staged_pair"])
    manifest=json.loads((root/"stage_manifest.json").read_text())
    selected=list(manifest["selected_ids"])
    high=set(manifest["highest_pressure_change_ids"])
    cases={label:evaluate(label,root,PARENTS[label],original["cases"][label],selected,high)
           for label in ("aligned","jitter")}
    result={
        "force_pair_receipt":str(args.force_pair_receipt.resolve()),
        "gamma":GAMMA,
        "cases":cases,
        "scope":("Exact softened circular-site Voronoi areas set adiabatic P_i(t); the native "
                 "t=6.5 face geometry and W0 force coefficients are held fixed. This is an "
                 "instantaneous force decomposition on 64 selected cells, not a counterfactual evolution."),
    }
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    for label in ("aligned","jitter"):
        print(label,json.dumps(cases[label]["all_selected"],indent=2))


if __name__=="__main__":
    main()
