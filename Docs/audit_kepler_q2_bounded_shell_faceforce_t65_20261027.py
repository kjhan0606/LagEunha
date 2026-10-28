"""Decompose an accepted W0 Kepler force probe by initial polar shell."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_bounded_faceforce_pair_t6501_20261027 import rows
from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
N = 55351
CASES = {
    "aligned": "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json",
    "jitter": "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json",
}


def evaluate(label, run_root, pair_receipt, selected):
    initial_receipt = json.loads((DOCS / CASES[label]).read_text())
    initial, _ = ordered(Path(initial_receipt["initial"]), N)
    r0 = np.hypot(initial["x"].astype(float) - 3., initial["y"].astype(float) - 3.)
    by_id = {int(row["id"]): row for row in pair_receipt["cases"][label]["cells"]}
    assert set(by_id) == set(selected)
    directory = run_root / label
    input_path = directory / "face_input_r0_call1.tsv"
    face_path = directory / "face_dump_r0_call1.tsv"
    assert digest(input_path) == pair_receipt["cases"][label]["face_input_sha256"]
    assert digest(face_path) == pair_receipt["cases"][label]["face_dump_sha256"]
    inputs = {int(row["index"]): row for row in rows(input_path)}
    categories = ("same_shell", "cross_shell", "core_interface")
    forces = {key: {i: np.zeros(2) for i in selected} for key in categories}
    counts = {key: 0 for key in categories}
    lengths = {key: 0. for key in categories}
    max_same_dr = 0.
    min_cross_dr = np.inf
    for row in rows(face_path):
        i = int(row["id_i"])
        j = int(row["id_j"])
        assert i in selected and int(row["wall"]) == 0
        dr = abs(r0[i] - r0[j])
        if min(r0[i], r0[j]) < .5 - 1e-8:
            category = "core_interface"
        elif dr < 1e-9:
            category = "same_shell"
            max_same_dr = max(max_same_dr, dr)
        else:
            category = "cross_shell"
            min_cross_dr = min(min_cross_dr, dr)
        a, b = inputs[int(row["i"])], inputs[int(row["j"])]
        assert int(a["id"]) == i and int(b["id"]) == j
        xi = np.array([float(a["x"]), float(a["y"])])
        xj = np.array([float(b["x"]), float(b["y"])])
        x1 = np.array([float(row["x1"]), float(row["y1"])])
        x2 = np.array([float(row["x2"]), float(row["y2"])])
        ds = np.array([float(row["dSx"]), float(row["dSy"])])
        length = np.linalg.norm(ds)
        coeff = .5*ds + length*.5*(x1+x2-xi-xj)/(np.linalg.norm(xj-xi)+1e-30)
        forces[category][i] += (float(a["P"])-float(b["P"]))*coeff
        counts[category] += 1
        lengths[category] += length
    ids = list(selected)
    mass = np.array([float(by_id[i]["mass"]) for i in ids])
    weights = mass / mass.sum()
    observed = np.array([by_id[i]["nongrav_acceleration"] for i in ids])
    original = np.array([by_id[i]["face_pressure_acceleration"] for i in ids])
    arrays = {key: np.array([forces[key][i]/mass[k] for k, i in enumerate(ids)])
              for key in categories}
    total = sum(arrays.values())
    max_reconstruction_error = float(np.max(np.abs(total-original)))
    assert max_reconstruction_error < 1e-15
    obs_norm2 = float(np.sum(weights*np.sum(observed*observed, axis=1)))
    out = {}
    for key, acc in arrays.items():
        norm2 = float(np.sum(weights*np.sum(acc*acc, axis=1)))
        dot = float(np.sum(weights*np.sum(acc*observed, axis=1)))
        out[key] = {
            "faces": counts[key],
            "total_face_length": lengths[key],
            "acceleration_rms": norm2**.5,
            "projection_on_observed_fraction": dot/obs_norm2,
            "vector_cosine_to_observed": dot/(norm2*obs_norm2)**.5 if norm2 else None,
        }
    return {
        "sample_cells": len(ids),
        "observed_nongrav_acceleration_rms": obs_norm2**.5,
        "reconstructed_total_acceleration_rms": float(np.sum(weights*np.sum(total*total, axis=1)))**.5,
        "max_reconstruction_error": max_reconstruction_error,
        "max_same_shell_initial_radial_difference": max_same_dr,
        "min_cross_shell_initial_radial_difference": min_cross_dr,
        "categories": out,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    receipt = json.loads((DOCS / "KEPLER_Q2_BOUNDED_FACEFORCE_PAIR_T6501_2026-10-27.json").read_text())
    root = Path(receipt["staged_pair"])
    selected = json.loads((root / "stage_manifest.json").read_text())["selected_ids"]
    result = {
        "source_receipt": str(DOCS / "KEPLER_Q2_BOUNDED_FACEFORCE_PAIR_T6501_2026-10-27.json"),
        "cases": {label: evaluate(label, root, receipt, selected) for label in CASES},
        "scope": "64 sampled inner disk cells at t=6.5, one W0 pressure-force call; initial ring grouping, not an irreversible heat or full-disk audit.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
