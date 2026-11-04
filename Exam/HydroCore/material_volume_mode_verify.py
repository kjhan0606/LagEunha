"""Independent bounded verification of surprising frozen acoustic stiffness.

Reassemble from native face records. Use 32 power iterations per operator,
inside the previous stage's combined 256-matvec/15-second limits. No coefficient
solve or fluid evolution. Also check Fourier edge integration on exact squares.
"""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import time

import numpy as np
from scipy import sparse

from rotating_shear_assessment import read_native


def fourier_square_check():
    n = np.array([[1., 0], [-1, 0], [0, 1], [0, -1]])
    t = np.column_stack([-n[:, 1], n[:, 0]])
    side, center = 0.13, np.array([0.27, 0.39])
    face_center = center + 0.5*side*n
    rows = []
    for k in (np.array([1., 2.]), np.array([25., 0]), np.array([100., 100.])):
        face_integral = side*np.exp(1j*(face_center@k))*np.sinc((t@k)*side/(2*np.pi))
        measured = np.sum((n@k)*face_integral)/(1j*np.dot(k, k)*side**2)
        exact = np.exp(1j*np.dot(center, k))*np.prod(np.sinc(k*side/(2*np.pi)))
        error = abs(measured-exact)
        assert error < 1e-13
        rows.append({"wave": k.tolist(), "absolute_error": float(error)})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir()
    native_input = args.root / "material_volume_native_cost_v3/material_volume_native_input.bin"
    with native_input.open("rb") as f:
        magic, n, nf = struct.unpack("<QQQ", f.read(24))
        assert magic == 0x4D56434F53543131
        cells = np.fromfile(f, "<f8", count=n*6).reshape(n, 6)
        faces = np.fromfile(f, dtype=[("i", "<i4"), ("j", "<i4"), ("r", "<f8", (2,)),
                                     ("c", "<f8", (2,)), ("A", "<f8")], count=nf)
        assert not f.read()
    psi = np.load(args.root / "material_volume_reference_v1/psi_reference.npy")
    i, j, r, centroid, area = (faces[k] for k in ["i", "j", "r", "c", "A"])
    d = np.linalg.norm(r, axis=1)
    M = np.zeros((n, 2, 2))
    np.add.at(M, i, (area/d**3)[:, None, None]*np.einsum("ka,kb->kab", r, r))
    inverse = np.linalg.inv(M)
    ell = (area/d**3)[:, None]*np.einsum("kab,kb->ka", inverse[i], r)
    tensors = np.array([[[1., 0], [0, -1]], [[0., 1], [1, 0]]])/np.sqrt(2)
    rv = np.einsum("kb,bad,kd->ka", psi[i], tensors, ell)
    bj = area[:, None]*(r-centroid)/d[:, None]

    def operator(vector):
        rows, cols, data = [], [], []
        for a in (0, 1):
            rows.extend([i, i]); cols.extend([2*j+a, 2*i+a]); data.extend([vector[:, a], -vector[:, a]])
        return sparse.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(n, 2*n))

    B, R = operator(bj), operator(rv)
    conductance = sparse.csr_matrix((area/(2*d), (i, j)), shape=(n, n))
    conductance = (conductance+conductance.T)/2
    L = sparse.diags(np.asarray(conductance.sum(axis=1)).ravel())-conductance
    C = (B+L@R).tocsr()
    state, snapshot_time = read_native(args.root / "rotating_shear_suite_v1/gresho_128x128/khout.000543.dat")
    assert snapshot_time == .5 and np.array_equal(state["x"], cells[:, 0])
    gamma = float(np.float32(5/3))
    sqrtH = np.sqrt(gamma*cells[:, 3]/cells[:, 2])
    inverse_sqrt_mass = np.repeat(1/np.sqrt(state["m"]), 2)
    previous = json.loads((args.root / "material_volume_modes_v2/results/native_probes.json").read_text())["positive_acoustic_block"]
    results = {}
    for label, X in (("B", B), ("C", C)):
        mark = time.monotonic()
        count = 0

        def action(v):
            nonlocal count
            count += 1
            assert count+previous[label]["matvecs"] <= 256
            assert time.monotonic()-mark+previous[label]["seconds"] < 15
            # Independent weighted rectangular factor gives A = T T^T.
            factor_action = inverse_sqrt_mass*(X.T@(sqrtH*v))
            return sqrtH*(X@(inverse_sqrt_mass*factor_action))

        v = np.cos(0.37*np.arange(n))
        v /= np.linalg.norm(v)
        for _ in range(32):
            v = action(v)
            v /= np.linalg.norm(v)
        Av = action(v)
        value = float(v@Av)
        residual = float(np.linalg.norm(Av-value*v))
        assert value > 0 and value <= previous[label]["squared_frequency"]*(1+1e-10)
        shape = v*v
        pressure_participation = 1/np.sum(shape*shape)
        velocity = inverse_sqrt_mass*(X.T@(sqrtH*v))/np.sqrt(value)
        kinetic = np.sum(velocity.reshape(n, 2)**2, axis=1)
        kinetic /= kinetic.sum()
        top = np.argsort(shape)[-8:][::-1]
        minimum_distance = np.full(n, np.inf)
        np.minimum.at(minimum_distance, i, d)
        max_rnorm = np.zeros(n)
        np.maximum.at(max_rnorm, i, np.linalg.norm(rv, axis=1))
        record = {"power_iterations": 32, "additional_matvecs": count,
                  "combined_stage_matvecs": count+previous[label]["matvecs"],
                  "combined_diagnostic_seconds": time.monotonic()-mark+previous[label]["seconds"],
                  "squared_frequency_Rayleigh": value,
                  "relative_difference_from_Ritz": value/previous[label]["squared_frequency"]-1,
                  "relative_power_residual": residual/value,
                  "pressure_mode_participation_cells": float(pressure_participation),
                  "kinetic_mode_participation_cells": float(1/np.sum(kinetic*kinetic)),
                  "top_pressure_cells": [{"index": int(k), "id": int(state["id"][k]),
                                          "x": float(cells[k, 0]), "y": float(cells[k, 1]),
                                          "pressure_mode_norm_fraction": float(shape[k]),
                                          "volume_over_mean": float(n*cells[k, 2]),
                                          "nearest_face_neighbor_distance_over_h": float(np.sqrt(n)*minimum_distance[k]),
                                          "largest_R_face_coefficient_over_h": float(np.sqrt(n)*max_rnorm[k])} for k in top]}
        results[label] = record
    assert abs(results["C"]["relative_difference_from_Ritz"]) < 1e-6
    result = {"scope": __doc__, "operator_results": results,
              "Fourier_square_checks": fourier_square_check(),
              "accepted_steps": 0, "new_coefficient_solves": 0,
              "full_background_stability_qualified": False,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
