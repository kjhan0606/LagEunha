"""Minimum uniform-pressure stiffness control using analytical inequalities.

No new eigensolve, coefficient solve, or PDE step. Reuse the independently
verified C Rayleigh value. B uses induced rectangular-matrix norm bounds,
not an assumed Ritz upper bound. Floating arithmetic is not interval certified.
"""
from pathlib import Path
import argparse
import hashlib
import json
import struct

import numpy as np
from scipy import sparse

from rotating_shear_assessment import read_native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir()
    with (args.root/"material_volume_native_cost_v3/material_volume_native_input.bin").open("rb") as f:
        magic, n, nf = struct.unpack("<QQQ", f.read(24))
        assert magic == 0x4D56434F53543131
        cells = np.fromfile(f, "<f8", count=n*6).reshape(n, 6)
        faces = np.fromfile(f, dtype=[("i", "<i4"), ("j", "<i4"), ("r", "<f8", (2,)),
                                     ("c", "<f8", (2,)), ("A", "<f8")], count=nf)
        assert not f.read()
    i, j, r, centroid, area = (faces[k] for k in ["i", "j", "r", "c", "A"])
    d = np.linalg.norm(r, axis=1)
    psi = np.load(args.root/"material_volume_reference_v1/psi_reference.npy")
    Mfit = np.zeros((n, 2, 2))
    np.add.at(Mfit, i, (area/d**3)[:, None, None]*np.einsum("ka,kb->kab", r, r))
    ell = (area/d**3)[:, None]*np.einsum("kab,kb->ka", np.linalg.inv(Mfit)[i], r)
    S = np.array([[[1., 0], [0, -1]], [[0., 1], [1, 0]]])/np.sqrt(2)
    Rface = np.einsum("kb,bad,kd->ka", psi[i], S, ell)

    def matrix(vector):
        row = np.repeat(i, 2)
        neighbor = 2*np.repeat(j, 2)+np.tile([0, 1], nf)
        own = 2*np.repeat(i, 2)+np.tile([0, 1], nf)
        return sparse.csr_matrix((np.r_[vector.ravel(), -vector.ravel()],
                                  (np.tile(row, 2), np.r_[neighbor, own])), shape=(n, 2*n))

    B = matrix(area[:, None]*(r-centroid)/d[:, None])
    adjacency = sparse.csr_matrix((area/(2*d), (i, j)), shape=(n, n))
    adjacency = (adjacency+adjacency.T)/2
    L = sparse.diags(np.asarray(adjacency.sum(axis=1)).ravel())-adjacency
    C = B+L@matrix(Rface)
    state, t = read_native(args.root/"rotating_shear_suite_v1/gresho_128x128/khout.000543.dat")
    assert t == .5 and np.array_equal(state["V"], cells[:, 2]) and np.all(state["m"] == 1/n)
    gamma = float(np.float32(5/3))
    P0 = 1/gamma
    V = cells[:, 2]
    U0 = P0*V/(gamma-1)
    assert np.all(U0 > 0) and np.all(V > 0)
    factor = sparse.diags(np.sqrt(gamma*P0/V)) @ B @ sparse.diags(np.repeat(1/np.sqrt(state["m"]), 2))
    absolute = abs(factor)
    rows = np.asarray(absolute.sum(axis=1)).ravel()
    columns = np.asarray(absolute.sum(axis=0)).ravel()
    norm_product = rows.max()*columns.max()
    # |T T^T| <= |T| |T|^T entrywise; row sums supply a second upper bound.
    row_sum_upper = np.max(absolute@columns)
    pressure = cells[:, 3]
    prior = json.loads((args.root/"material_volume_modes_v2/verification/result.json").read_text())
    rayleigh = prior["operator_results"]["C"]["squared_frequency_Rayleigh"]
    # A(P) and T(P) have the same nonzero spectrum as the velocity block.
    # K(P) <= (Pmax/P0) K(P0), so lambda_C(P0) >= P0/Pmax * Rayleigh_C(P).
    padding = 1e-9
    C_lower = (1-padding)*P0/pressure.max()*rayleigh
    B_upper = (1+padding)*min(norm_product, row_sum_upper)
    ratio_lower = np.sqrt(C_lower/B_upper)
    constant_force = {label: float(np.max(abs(X.T@np.ones(n)))) for label, X in [("B", B), ("C", C)]}
    relative_force = {label: constant_force[label]/float(np.max(np.asarray(abs(X).sum(axis=0)))) for label, X in [("B", B), ("C", C)]}
    assert max(relative_force.values()) < 1e-12
    assert C_lower > B_upper
    result = {"scope": __doc__, "cells": n, "diagnostic_equilibrium": {"velocity": 0, "P0": P0,
              "U0_min": float(U0.min()), "fixed_original_masses": True, "weight": 0,
              "uniform_pressure_force_max": constant_force, "relative_constant_pressure_force_defect": relative_force},
              "original_pressure_range": [float(pressure.min()), float(pressure.max())],
              "B_uniform_lambda_upper_norm1_times_norminf": float((1+padding)*norm_product),
              "B_uniform_lambda_upper_absolute_product_row_sum": float((1+padding)*row_sum_upper),
              "B_uniform_lambda_upper_used": float(B_upper),
              "C_original_pressure_Rayleigh_reused": rayleigh,
              "C_uniform_lambda_lower_from_pressure_comparison": float(C_lower),
              "uniform_frequency_ratio_lower_estimate": float(ratio_lower),
              "arithmetic_padding_fraction": padding, "interval_arithmetic_certified": False,
              "B_Ritz_value_used_as_upper_bound": False,
              "new_spectral_matvecs": 0, "new_coefficient_solves": 0, "accepted_steps": 0,
              "new_native_spectrum_needed": False,
              "interpretation": "Large stiffness persists in a uniform-pressure rest equilibrium. C^T1 identically zero removes the background geometric force term in the inviscid pressure linearization. This is a pressure-subsystem result, not full nonlinear/AV/KDK or GIZMO runtime validation.",
              "current_psi_times_strain_candidate_retired": True,
              "retirement_scope": "Current local psi_i times fitted strain construction only; not every C=B+LR model.",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
