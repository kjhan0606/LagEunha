"""Finite algebra proof on exported native geometry; never evolves a fluid.

The psi coefficient construction is a reference dense eigensolve, not a
production implementation or a performance measurement. Production PDE,
pressure force, weights, and time integration remain untouched.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess

import numpy as np


def basis(dim):
    matrices = []
    for k in range(1, dim):
        s = np.zeros((dim, dim))
        for a in range(k):
            s[a, a] = 1 / np.sqrt(k * (k + 1))
        s[k, k] = -k / np.sqrt(k * (k + 1))
        matrices.append(s)
    for a in range(dim):
        for b in range(a + 1, dim):
            s = np.zeros((dim, dim))
            s[a, b] = s[b, a] = 1 / np.sqrt(2)
            matrices.append(s)
    return np.array(matrices).reshape((-1, dim, dim))


def load_geometry(executable, directory, dim, side, fraction=0):
    directory.mkdir()
    command = [str(executable), str(dim), str(side), str(fraction)]
    (directory / "command.json").write_text(json.dumps(command) + "\n")
    with (directory / "geometry.json").open("x") as out, (directory / "log").open("x") as log:
        p = subprocess.Popen(command, stdout=out, stderr=log)
        (directory / "process.json").write_text(json.dumps({"pid": p.pid, "accepted_steps": 0}) + "\n")
        code = p.wait()
    assert code == 0, (directory, code)
    return json.loads((directory / "geometry.json").read_text())


def assemble(g):
    n, dim = g["cells"], g["dimension"]
    V, x, W = np.array(g["volume"]), np.array(g["positions"]), np.array(g["weights"])
    assert np.all(V > 0) and abs(V.sum() - 1) < 1e-10
    inverse = np.array(g["native_inverse_fit"])
    B = np.zeros((n, n * dim))
    L = np.zeros((n, n))
    H = np.zeros((n, dim, dim, n * dim))
    K = np.zeros((n, dim, dim))
    moment = np.zeros_like(K)
    cut = np.zeros(n, dtype=bool)
    neighbors = [set() for _ in range(n)]
    for f in g["faces"]:
        i, j = f["i"], f["j"]
        r, bi, bj = np.array(f["r"]), np.array(f["bi"]), np.array(f["bj"])
        distance = np.linalg.norm(r)
        B[i, dim * j:dim * (j + 1)] += bj
        B[i, dim * i:dim * (i + 1)] -= bj
        coeff = f["area"] / (2 * distance)
        L[i, i] += coeff
        L[i, j] -= coeff
        K[i] += np.outer(bj, r) / V[i]
        moment[i] += np.outer(bi, r) / V[i]
        ell = f["area"] / distance**(dim + 1) * (inverse[i] @ r)
        for a in range(dim):
            H[i, a, :, dim * j + a] += ell
            H[i, a, :, dim * i + a] -= ell
        neighbors[i].add(j)
        cut[i] |= np.max(abs(r - (x[j] - x[i]))) > 1e-9
    assert np.max(abs(moment - np.eye(dim))) < 1e-10
    assert np.max(abs(L - L.T)) < 1e-10
    S = basis(dim)
    rhs = np.einsum("iab,sab->is", V[:, None, None] * (np.eye(dim) - K), S)
    compatibility = np.sum(rhs, axis=0)
    threshold = 4096 * np.finfo(float).eps * max(1.0, float(np.sum(abs(rhs))))
    assert np.max(abs(compatibility), initial=0) < threshold, ("incompatible RHS", compatibility)
    eigenvalues, vectors = np.linalg.eigh(L)
    assert abs(eigenvalues[0]) < 1e-10 * eigenvalues[-1]
    assert eigenvalues[1] > 1e-10 * eigenvalues[-1]
    psi = vectors[:, 1:] @ ((vectors[:, 1:].T @ rhs) / eigenvalues[1:, None])
    psi -= psi.mean(axis=0)
    solve_error = np.linalg.norm(L @ psi - rhs) / max(np.linalg.norm(rhs), 1e-30) if len(S) else 0
    assert solve_error < 1e-10
    strains = np.einsum("sab,iabk->isk", S, H)
    divergence = np.einsum("iaak->ik", H)

    def make_R(weights):
        centered = weights - weights.mean()
        return np.einsum("is,isk->ik", psi, strains) + (2 / dim) * centered[:, None] * divergence

    R = make_R(W)
    C = B + L @ R
    gauge_error = np.max(abs(make_R(W + 0.37) - R))
    assert gauge_error < 1e-10
    core = np.array([not cut[i] and not any(cut[j] for j in neighbors[i]) for i in range(n)])
    return {"V": V, "x": x, "W": W, "B": B, "L": L, "H": H, "K": K,
            "S": S, "psi": psi, "R": R, "C": C, "core": core,
            "solve_error": float(solve_error), "compatibility_max": float(np.max(abs(compatibility), initial=0)),
            "weight_common_shift_R_error": float(gauge_error)}


def algebra(g, a):
    dim, n = g["dimension"], g["cells"]
    V, x, C, B, L = (a[key] for key in ("V", "x", "C", "B", "L"))
    pressure = 1 + 0.3 * np.cos(0.9 * np.arange(n))
    velocity = np.sin(0.31 * np.arange(n * dim))
    force = C.T @ pressure
    base_force = np.zeros((n, dim))
    for f in g["faces"]:
        base_force[f["i"]] += (pressure[f["i"]] - pressure[f["j"]]) * np.array(f["bi"])
    base_adjoint_error = np.max(abs(base_force.ravel() - B.T @ pressure))
    constant_pressure_error = np.max(abs(C.T @ np.ones(n)))
    translation_error = max(np.max(abs(C @ np.tile(np.eye(dim)[d], n))) for d in range(dim))
    momentum_error = np.max(abs(force.reshape(n, dim).sum(axis=0)))
    power_error = abs(velocity @ force - pressure @ (C @ velocity))
    affine_error = 0
    cut_rotation = []
    strains = [np.eye(dim), np.array([[0.23 * np.sin(0.7 * (i + 1) * (j + 2))
                                     for j in range(dim)] for i in range(dim)])]
    if dim > 1:
        rotation = np.zeros((dim, dim))
        rotation[0, 1], rotation[1, 0] = -1, 1
        strains.append(rotation)
    for gradient in strains:
        local_wdot = a["psi"] @ np.einsum("sab,ab->s", a["S"], gradient)
        local_wdot += (2 / dim) * (a["W"] - a["W"].mean()) * np.trace(gradient)
        patch = V * np.einsum("iab,ab->i", a["K"], gradient) + L @ local_wdot
        affine_error = max(affine_error, float(np.max(abs(patch / V - np.trace(gradient)))))
        owner_velocity = (x @ gradient.T).ravel()
        cut_rate = patch - C @ owner_velocity
        cut_rotation.append({"gradient": gradient.tolist(),
                             "owner_pressure_power": float(owner_velocity @ force),
                             "cut_pressure_power": float(pressure @ cut_rate),
                             "thermal_pressure_power": float(-pressure @ patch),
                             "closed_power_residual": float(owner_velocity @ force + pressure @ cut_rate - pressure @ patch)})
    gradient_p = 0.31 * np.arange(1, dim + 1)
    linear_force = (C.T @ (1 + x @ gradient_p)).reshape(n, dim)
    patch_error = float(np.max(abs(linear_force[a["core"]] / V[a["core"], None] + gradient_p))) if a["core"].any() else None
    assert max(base_adjoint_error, constant_pressure_error, translation_error,
               momentum_error, power_error, affine_error) < 1e-9
    if patch_error is not None:
        assert patch_error < 1e-9
    if g["side"] >= 6:
        assert a["core"].any(), "No interior pressure-patch samples"
    return {"dimension": dim, "side": g["side"], "cells": n,
            "rhs_compatibility_max": a["compatibility_max"], "solve_relative_residual": a["solve_error"],
            "constant_pressure_force_error": float(constant_pressure_error),
            "translation_volume_rate_error": float(translation_error), "momentum_error": float(momentum_error),
            "work_adjoint_error": float(power_error), "base_native_pressure_adjoint_error": float(base_adjoint_error),
            "affine_transport_patch_error": affine_error, "interior_pressure_samples": int(a["core"].sum()),
            "interior_linear_pressure_error": patch_error,
            "weight_common_shift_R_error": a["weight_common_shift_R_error"],
            "one_dimensional_W0_R_max": float(np.max(abs(a["R"]))) if dim == 1 else None,
            "psi_gauge": "arithmetic mean zero", "periodic_cut_power_budgets": cut_rotation,
            "necessary_algebra_checks_pass": True}


def nonaffine(g, a):
    n, dim = g["cells"], g["dimension"]
    v = np.zeros((n, dim))
    v[:, 0] = np.sin(2 * np.pi * a["x"][:, 1])
    v = v.ravel()
    G = np.einsum("iabk,k->iab", a["H"], v)
    div = np.trace(G, axis1=1, axis2=2)
    scalar_strains = np.einsum("sab,iab->is", a["S"], G)
    commutator = np.zeros(n)
    for k in range(len(a["S"])):
        p, s = a["psi"][:, k], scalar_strains[:, k]
        commutator += a["L"] @ (p * s) - s * (a["L"] @ p)
    w = a["W"] - a["W"].mean()
    commutator += (2 / dim) * (a["L"] @ (w * div) - div * (a["L"] @ w))
    curvature = a["B"] @ v - a["V"] * np.einsum("iab,iab->i", a["K"], G)
    residual = a["C"] @ v - a["V"] * div
    decomposition_error = np.max(abs(residual - curvature - commutator))
    assert decomposition_error < 1e-10
    return {"dimension": dim, "side": g["side"], "cells": n,
            "fixed_physical_wavelength": 1.0,
            "base_divergence_gap_L1": float(np.sum(abs(a["B"] @ v - a["V"] * div))),
            "candidate_divergence_gap_L1": float(np.sum(abs(residual))),
            "candidate_physical_divergence_L1": float(np.sum(abs(a["C"] @ v))),
            "commutator_L1": float(np.sum(abs(commutator))),
            "curvature_remainder_L1": float(np.sum(abs(curvature))),
            "decomposition_max_error": float(decomposition_error)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir()
    registration = {"scope": "Necessary algebra and nonaffine consistency only, actual project geometry, zero PDE steps. Dense reference solves are not a performance result.",
                    "algebra_cases": [[1, 16], [1, 32], [2, 3], [2, 6], [3, 3], [3, 6]],
                    "nonaffine_cases": [[2, 6], [2, 12], [2, 24]],
                    "gauge": "arithmetic mean zero for psi and W",
                    "gate": "First substantive failure stops candidate before weighted geometry/cost/PDE. Nonaffine gap and physical-divergence norms must decrease with refinement."}
    (args.output / "registration.json").write_text(json.dumps(registration, indent=2) + "\n")
    algebra_rows, nonaffine_rows = [], []
    try:
        for dim, side in registration["algebra_cases"]:
            g = load_geometry(args.executable, args.output / f"algebra_D{dim}_N{side}", dim, side)
            row = algebra(g, assemble(g))
            algebra_rows.append(row)
            print("algebra", dim, side, row["affine_transport_patch_error"], flush=True)
        (args.output / "algebra.json").write_text(json.dumps(algebra_rows, indent=2) + "\n")
        for dim, side in registration["nonaffine_cases"]:
            g = load_geometry(args.executable, args.output / f"nonaffine_D{dim}_N{side}", dim, side)
            row = nonaffine(g, assemble(g))
            nonaffine_rows.append(row)
            print("nonaffine", json.dumps(row), flush=True)
        (args.output / "nonaffine.json").write_text(json.dumps(nonaffine_rows, indent=2) + "\n")
        for key in ["candidate_divergence_gap_L1", "candidate_physical_divergence_L1"]:
            assert all(b[key] < a[key] for a, b in zip(nonaffine_rows, nonaffine_rows[1:])), key
        status = {"status": "necessary algebra/nonaffine gates pass; weighted FD and native cost still required", "accepted_steps": 0}
    except Exception as error:
        status = {"status": "candidate stopped at finite gate", "error": repr(error),
                  "completed_algebra_cases": algebra_rows, "completed_nonaffine_cases": nonaffine_rows,
                  "accepted_steps": 0}
        (args.output / "status.json").write_text(json.dumps(status, indent=2) + "\n")
        raise
    (args.output / "status.json").write_text(json.dumps(status, indent=2) + "\n")
