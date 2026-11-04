"""Finite eligibility screen of a midpoint-offset Gauss/adjoint operator.

Analysis of existing original-native geometry only. No fluid evolution, new
mesh, global weight/coefficient solve, or production change. The pressure gate
precedes Bloch/native stiffness/cost work. Pressure representation is fixed
before measurements and is never selected from the better-looking result.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json

import numpy as np
from scipy import sparse


def load_geometry(path):
    g = json.loads(path.read_text())
    return {"dimension": g["dimension"], "side": g["side"],
            "V": np.array(g["volume"]), "x": np.array(g["positions"]),
            "inverse": np.array(g["native_inverse_fit"]),
            "i": np.array([f["i"] for f in g["faces"]]),
            "j": np.array([f["j"] for f in g["faces"]]),
            "r": np.array([f["r"] for f in g["faces"]]),
            "c": np.array([f["c"] for f in g["faces"]]),
            "area": np.array([f["area"] for f in g["faces"]]),
            "weights": np.array(g["weights"]),
            "geometry_path": str(path), "geometry_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def tile_geometry(g, repeats):
    d, primitive = g["dimension"], len(g["V"])
    # Tile and native index conventions both place x first in the linear index.
    tiles = np.array(list(itertools.product(range(repeats), repeat=d)))[:, ::-1]
    shifts = g["r"]-(g["x"][g["j"]]-g["x"][g["i"]])
    assert np.max(abs(shifts-np.rint(shifts))) < 1e-12
    shifts = np.rint(shifts).astype(int)
    scale = 1/repeats
    arrays = {k: [] for k in ["V", "x", "inverse", "i", "j", "r", "c", "area", "weights"]}
    for tile_index, tile in enumerate(tiles):
        neighbor_tiles = (tile+shifts) % repeats
        linear = sum(neighbor_tiles[:, a]*repeats**a for a in range(d))
        arrays["i"].append(tile_index*primitive+g["i"])
        arrays["j"].append(linear*primitive+g["j"])
        arrays["V"].append(g["V"]*scale**d)
        arrays["x"].append((g["x"]+tile)*scale)
        arrays["inverse"].append(g["inverse"])
        arrays["r"].append(g["r"]*scale)
        arrays["c"].append(g["c"]*scale)
        arrays["area"].append(g["area"]*scale**(d-1))
        arrays["weights"].append(g["weights"]*scale**2)
    out = {k: np.concatenate(v) for k, v in arrays.items()}
    center_tile = np.full(d, repeats//2)
    center_linear = sum(center_tile[a]*repeats**a for a in range(d))
    out.update(dimension=d, side=g["side"]*repeats, primitive_cells=primitive,
               central_pattern=np.arange(center_linear*primitive, (center_linear+1)*primitive),
               repeats=repeats, geometry_path=g["geometry_path"], geometry_sha256=g["geometry_sha256"])
    return out


def operators(g):
    d, n = g["dimension"], len(g["V"])
    i, j, r, c, A = (g[k] for k in ["i", "j", "r", "c", "area"])
    distance = np.linalg.norm(r, axis=1)
    normal = r/distance[:, None]
    offset = c-r/2
    ell = (A/distance**(d+1))[:, None]*np.einsum("kab,kb->ka", g["inverse"][i], r)
    H = [sparse.csr_matrix((np.r_[ell[:, b], -ell[:, b]],
                            (np.tile(i, 2), np.r_[j, i])), shape=(n, n)) for b in range(d)]
    B_columns, D_columns, extra_columns = [], [], []
    for a in range(d):
        bj = A*(r[:, a]-c[:, a])/distance
        B_columns.append(sparse.csr_matrix((np.r_[bj, -bj], (np.tile(i, 2), np.r_[j, i])), shape=(n, n)))
        flux = A*normal[:, a]/2
        midpoint = sparse.csr_matrix((np.tile(flux, 2), (np.tile(i, 2), np.r_[i, j])), shape=(n, n))
        extra = sparse.csr_matrix((n, n))
        for b in range(d):
            q = flux*offset[:, b]
            average_gradient = sparse.csr_matrix((np.tile(q, 2), (np.tile(i, 2), np.r_[i, j])), shape=(n, n))
            midpoint = midpoint+average_gradient@H[b]
            qold = A*normal[:, a]*r[:, b]/4
            gradient_difference = sparse.csr_matrix((np.r_[qold, -qold], (np.tile(i, 2), np.r_[i, j])), shape=(n, n))
            extra = extra+gradient_difference@H[b]
        D_columns.append(midpoint)
        extra_columns.append(extra)
    # hstack groups components; permute columns to cell-major native ordering.
    permutation = np.arange(n*d).reshape(d, n).T.ravel()
    B = sparse.hstack(B_columns, format="csr")[:, permutation]
    D = sparse.hstack(D_columns, format="csr")[:, permutation]
    extra = sparse.hstack(extra_columns, format="csr")[:, permutation]
    cut = np.zeros(n, dtype=bool)
    shift = r-(g["x"][j]-g["x"][i])
    np.logical_or.at(cut, i, np.max(abs(shift), axis=1) > 1e-10)
    unsafe = cut.copy()
    np.logical_or.at(unsafe, i, cut[j])
    center = np.zeros((n, d))
    np.add.at(center, i, (A*np.einsum("ka,ka->k", c, normal))[:, None]*c)
    center /= ((d+1)*g["V"])[:, None]
    return {"B": B, "D": D, "extra": extra, "H": H, "ell": ell,
            "normal": normal, "offset": offset, "shift": shift,
            "core": ~unsafe, "centroid": g["x"]+center}


def gradients(a, velocity):
    n, d = velocity.shape
    return np.stack([np.column_stack([a["H"][b]@velocity[:, k] for b in range(d)])
                     for k in range(d)], axis=1).reshape(n, d, d)


def direct_face_divergence(g, a, velocity, old=False):
    i, j = g["i"], g["j"]
    G = gradients(a, velocity)
    vf = (velocity[i]+velocity[j])/2+np.einsum("kab,kb->ka", (G[i]+G[j])/2, a["offset"])
    if old:
        vf += np.einsum("kab,kb->ka", (G[i]-G[j])/4, g["r"])
    return np.bincount(i, weights=g["area"]*np.einsum("ka,ka->k", a["normal"], vf), minlength=len(g["V"]))


def direct_transpose(g, a, pressure):
    d, n = g["dimension"], len(g["V"])
    i, j = g["i"], g["j"]
    q = (pressure[i]*g["area"]/2)[:, None]*a["normal"]
    force = np.zeros((n, d))
    np.add.at(force, i, q); np.add.at(force, j, q)
    moment = np.zeros((n, d, d))
    term = np.einsum("ka,kb->kab", q, a["offset"])
    np.add.at(moment, i, term); np.add.at(moment, j, term)
    gradient_force = np.einsum("kab,kb->ka", moment[i], a["ell"])
    np.add.at(force, j, gradient_force); np.add.at(force, i, -gradient_force)
    return force


def algebra(g, a):
    d, n = g["dimension"], len(g["V"])
    velocity = np.sin(.37*np.arange(n*d)).reshape(n, d)
    pressure = 2+.3*np.cos(.71*np.arange(n))
    D, B = a["D"], a["B"]
    direct = direct_face_divergence(g, a, velocity)
    force = direct_transpose(g, a, pressure)
    delta = direct_face_divergence(g, a, velocity, old=True)-direct
    translation = max(np.max(abs(D@np.tile(np.eye(d)[b], n))) for b in range(d))
    checks = {"matrix_vs_direct_face_max": float(np.max(abs(D@velocity.ravel()-direct))),
              "transpose_vs_direct_face_max": float(np.max(abs((D.T@pressure).reshape(n, d)-force))),
              "historical_removed_term_max_defect": float(np.max(abs(delta-a["extra"]@velocity.ravel()))),
              "translation_max": float(translation),
              "constant_pressure_force_max": float(np.max(abs(D.T@np.ones(n)))),
              "force_work_adjoint_error": float(abs(np.sum(velocity*force)-pressure@direct)),
              "weight_rhs_compatibility_sum": float(np.sum((D-B)@velocity.ravel()))}
    affine = []
    matrices = [("dilation", np.eye(d)), ("general", np.sin(.41*np.arange(d*d)).reshape(d, d))]
    if d > 1:
        rotation = np.zeros((d, d));rotation[0, 1]=-1;rotation[1, 0]=1
        shear = np.zeros((d, d));shear[0, 1]=1
        matrices.extend([("rotation", rotation), ("shear", shear)])
    for name, G in matrices:
        v = g["x"]@G.T
        owner_G = gradients(a, v)
        correction_G = np.zeros_like(owner_G)
        np.add.at(correction_G, g["i"], np.einsum("ka,kb->kab", a["shift"]@G.T, a["ell"]))
        gradient_error = np.max(abs(owner_G+correction_G-G))
        vf = (v[g["i"]]+v[g["j"]]+a["shift"]@G.T)/2+g["c"]@G.T-g["r"]@G.T/2
        exact = np.bincount(g["i"], weights=g["area"]*np.einsum("ka,ka->k", a["normal"], vf), minlength=n)
        cut_vf = (a["shift"]@G.T)/2+np.einsum("kab,kb->ka", (correction_G[g["i"]]+correction_G[g["j"]])/2, a["offset"])
        cut_rate = np.bincount(g["i"], weights=g["area"]*np.einsum("ka,ka->k", a["normal"], cut_vf), minlength=n)
        affine.append({"name": name, "relative_volume_rate_max_error": float(np.max(abs(exact/g["V"]-np.trace(G)))),
                       "continued_gradient_max_error": float(gradient_error),
                       "owner_plus_cut_rate_max_error": float(np.max(abs(D@v.ravel()+cut_rate-exact))),
                       "pressure_cut_work_ledger_error": float(abs(v.ravel()@(D.T@pressure)+pressure@cut_rate-pressure@exact))})
    scale = max(1., float(np.max(abs(D.data), initial=0)))
    assert max(abs(x) for x in checks.values()) < 1e-10*scale, checks
    assert max(x[k] for x in affine for k in x if k != "name") < 1e-10*scale, affine
    reduction = float(np.max(abs((D-B).data), initial=0)) if d == 1 and np.all(g["weights"] == 0) else None
    if reduction is not None:
        assert reduction < 1e-12
    return {"dimension": d, "cells": n, "side": g["side"], "checks": checks,
            "affine_and_periodic_cut": affine, "one_dimensional_W0_D_minus_B_max": reduction,
            "historical_old_formula_minus_B_max": float(np.max(abs((D+a["extra"]-B).data), initial=0)),
            "algebra_pass": True}


def pressure_patch(g, a, selection):
    d = g["dimension"]
    assert len(selection) and a["core"][selection].all(), "Pressure samples must avoid both reconstruction layers of the seam."
    V = g["V"][selection]
    rows = []
    for axis in range(d):
        gradient = np.eye(d)[axis]
        for representation, position in [("generator_sample", g["x"]), ("exact_cell_average", a["centroid"])]:
            pressure = 2+position@gradient
            for label, X in [("B", a["B"]), ("D", a["D"])]:
                force = (X.T@pressure).reshape(-1, d)[selection]
                error = np.linalg.norm(force/V[:, None]+gradient, axis=1)
                rows.append({"axis": axis, "representation": representation, "operator": label,
                             "force_density_error_L1": float(V@error/V.sum()),
                             "force_density_error_max": float(error.max())})
    return {"dimension": d, "side": g["side"], "cells": len(g["V"]),
            "pressure_samples": len(selection), "patches": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir()
    registration = {"scope": __doc__, "algebra_cases": [[1, 16], [1, 32], [2, 3], [2, 6], [3, 3], [3, 6]],
                    "pressure_pattern_repeats": [3, 6], "pressure_samples": "One complete central primitive pattern, guarded through both reconstruction layers.",
                    "primary_pressure_representation": "exact cell average consistent with cell thermodynamic P=(gamma-1)U/V",
                    "secondary_pressure_representation": "generator sample, diagnosed separately without changing the physical gate",
                    "pressure_stop_rule": "For any dimension/axis, a primary D force-density error above 1e-9 retaining >=99% of its value under factor-two pattern refinement rejects the operator stage. Record inherited B errors separately. Algebra failure also stops. No tolerance or representation selection after measurement.",
                    "later_phases_only_if_pass": "Small Bloch and N128 known-wave/uniform-P acoustic screen; fixed 256 matvec/15 s limits and omega_D/omega_B<=2 engineering ceiling. No weight solve or PDE in this bundle.",
                    "new_global_coefficient_or_weight_solves": 0, "accepted_PDE_steps": 0}
    (args.output/"registration.json").write_text(json.dumps(registration, indent=2)+"\n")
    source = args.root/"material_volume_candidate_v1"
    algebra_rows = []
    for dim, side in registration["algebra_cases"]:
        g = load_geometry(source/f"algebra_D{dim}_N{side}/geometry.json")
        a = operators(g)
        algebra_rows.append(algebra(g, a))
    (args.output/"algebra.json").write_text(json.dumps(algebra_rows, indent=2)+"\n")
    pressure_rows, scaling_controls = [], []
    for dim, primitive in [(1, 16), (2, 3), (3, 3)]:
        g = load_geometry(source/f"algebra_D{dim}_N{primitive}/geometry.json")
        native_double = load_geometry(source/f"algebra_D{dim}_N{2*primitive}/geometry.json")
        doubled = tile_geometry(g, 2)
        index = np.arange(len(doubled["V"]))
        cell = np.column_stack([(index % len(g["V"]) // primitive**a) % primitive for a in range(dim)])
        tile = np.column_stack([(index // len(g["V"]) // 2**a) % 2 for a in range(dim)])
        native_index = sum((cell[:, a]+primitive*tile[:, a])*(2*primitive)**a for a in range(dim))
        xerr = np.max(abs(doubled["x"]-native_double["x"][native_index]))
        verr = np.max(abs(doubled["V"]-native_double["V"][native_index]))
        assert max(xerr, verr) < 1e-12
        scaling_controls.append({"dimension": dim, "native_doubled_position_max_error": float(xerr), "native_doubled_volume_max_error": float(verr)})
        for repeats in registration["pressure_pattern_repeats"]:
            tiled = tile_geometry(g, repeats)
            op = operators(tiled)
            pressure_rows.append(pressure_patch(tiled, op, tiled["central_pattern"]))
    (args.output/"pressure_refinement.json").write_text(json.dumps(pressure_rows, indent=2)+"\n")
    failures = []
    for dim in (1, 2, 3):
        cases = [r for r in pressure_rows if r["dimension"] == dim]
        for axis in range(dim):
            errors = [next(p["force_density_error_L1"] for p in r["patches"] if p["axis"] == axis and p["representation"] == "exact_cell_average" and p["operator"] == "D") for r in cases]
            if errors[0] > 1e-9 and errors[1] >= .99*errors[0]:
                failures.append({"dimension": dim, "axis": axis, "coarse_error": errors[0], "fine_error": errors[1], "ratio": errors[1]/errors[0]})
    result = {"scope": __doc__, "algebra_pass": True, "scaling_controls": scaling_controls,
              "primary_pressure_nonconvergence_failures": failures, "pressure_gate_pass": not failures,
              "decision": "REJECT_AT_OPERATOR_PRESSURE_GATE" if failures else "PRIMARY_GATES_PASS_LATER_PHASES_REQUIRED",
              "weight_solve_released": False, "candidate_PDE_released": False,
              "Bloch_native_acoustic_phase_executed": False,
              "new_global_solves": 0, "new_native_meshes": 0, "accepted_steps": 0,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "geometry_sources": [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(source.glob("algebra*/geometry.json"))]}
    (args.output/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
