"""Bloch and known-wave response of frozen native operators, zero PDE steps.

Reuse the original irregular unit geometries and the saved N128 coefficient
reference. No new native coefficient solve, limiter, or production change.
Uniform-pressure acoustic spectra concern the positive adjoint block only.
"""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, eigsh, ArpackNoConvergence

from material_volume_candidate_audit import assemble
from rotating_shear_assessment import read_native


def unit_symbol(geometry, a, wave):
    dim, n = geometry["dimension"], geometry["cells"]
    B = np.zeros((n, n*dim), dtype=complex)
    L = np.zeros((n, n), dtype=complex)
    R = np.zeros_like(B)
    for f in geometry["faces"]:
        i, j = f["i"], f["j"]
        r, bj = np.array(f["r"]), np.array(f["bj"])
        shift = r - (a["x"][j]-a["x"][i])
        assert np.max(abs(shift-np.rint(shift))) < 1e-12
        phase = np.exp(1j*np.dot(wave, shift))
        distance = np.linalg.norm(r)
        ell = f["area"]/distance**(dim+1) * (np.array(geometry["native_inverse_fit"])[i] @ r)
        rv = np.einsum("b,bad,d->a", a["psi"][i], a["S"], ell)
        B[i, dim*j:dim*(j+1)] += bj*phase
        B[i, dim*i:dim*(i+1)] -= bj
        R[i, dim*j:dim*(j+1)] += rv*phase
        R[i, dim*i:dim*(i+1)] -= rv
        coefficient = f["area"]/(2*distance)
        L[i, i] += coefficient
        L[i, j] -= coefficient*phase
    assert np.max(abs(L-L.conj().T)) < 1e-11
    return B, B + L @ R


def acoustic(X, V, wave, x):
    dim = x.shape[1]
    # rho0=1, gamma P0=1. Therefore cs=1, M=diag(V), H=diag(1/V).
    weighted = X / np.sqrt(np.repeat(V, dim))[None, :] / np.sqrt(V)[:, None]
    A = weighted @ weighted.conj().T
    values, directions = np.linalg.eigh(A)
    assert np.max(abs(A-A.conj().T)) < 1e-10
    assert values[0] > -4096*np.finfo(float).eps*max(1, values[-1])
    profile = np.sqrt(V) * np.exp(1j*(x @ wave))
    profile /= np.linalg.norm(profile)
    overlaps = abs(directions.conj().T @ profile)**2
    branch = int(np.argmax(overlaps))
    return {"lowest_squared_frequencies": values[:min(5, len(values))].tolist(),
            "maximum_squared_frequency": float(values[-1]),
            "numerical_pressure_null_modes": int(np.sum(abs(values) < 4096*np.finfo(float).eps*max(1, values[-1]))),
            "pressure_plane_wave_branch_overlap": float(overlaps[branch]),
            "pressure_plane_wave_branch_squared_frequency": float(values[branch]),
            "pressure_plane_wave_branch_speed": float(np.sqrt(max(0, values[branch]))/np.linalg.norm(wave)) if np.linalg.norm(wave) else None}


def small_patterns(root):
    symbols, refinement, controls = [], [], []
    for dim, primitive in ((1, 16), (2, 3), (3, 3)):
        geometry = json.loads((root / f"algebra_D{dim}_N{primitive}/geometry.json").read_text())
        a = assemble(geometry)
        n, V, x = geometry["cells"], a["V"], a["x"]
        B0, C0 = unit_symbol(geometry, a, np.zeros(dim))
        assert np.max(abs(B0-a["B"])) < 1e-11
        assert np.max(abs(C0-a["C"])) < 1e-11
        # Reconstruct the archived finite-matrix definition; no new native graph.
        control = {"dimension": dim, "cells": n, "zero_wave_B_error": float(np.max(abs(B0-a["B"]))),
                   "zero_wave_C_error": float(np.max(abs(C0-a["C"]))),
                   "constant_pressure_B_force_max": float(np.max(abs(B0.conj().T @ np.ones(n)))),
                   "constant_pressure_C_force_max": float(np.max(abs(C0.conj().T @ np.ones(n))))}
        if dim == 1:
            control["one_dimensional_W0_C_minus_B_max"] = float(np.max(abs(C0-B0)))
            assert control["one_dimensional_W0_C_minus_B_max"] == 0
        else:
            repeated = json.loads((root / f"algebra_D{dim}_N6/geometry.json").read_text())
            repeated_a = assemble(repeated)
            index = np.arange(6**dim)
            tile = sum(((index // (6**d)) % 6 % 3) * 3**d for d in range(dim))
            mismatch = np.max(abs(repeated_a["psi"] - 0.25*a["psi"][tile]))
            assert mismatch < 1e-11
            control["repeated_pattern_scaled_coefficient_max_error"] = float(mismatch)
        controls.append(control)
        directions = [("axis", np.eye(dim)[0])]
        if dim > 1:
            directions.append(("diagonal", np.ones(dim)/np.sqrt(dim)))
        for name, direction in directions:
            for magnitude in (0.0, 0.01, 0.02, 0.04, np.pi/4, np.pi/2, np.pi):
                wave = magnitude*direction
                B, C = unit_symbol(geometry, a, wave)
                row = {"dimension": dim, "direction": name, "wave_norm": magnitude,
                       "B": acoustic(B, V, wave, x), "C": acoustic(C, V, wave, x)}
                if magnitude:
                    plane = np.exp(1j*(x @ wave))
                    longitudinal = (plane[:, None]*direction).ravel()
                    for label, X in (("B", B), ("C", C)):
                        truth = 1j*magnitude*V*plane
                        row[label]["longitudinal_nodal_divergence_relative_L1"] = float(np.sum(abs(X@longitudinal-truth))/(magnitude*V.sum()))
                    if dim > 1:
                        transverse_direction = np.zeros(dim)
                        transverse_direction[0], transverse_direction[1] = -direction[1], direction[0]
                        if np.linalg.norm(transverse_direction) == 0:
                            transverse_direction[1] = 1
                        transverse_direction /= np.linalg.norm(transverse_direction)
                        transverse = (plane[:, None]*transverse_direction).ravel()
                        for label, X in (("B", B), ("C", C)):
                            row[label]["transverse_false_compression_normalized_L1"] = float(np.sum(abs(X@transverse))/(magnitude*V.sum()))
                symbols.append(row)
            if dim > 1:
                physical_wave = 2*np.pi*(np.eye(dim)[0] if name == "axis" else np.ones(dim))
                normal = physical_wave/np.linalg.norm(physical_wave)
                transverse = np.zeros(dim)
                transverse[0], transverse[1] = -normal[1], normal[0]
                transverse /= np.linalg.norm(transverse)
                for side in (12, 24, 48, 96):
                    scale = 3/side
                    wave = scale*physical_wave
                    B, C = unit_symbol(geometry, a, wave)
                    velocity = (np.exp(1j*(x @ wave))[:, None]*transverse).ravel()
                    refinement.append({"dimension": dim, "direction": name, "global_side": side,
                                       "unit_pattern_physical_scale": scale,
                                       "physical_wave": physical_wave.tolist(),
                                       "B_transverse_false_compression_normalized_L1": float(np.sum(abs(B@velocity))/(np.linalg.norm(wave)*V.sum())),
                                       "C_transverse_false_compression_normalized_L1": float(np.sum(abs(C@velocity))/(np.linalg.norm(wave)*V.sum()))})
    return {"controls": controls, "symbols": symbols, "fixed_wavelength_refinement": refinement,
            "scope": "Bloch phase uses periodic image shifts and conjugate adjoints. The acoustic plane branch is identified by overlap. Highest frequencies are finite sampled values, not global symbol maxima. Refinement rescales the same repeated irregular native unit pattern analytically."}


def bounded_acoustic(X, V, mass, pressure, gamma):
    sqrtH = np.sqrt(gamma*pressure/V)
    deadline = time.monotonic()+15
    count = 0

    def action(v):
        nonlocal count
        count += 1
        if count > 256 or time.monotonic() > deadline:
            raise RuntimeError("fixed diagnostic matvec/time cap reached")
        return sqrtH*(X @ ((X.T @ (sqrtH*v))/mass))

    operator = LinearOperator((len(V), len(V)), matvec=action, dtype=float)
    result = {"maximum_matvecs": 256, "wall_budget_seconds": 15,
              "full_background_stability_or_KDK_dt_bound": False}
    mark = time.monotonic()
    try:
        values, vectors = eigsh(operator, k=1, which="LA", v0=np.cos(0.17*np.arange(len(V))),
                                tol=1e-8, maxiter=32, ncv=16)
        residual = np.linalg.norm(action(vectors[:, 0])-values[0]*vectors[:, 0])
        result.update(status="converged_ritz_estimate", squared_frequency=float(values[0]),
                      frequency=float(np.sqrt(values[0])), ritz_absolute_residual=float(residual))
    except (RuntimeError, ArpackNoConvergence) as error:
        result.update(status="diagnostic_cap_or_nonconvergence", reason=str(error))
    result.update(matvecs=count, seconds=time.monotonic()-mark)
    return result


def native_probes(root):
    input_path = root / "material_volume_native_cost_v3/material_volume_native_input.bin"
    with input_path.open("rb") as stream:
        magic, n, nf = struct.unpack("<QQQ", stream.read(24))
        assert magic == 0x4D56434F53543131
        cells = np.fromfile(stream, "<f8", count=6*n).reshape(n, 6)
        faces = np.fromfile(stream, dtype=np.dtype([
            ("i", "<i4"), ("j", "<i4"), ("rx", "<f8"), ("ry", "<f8"),
            ("cx", "<f8"), ("cy", "<f8"), ("area", "<f8")]), count=nf)
        assert not stream.read()
    with (root / "material_volume_error_propagation_v1/coefficients.bin").open("rb") as stream:
        magic, count, mask = struct.unpack("<QQQ", stream.read(24))
        assert magic == 0x4D56434150543131 and count == n and mask == 3
        stream.seek(4*n*8, 1)
        inverse = np.fromfile(stream, "<f8", count=4*n).reshape(n, 2, 2)
        assert not stream.read()
    psi = np.load(root / "material_volume_reference_v1/psi_reference.npy")
    state, time_value = read_native(root / "rotating_shear_suite_v1/gresho_128x128/khout.000543.dat")
    assert time_value == 0.5 and np.array_equal(cells[:, 0], state["x"]) and np.array_equal(cells[:, 2], state["V"])
    V, x = cells[:, 2], cells[:, :2]
    i, j = faces["i"], faces["j"]
    r = np.column_stack([faces["rx"], faces["ry"]])
    c = np.column_stack([faces["cx"], faces["cy"]])
    area = faces["area"]
    distance = np.linalg.norm(r, axis=1)
    bj = area[:, None]*(r-c)/distance[:, None]
    ell = (area/distance**3)[:, None]*np.einsum("kab,kb->ka", inverse[i], r)
    kvec = np.column_stack([psi[i, 0]*ell[:, 0]+psi[i, 1]*ell[:, 1],
                            psi[i, 1]*ell[:, 0]-psi[i, 0]*ell[:, 1]])/np.sqrt(2)

    def matrix(vectors):
        row = np.repeat(i, 2)
        column = 2*np.repeat(j, 2)+np.tile([0, 1], nf)
        own = 2*np.repeat(i, 2)+np.tile([0, 1], nf)
        return sparse.csr_matrix((np.concatenate([vectors.ravel(), -vectors.ravel()]),
                                  (np.tile(row, 2), np.concatenate([column, own]))), shape=(n, 2*n))

    B, R = matrix(bj), matrix(kvec)
    weights = area/(2*distance)
    conductance = sparse.csr_matrix((weights, (i, j)), shape=(n, n))
    conductance = (conductance+conductance.T)*0.5
    L = sparse.diags(np.asarray(conductance.sum(axis=1)).ravel())-conductance
    C = (B+L@R).tocsr()
    native = np.genfromtxt(root / "material_volume_native_cost_v3/unified_transport.csv", delimiter=",", names=True)
    base_rate_error = float(np.max(abs(B@cells[:, 4:6].ravel()/V-native["div_geometric"])))
    assert base_rate_error < 1e-10
    prior = json.loads((root / "material_volume_error_propagation_v1/analysis_v2/propagation.json").read_text())
    candidate_rate = float(np.sum(abs(C@cells[:, 4:6].ravel()))/V.sum())
    assert abs(candidate_rate-prior["velocity_cases"][0]["reference_candidate_relative_volume_rate_L1"]) < 1e-10
    normal = r/distance[:, None]
    tangent = np.column_stack([-normal[:, 1], normal[:, 0]])
    face_center = x[i]+c

    def average_plane(wave):
        # Exact line integral on each straight native polygon edge.
        face_integral = area*np.exp(1j*(face_center@wave))*np.sinc((tangent@wave)*area/(2*np.pi))
        integral = np.zeros(n, dtype=complex)
        np.add.at(integral, i, (normal@wave)*face_integral)
        average = integral/(1j*np.dot(wave, wave)*V)
        assert np.max(abs(average)) < 1+1e-10
        return average

    probes = []
    for direction, axis in (("axis", np.array([1, 0])), ("diagonal", np.array([1, 1]))):
        for mode in ((4, 8, 16, 64) if direction == "axis" else (4, 8, 16)):
            wave = 2*np.pi*mode*axis
            norm = np.linalg.norm(wave)
            longitudinal = wave/norm
            transverse = np.array([-longitudinal[1], longitudinal[0]])
            plane = np.exp(1j*(x@wave))
            averaged = average_plane(wave)
            pressure_representations = [("generator_sample", plane.imag), ("exact_cell_average", averaged.imag)]
            entry = {"direction": direction, "integer_wave_components": (mode*axis).tolist(),
                     "wavelength": float(2*np.pi/norm), "amplitude": 1,
                     "high_frequency_diagnostic_only": mode == 64,
                     "volume_normalization": "|k| amplitude total_volume; transverse truth exactly zero",
                     "pressure_force_reference": "Exact polygon-integrated pressure gradient. Generator and exact-cell-average pressure inputs reported separately.",
                     "responses": {}}
            for label, X in (("B", B), ("C", C)):
                measured = {}
                for kind, vector in (("transverse", transverse), ("longitudinal", longitudinal)):
                    velocity = plane.imag[:, None]*vector
                    truth = (np.dot(wave, vector)*V*averaged.real) if kind == "longitudinal" else np.zeros(n)
                    error = X@velocity.ravel()-truth
                    measured[kind+"_integrated_divergence_error_normalized_L1"] = float(np.sum(abs(error))/(norm*V.sum()))
                    measured[kind+"_integrated_divergence_error_normalized_max"] = float(np.max(abs(error)/V)/norm)
                for representation, pressure in pressure_representations:
                    force = (X.T@pressure).reshape(n, 2)
                    truth = -V[:, None]*averaged.real[:, None]*wave
                    error = np.linalg.norm(force-truth, axis=1)
                    measured[representation+"_integrated_pressure_force_error_normalized_L1"] = float(error.sum()/(norm*V.sum()))
                entry["responses"][label] = measured
            probes.append(entry)
    gamma = float(np.float32(5/3))
    stiffness = {label: bounded_acoustic(X, V, state["m"], cells[:, 3], gamma) for label, X in (("B", B), ("C", C))}
    if all(row["status"] == "converged_ritz_estimate" for row in stiffness.values()):
        stiffness["positive_block_frequency_ratio_C_over_B"] = stiffness["C"]["frequency"]/stiffness["B"]["frequency"]
    return {"cells": n, "checkpoint_time": time_value, "B_nonzeros": B.nnz, "R_nonzeros": R.nnz, "C_nonzeros": C.nnz,
            "operator_reproduction_max_base_rate_error": base_rate_error,
            "reference_candidate_actual_velocity_volume_rate_L1": candidate_rate,
            "known_wave_probes": probes, "positive_acoustic_block": stiffness,
            "pressure_integration_scope": "Analytical Fourier edge integrals and divergence theorem on original 2D straight-edge cells. Analysis only; no runtime subcell quadrature.",
            "full_stability_qualified": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir()
    registration = {"scope": __doc__, "primitive_dimensions": [1, 2, 3],
                    "small_symbol_wave_norms": [0, 0.01, 0.02, 0.04, np.pi/4, np.pi/2, np.pi],
                    "directions": ["axis", "diagonal"], "fixed_wavelength_refinement_sides": [12, 24, 48, 96],
                    "native_modes": [4, 8, 16], "native_high_frequency_diagnostic_axis_mode": 64,
                    "stiffness_max_matvecs_per_operator": 256, "stiffness_seconds_per_operator": 15,
                    "stop_rule": "Hold physics and solver tuning if long-wave speed has a wrong limit, transverse error does not decrease, or response is ambiguous. Larger compressive response alone is not rejection. No limiter, PDE, coefficient tuning or new native global reference solve.",
                    "previous_CPU_cost_and_tolerance_failure_preserved": True}
    (args.output/"registration.json").write_text(json.dumps(registration, indent=2)+"\n")
    small = small_patterns(args.root / "material_volume_candidate_v1")
    (args.output/"small_patterns.json").write_text(json.dumps(small, indent=2)+"\n")
    print("Small periodic symbols and analytic repeated-pattern refinement completed.", flush=True)
    native = native_probes(args.root)
    (args.output/"native_probes.json").write_text(json.dumps(native, indent=2)+"\n")
    result = {"scope": __doc__, "small_symbol_cases": len(small["symbols"]),
              "fixed_wavelength_refinement_cases": len(small["fixed_wavelength_refinement"]),
              "native_probe_cases": len(native["known_wave_probes"]),
              "accepted_steps": 0, "new_native_global_reference_solves": 0,
              "coefficient_solver_released": False, "production_physics_changed": False,
              "full_stability_qualified": False,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.output/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)
    print(json.dumps(native["positive_acoustic_block"], indent=2), flush=True)


if __name__ == "__main__":
    main()
