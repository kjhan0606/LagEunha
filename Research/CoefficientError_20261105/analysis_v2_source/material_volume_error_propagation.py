"""Frozen coefficient-error propagation, never a fluid evolution program.

The native C probe keeps its original iteration/tolerance/shared-time caps.
The independent sparse reference is reused, not solved again. Positive graph
conductances define analytical residual-flow upper estimates. Extended-precision
evaluation monitors roundoff, but is not an interval-arithmetic proof.
"""
from pathlib import Path
import argparse
import hashlib
import json
import struct

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra

from material_volume_candidate_audit import assemble
from rotating_shear_assessment import read_native

LD = np.longdouble


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Graph:
    """Undirected positive Laplacian with one deterministic resistance tree."""

    def __init__(self, conductance):
        conductance = ((conductance + conductance.T) * 0.5).tocsr()
        assert np.max(abs(conductance.diagonal()), initial=0) == 0
        assert np.all(conductance.data > 0)
        self.n = conductance.shape[0]
        edges = sparse.triu(conductance, k=1).tocoo()
        self.i, self.j = edges.row, edges.col
        self.c = edges.data.astype(LD)
        resistance = conductance.copy()
        resistance.data = 1 / resistance.data
        distance, parent = dijkstra(resistance, indices=0, return_predecessors=True)
        assert np.isfinite(distance).all()
        self.parent = parent
        self.order = np.argsort(distance, kind="stable")
        self.tree_c = np.zeros(self.n, dtype=LD)
        self.path_resistance = np.zeros(self.n, dtype=LD)
        for node in self.order[1:]:
            ancestor = int(parent[node])
            assert ancestor >= 0 and distance[ancestor] < distance[node]
            self.tree_c[node] = conductance[node, ancestor]
            self.path_resistance[node] = (
                self.path_resistance[ancestor] + 1 / self.tree_c[node])
        self.path_mean = self.path_resistance.mean(dtype=LD)
        self.conductance = conductance

    def action(self, x):
        x = np.asarray(x, dtype=LD)
        edge = self.c * (x[self.i] - x[self.j])
        out = np.zeros(self.n, dtype=LD)
        np.add.at(out, self.i, edge)
        np.add.at(out, self.j, -edge)
        return out

    def estimate(self, x, rhs):
        """Le=r, e mean zero. Energy <= feasible residual-flow energy.

        |e_i|² <= E (R_path(i,root) + mean R_path(j,root)).
        A roundoff root closure makes the routed residual exactly compatible
        to working precision. The closure and gauge correction are reported.
        """
        x = np.asarray(x, dtype=LD)
        gauge_mean = x.mean(dtype=LD)
        x = x - gauge_mean
        residual = self.action(x) - rhs
        imbalance = residual.sum(dtype=LD)
        residual[0] -= imbalance
        flow = residual.copy()
        for node in self.order[:0:-1]:
            flow[self.parent[node]] += flow[node]
        energy = np.sum(flow[1:]**2 / self.tree_c[1:], dtype=LD)
        divergence = flow.copy()
        divergence[0] = 0
        for node in self.order[1:]:
            divergence[self.parent[node]] -= flow[node]
        flow_defect = np.max(abs(divergence - residual))
        assert flow_defect < 1e-14 * max(LD(1e-20), np.max(abs(residual)))
        # Conservative padding for ordinary positive arithmetic, not intervals.
        energy *= 1 + LD(4096) * np.finfo(LD).eps
        point = np.sqrt(energy * (self.path_resistance + self.path_mean))
        return x, residual, point, energy, {
            "residual_L2": float(np.linalg.norm(residual)),
            "residual_compatibility_root_closure": float(imbalance),
            "tree_flow_divergence_max_defect": float(flow_defect),
            "energy_upper_estimate": float(energy),
            "coefficient_pointwise_upper_max": float(point.max()),
            "coefficient_pointwise_upper_mean": float(point.mean()),
            "gauge_mean_removed": float(gauge_mean),
        }


def reference_flow(x):
    dx, dy = x[:, 0] - 0.5, x[:, 1] - 0.5
    radius = np.hypot(dx, dy)
    a = np.where(radius <= 0.2, 5.0,
                 np.where(radius <= 0.4, 2 / np.maximum(radius, 1e-100) - 5, 0))
    velocity = np.column_stack([-a * dy, a * dx])
    pressure = np.where(radius <= 0.2, 5 + 12.5 * radius**2,
                        np.where(radius <= 0.4,
                                 9 + 12.5 * radius**2 - 20 * radius
                                 + 4 * np.log(np.maximum(5 * radius, 1e-100)),
                                 3 + 4 * np.log(2)))
    gradient = np.column_stack([a*a*dx, a*a*dy])
    return radius, velocity, pressure, gradient


def small_mode_checks(root):
    results = []
    for dim in (2, 3):
        geometry = json.loads((root / f"algebra_D{dim}_N6/geometry.json").read_text())
        a = assemble(geometry)
        L = (a["L"] + a["L"].T) * 0.5
        values, vectors = np.linalg.eigh(L)
        Q = np.einsum("sab,iabk->isk", a["S"], a["H"])
        x = a["x"]
        velocity = np.column_stack([np.sin(2*np.pi*x[:, (d+1) % dim])
                                    * (1 if d == 0 else 0.3) for d in range(dim)]).ravel()
        pressure = 2 + 0.3 * np.cos(2*np.pi*x[:, 0]) + 0.2 * np.sin(2*np.pi*x[:, 1])
        s = np.einsum("ibk,k->ib", Q, velocity)
        graph = Graph(sparse.csr_matrix(np.diag(np.diag(L)) - L))
        for b in range(len(a["S"])):
            for mode, k in (("lowest_positive", 1), ("largest", len(values)-1)):
                e = vectors[:, k] / values[k]
                e -= e.mean()
                e /= np.linalg.norm(L @ e)
                r = L @ e
                volume = L @ (e * s[:, b])
                residual_only = s[:, b] * r
                commutator = volume - residual_only
                force = Q[:, b, :].T @ (e * (L @ pressure))
                affine_volume = L @ e  # Unit, constant affine strain component.
                _, _, point, energy, tree = graph.estimate(e, np.zeros(len(e), dtype=LD))
                observed_energy = LD(e) @ graph.action(e)
                assert observed_energy <= energy * (1 + 1e-12)
                assert np.max(abs(e) - point) < 1e-11
                work = abs(velocity @ force - pressure @ volume)
                assert abs(np.linalg.norm(r) - 1) < 1e-10
                assert np.max(abs(affine_volume-r)) < 1e-10
                assert work < 1e-9
                results.append({"dimension": dim, "cells": len(e), "basis": b,
                                "mode": mode, "eigenvalue": float(values[k]),
                                "residual_L2": float(np.linalg.norm(r)),
                                "coefficient_error_L2": float(np.linalg.norm(e)),
                                "volume_rate_error_L2": float(np.linalg.norm(volume)),
                                "residual_only_volume_term_L2": float(np.linalg.norm(residual_only)),
                                "commutator_L2": float(np.linalg.norm(commutator)),
                                "force_error_L2": float(np.linalg.norm(force)),
                                "force_work_adjoint_error": float(work),
                                "observed_energy": float(observed_energy),
                                "tree_estimate": tree})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("capture", type=Path)
    parser.add_argument("--output", type=Path, help="New analysis directory; preserve any earlier analysis.")
    args = parser.parse_args()
    root, capture = args.root, args.capture
    out = args.output if args.output is not None else capture
    if args.output is not None:
        out.mkdir()
    result_path = out / "propagation.json"
    assert not result_path.exists()
    native_input = root / "material_volume_native_cost_v3/material_volume_native_input.bin"
    with native_input.open("rb") as stream:
        magic, n, nf = struct.unpack("<QQQ", stream.read(24))
        assert magic == 0x4D56434F53543131
        cells = np.fromfile(stream, "<f8", count=6*n).reshape(n, 6)
        faces = np.fromfile(stream, dtype=np.dtype([
            ("i", "<i4"), ("j", "<i4"), ("rx", "<f8"), ("ry", "<f8"),
            ("cx", "<f8"), ("cy", "<f8"), ("area", "<f8")]), count=nf)
        assert not stream.read()
    with (capture / "coefficients.bin").open("rb") as stream:
        magic, count, mask = struct.unpack("<QQQ", stream.read(24))
        assert magic == 0x4D56434150543131 and count == n
        assert mask == 3, "Both fields are needed for the full-operator audit; never substitute an unattempted field."
        captured = np.fromfile(stream, "<f8", count=4*n).reshape(2, 2, n)
        inverse = np.fromfile(stream, "<f8", count=4*n).reshape(n, 2, 2)
        assert not stream.read()
    rhs = captured[:, 0, :].T.astype(LD)
    rhs_compatibility = rhs.sum(axis=0, dtype=LD)
    rhs -= rhs.mean(axis=0, dtype=LD)
    rhs[0] -= rhs.sum(axis=0, dtype=LD)
    approximate = captured[:, 1, :].T.astype(LD)
    reference = np.load(root / "material_volume_reference_v1/psi_reference.npy").astype(LD)
    assert reference.shape == approximate.shape == (n, 2)
    i, j = faces["i"], faces["j"]
    r = np.column_stack([faces["rx"], faces["ry"]])
    c = np.column_stack([faces["cx"], faces["cy"]])
    distance = np.linalg.norm(r, axis=1)
    directed_c = faces["area"] / (2*distance)
    raw_conductance = sparse.csr_matrix((directed_c, (i, j)), shape=(n, n))
    graph = Graph(raw_conductance)
    canonical = sparse.diags(np.asarray(graph.conductance.sum(axis=1)).ravel()) - graph.conductance
    raw = sparse.diags(np.bincount(i, weights=directed_c, minlength=n)) - raw_conductance
    defect = canonical - raw
    estimates = {"approximate": [], "reference": []}
    residuals, points, energies = {}, {}, {}
    for label, coefficients in (("approximate", approximate), ("reference", reference)):
        residuals[label] = np.zeros((n, 2), dtype=LD)
        points[label] = np.zeros((n, 2), dtype=LD)
        energies[label] = np.zeros(2, dtype=LD)
        for b in range(2):
            centered, residuals[label][:, b], points[label][:, b], energies[label][b], info = graph.estimate(coefficients[:, b], rhs[:, b])
            coefficients[:, b] = centered
            info["relative_residual"] = float(np.linalg.norm(residuals[label][:, b]) / np.linalg.norm(rhs[:, b]))
            estimates[label].append(info)
    e = approximate - reference
    e_residual = residuals["approximate"] - residuals["reference"]
    assert np.max(abs(e) - points["approximate"] - points["reference"]) <= 1e-20
    state, time = read_native(root / "rotating_shear_suite_v1/gresho_128x128/khout.000543.dat")
    assert time == 0.5
    assert np.array_equal(state["x"], cells[:, 0]) and np.array_equal(state["y"], cells[:, 1])
    assert np.array_equal(state["V"], cells[:, 2])
    assert np.all(state["m"] == 1/n) and np.all(state["W"] == 0)
    gamma = float(np.float32(5/3))
    assert np.max(abs((gamma-1)*state["U"]/state["V"] - cells[:, 3])) == 0
    V, mass = cells[:, 2], state["m"]
    radius, exact_v, exact_p, gradient_p = reference_flow(cells[:, :2])
    ell = (faces["area"] / distance**3)[:, None] * np.einsum("kab,kb->ka", inverse[i], r)
    S = np.array([[[1, 0], [0, -1]], [[0, 1], [1, 0]]]) / np.sqrt(2)
    qface = np.einsum("bad,kd->kba", S, ell).astype(LD)
    bj = faces["area"][:, None] * (r-c) / distance[:, None]
    bi = faces["area"][:, None] * c / distance[:, None]
    gi, gj, gc = graph.i, graph.j, graph.c

    def strain(velocity):
        values = np.einsum("kba,ka->kb", qface, velocity[j]-velocity[i])
        fitted = np.zeros((n, 2), dtype=LD)
        np.add.at(fitted, i, values)
        return fitted

    def base_volume(velocity):
        return np.bincount(i, weights=np.einsum("ka,ka->k", bj, velocity[j]-velocity[i]), minlength=n).astype(LD)

    def base_force(pressure):
        result = np.zeros((n, 2))
        np.add.at(result, i, (pressure[i]-pressure[j])[:, None]*bi)
        return result

    def force(coefficients, pressure):
        lp = graph.action(pressure)
        flux = np.einsum("kba,kb->ka", qface, coefficients[i]*lp[i, None])
        result = np.zeros((n, 2), dtype=LD)
        np.add.at(result, j, flux)
        np.add.at(result, i, -flux)
        return result

    def force_bound(point, pressure):
        lp = abs(graph.action(pressure))
        qsum = np.zeros((n, 2, 2), dtype=LD)
        np.add.at(qsum, i, qface)
        bound = lp * np.sum(point * np.linalg.norm(qsum, axis=2), axis=1)
        incoming = np.sum(point[i] * lp[i, None] * np.linalg.norm(qface, axis=2), axis=1)
        np.add.at(bound, j, incoming)
        return bound

    pair_c = np.asarray(graph.conductance[i, j]).ravel().astype(LD)
    assert np.all(pair_c > 0)
    # Same saved-state native dt is a reporting scale, not a candidate stability claim.
    dt = 0.00070072076050564647
    h = 1/np.sqrt(n)
    regions = [("all", np.ones(n, dtype=bool)), ("rigid_core_r_lt_015", radius < 0.15),
               ("smooth_annulus_025_to_035", (radius > 0.25) & (radius < 0.35)),
               ("outer_r_gt_045", radius > 0.45)]

    def regional(values, weights=V):
        return [{"region": name, "samples": int(selection.sum()),
                 "weighted_mean": float(np.sum(weights[selection]*values[selection]) / weights[selection].sum()),
                 "maximum": float(values[selection].max())} for name, selection in regions]

    flows = []
    # Verify the captured common-fit inverse against independent saved native diagnostics.
    transport = np.genfromtxt(root / "material_volume_native_cost_v3/unified_transport.csv", delimiter=",", names=True)
    fit_errors = []
    for name, velocity in (("actual", cells[:, 4:6]), ("exact_nodal_gresho", exact_v)):
        s = strain(velocity)
        native_G = np.zeros((n, 2, 2), dtype=LD)
        np.add.at(native_G, i, np.einsum("ka,kb->kab", velocity[j]-velocity[i], ell))
        native_div = native_G[:, 0, 0] + native_G[:, 1, 1]
        column = "div_fit" if name == "actual" else "div_fit_exact_velocity"
        fit_errors.append(float(np.max(abs(native_div - transport[column]))))
        assert fit_errors[-1] < 1e-10
        dw = np.sum(e*s, axis=1)
        volume = graph.action(dw)
        residual_only = np.sum(s*e_residual, axis=1)
        commutator = np.zeros(n, dtype=LD)
        for b in range(2):
            difference = s[gi, b]-s[gj, b]
            np.add.at(commutator, gi, gc*difference*e[gj, b])
            np.add.at(commutator, gj, -gc*difference*e[gi, b])
        decomposition = float(np.max(abs(volume - residual_only - commutator)))
        assert decomposition < 1e-12
        volume_bounds, face_bounds = {}, {}
        for label in ("approximate", "reference"):
            point, energy = points[label], energies[label]
            bound = np.sum(abs(s*residuals[label]), axis=1)
            face_bound = np.zeros(nf, dtype=LD)
            for b in range(2):
                difference = abs(s[gi, b]-s[gj, b])
                np.add.at(bound, gi, gc*difference*point[gj, b])
                np.add.at(bound, gj, gc*difference*point[gi, b])
                adjacent = np.sqrt(energy[b]/pair_c)
                left = abs(s[i, b])*adjacent + abs(s[i, b]-s[j, b])*point[j, b]
                right = abs(s[j, b])*adjacent + abs(s[i, b]-s[j, b])*point[i, b]
                nodal = abs(s[i, b])*point[i, b] + abs(s[j, b])*point[j, b]
                face_bound += np.minimum(np.minimum(left, right), nodal)
            volume_bounds[label], face_bounds[label] = bound, face_bound
        assert np.max(abs(volume) - volume_bounds["approximate"] - volume_bounds["reference"]) < 1e-14
        plane = abs(dw[i]-dw[j]) / (2*distance)
        plane_bound = face_bounds["approximate"]/(2*distance)
        assert np.max(plane - plane_bound - face_bounds["reference"]/(2*distance)) < 1e-14
        base = base_volume(velocity)
        candidate = base + graph.action(np.sum(reference*s, axis=1))
        flows.append({"velocity": name,
                      "observed_relative_volume_rate_error": regional(abs(volume)/V),
                      "relative_volume_rate_upper_estimate": regional(volume_bounds["approximate"]/V),
                      "reference_relative_volume_rate_uncertainty_upper": regional(volume_bounds["reference"]/V),
                      "one_step_relative_thermal_error_observed": regional((gamma-1)*dt*abs(volume)/V),
                      "one_step_relative_thermal_error_upper_estimate": regional((gamma-1)*dt*volume_bounds["approximate"]/V),
                      "residual_only_volume_term_L1": float(np.sum(abs(residual_only))),
                      "commutator_volume_term_L1": float(np.sum(abs(commutator))),
                      "commutator_decomposition_max_error": decomposition,
                      "base_relative_volume_rate_L1": float(np.sum(abs(base))/V.sum()),
                      "reference_candidate_relative_volume_rate_L1": float(np.sum(abs(candidate))/V.sum()),
                      "base_signed_volume_rate": float(base.sum()),
                      "reference_candidate_signed_volume_rate": float(candidate.sum()),
                      "fitted_physical_divergence_L1": float(np.sum(V*abs(native_div))/V.sum()),
                      "base_volume_rate_fit_gap_L1": float(np.sum(abs(base-V*native_div))/V.sum()),
                      "reference_candidate_volume_rate_fit_gap_L1": float(np.sum(abs(candidate-V*native_div))/V.sum()),
                      "one_step_face_normal_displacement_over_h": {
                          "observed_area_weighted_mean": float(np.average(dt*plane/h, weights=faces["area"])),
                          "observed_maximum": float(np.max(dt*plane/h)),
                          "upper_area_weighted_mean": float(np.average(dt*plane_bound/h, weights=faces["area"])),
                          "upper_maximum": float(np.max(dt*plane_bound/h))}})

    pressures = []
    for name, pressure in (("actual", cells[:, 3]), ("exact_nodal_gresho", exact_p)):
        df = force(e, pressure)
        upper = force_bound(points["approximate"], pressure)
        ref_upper = force_bound(points["reference"], pressure)
        norm = np.linalg.norm(df, axis=1)
        assert np.max(norm-upper-ref_upper) < 1e-14
        checks = []
        for velocity in (cells[:, 4:6], exact_v):
            lhs = np.sum(velocity*df)
            rhs_work = pressure @ graph.action(np.sum(e*strain(velocity), axis=1))
            checks.append(float(abs(lhs-rhs_work)))
        assert max(checks) < 1e-12
        record = {"pressure": name,
                  "observed_acceleration_error": regional(norm/mass),
                  "acceleration_error_upper_estimate": regional(upper/mass),
                  "reference_acceleration_uncertainty_upper": regional(ref_upper/mass),
                  "one_step_velocity_error_observed": regional(dt*norm/mass),
                  "one_step_velocity_error_upper_estimate": regional(dt*upper/mass),
                  "force_work_adjoint_errors": checks,
                  "total_force_error_vector": np.asarray(df.sum(axis=0), dtype=float).tolist()}
        if name == "exact_nodal_gresho":
            native_force = base_force(pressure)
            corrected_force = native_force + force(reference, pressure)
            record["base_force_density_error_from_analytic_gradient"] = regional(np.linalg.norm(native_force/V[:, None]+gradient_p, axis=1))
            record["reference_candidate_force_density_error_from_analytic_gradient"] = regional(np.linalg.norm(corrected_force/V[:, None]+gradient_p, axis=1))
        pressures.append(record)

    velocity_noise = cells[:, 4:6] - exact_v
    pressure_noise = cells[:, 3] - exact_p
    base_noise_volume = base_volume(velocity_noise)
    candidate_noise_volume = base_noise_volume + graph.action(np.sum(reference*strain(velocity_noise), axis=1))
    base_noise_force = base_force(pressure_noise)
    candidate_noise_force = base_noise_force + force(reference, pressure_noise)
    frozen_noise = {
        "scope": "Linear response on a saved original-operator state. Noise means deviation from the nodal Gresho reference. Not a candidate trajectory or instability proof.",
        "base_velocity_noise_relative_volume_rate_L1": float(np.sum(abs(base_noise_volume))/V.sum()),
        "reference_candidate_velocity_noise_relative_volume_rate_L1": float(np.sum(abs(candidate_noise_volume))/V.sum()),
        "base_pressure_noise_force_density_L1": float(np.sum(np.linalg.norm(base_noise_force, axis=1))/V.sum()),
        "reference_candidate_pressure_noise_force_density_L1": float(np.sum(np.linalg.norm(candidate_noise_force, axis=1))/V.sum()),
    }

    np.savez_compressed(out / "propagation_arrays.npz", coefficient_difference=np.asarray(e, dtype=float),
                        approximate_error_pointwise_upper=np.asarray(points["approximate"], dtype=float),
                        reference_error_pointwise_upper=np.asarray(points["reference"], dtype=float),
                        residual_approximate=np.asarray(residuals["approximate"], dtype=float),
                        residual_reference=np.asarray(residuals["reference"], dtype=float))
    result = {"scope": __doc__, "cells": n, "directed_faces": nf, "checkpoint_time": time,
              "same_saved_native_dt_reporting_scale": dt, "candidate_dt_stability_qualified": False,
              "gamma": gamma, "mass_equal_1_over_n_verified": True,
              "raw_to_symmetric_graph_max_entry_defect": float(np.max(abs(defect.data), initial=0)),
              "canonicalization_scope": "Analysis graph only. Reciprocal positive conductances averaged; no face discarded. Runtime graph unchanged.",
              "rhs_mean_closure_before_analysis": np.asarray(rhs_compatibility, dtype=float).tolist(),
              "extended_precision_mantissa_bits": np.finfo(LD).nmant,
              "bounds_are_interval_arithmetic_certificates": False,
              "bound_roundoff_scope": "Analytical feasible-flow estimates in long double with monitored residual closure and positive-arithmetic padding. Not a computer-assisted interval proof.",
              "coefficient_error_estimates": estimates,
              "coefficient_difference_relative_L2": np.asarray(np.linalg.norm(e, axis=0)/np.linalg.norm(reference, axis=0), dtype=float).tolist(),
              "coefficient_difference_L2": np.asarray(np.linalg.norm(e, axis=0), dtype=float).tolist(),
              "common_native_fit_divergence_max_errors": fit_errors,
              "velocity_cases": flows, "pressure_cases": pressures,
              "frozen_noise_response": frozen_noise,
              "small_equal_residual_mode_checks": small_mode_checks(root / "material_volume_candidate_v1"),
              "scientific_error_budget_numerically_declared": False,
              "new_tolerance_selected": False, "production_solver_released": False,
              "accepted_steps": 0,
              "source_sha256": digest(__file__), "capture_sha256": digest(capture / "coefficients.bin"),
              "native_input_sha256": digest(native_input)}
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"cells": n, "coefficient_difference_relative_L2": result["coefficient_difference_relative_L2"],
                      "volume_rates": [{"velocity": r["velocity"], "observed": r["observed_relative_volume_rate_error"][0],
                                        "upper": r["relative_volume_rate_upper_estimate"][0]} for r in flows],
                      "accelerations": [{"pressure": r["pressure"], "observed": r["observed_acceleration_error"][0],
                                         "upper": r["acceleration_error_upper_estimate"][0]} for r in pressures],
                      "small_mode_checks": len(result["small_equal_residual_mode_checks"]),
                      "production_solver_released": False}, indent=2), flush=True)


if __name__ == "__main__":
    main()
