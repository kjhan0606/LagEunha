"""Independent original-registry/IC checks and moving guard-boundary controls.

Initialization checks explicitly perform no evolution. Only the four
contact boundary controls are trajectories; this is not a shock-suite pass.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from run_fv_native_suite import save


def run(binary, root, name, problem, target, guard, initial_only):
    case = root/name
    case.mkdir(exist_ok=False)
    env = dict(os.environ, TMPDIR=str(root), OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE",
               GFS_TEST_GUARD_MULTIPLIER=str(guard))
    env.pop("GFS_TEST_INIT_ONLY", None)
    if initial_only:
        env["GFS_TEST_INIT_ONLY"] = "1"
    command = [str(binary), "1", "200", problem, str(target)]
    with (case/"stdout.log").open("x") as out, (case/"stderr.log").open("x") as err:
        process = subprocess.run(command, cwd=case, env=env, stdout=out, stderr=err)
    save(case/"process.json", {"command": command, "exit_code": process.returncode,
                               "initialization_only": initial_only, "guard_multiplier": guard})
    if process.returncode:
        raise RuntimeError(f"Native boundary/registry check failed: {case}")
    return case


def segments(left, right, cuts):
    points = [left]+[c for c in cuts if left < c < right]+[right]
    return zip(points[:-1], points[1:])


def initial_check(case, spec):
    meta = json.loads((case/"initialization_check.json").read_text())
    assert meta["evolution_performed"] is False and not (case/"result.json").exists()
    assert meta["gamma"] == spec["gamma"]
    flags = meta["bc_flags"][0]
    lower, upper = (flags & 255, (flags >> 8) & 255) if flags & 0x10000 else (flags, flags)
    assert (lower == 1) == spec["lower_reflect"] and (upper == 1) == spec["upper_reflect"]
    if (case/"problem_spec.json").exists():
        stated = json.loads((case/"problem_spec.json").read_text())
        assert stated["physical_domain"] == spec["domain"]
        assert stated["original_endpoint"] == spec["endpoint"]
        assert stated["regions"] == spec["regions"]
        for k, state in enumerate(spec["states"]):
            assert state == [stated["rho"][k], stated["velocity"][k], stated["pressure"][k]]
    data = np.loadtxt(case/"initial.dat")
    cuts = [spec["cut"]] if spec["regions"] >= 2 else []
    if spec["regions"] == 3:
        cuts.append(spec["cut2"])
    node, weight = np.polynomial.legendre.leggauss(8)
    errors = np.zeros(4)
    for row in data:
        center, volume = row[1]+row[16], row[10]
        left, right = center-.5*volume, center+.5*volume
        mass = momentum = pressure = 0.
        for a, b in segments(left, right, cuts):
            x = .5*(a+b)+.5*(b-a)*node
            region = int(np.searchsorted(cuts, .5*(a+b), side="right"))
            rho, velocity, p = spec["states"][region]
            densities = 1+.2*np.sin(5*np.pi*x) if spec["problem"] == "shu" and region == 1 else np.full(8, rho)
            dm = .5*(b-a)*np.dot(weight, densities)
            mass += dm
            momentum += velocity*dm
            pressure += p*(b-a)
        energy = pressure/(spec["gamma"]-1)+.5*momentum**2/mass
        errors = np.maximum(errors, [abs(row[5]-mass)/mass,
                                     abs(row[6]-momentum)/(mass+abs(momentum)),
                                     abs(row[9]-energy)/energy,
                                     abs(row[19]-row[20]/volume**2)/(1+abs(row[19]))])
    assert np.max(errors) < 1e-11
    return {"problem": spec["problem"], "scope": "Original legacy setup plus independent split-Gauss actual-cell IC integration; no PDE endpoint",
            "cells": len(data), "maximum_normalized_mass_momentum_energy_mu_errors": errors.tolist(), "pass": True}


def contact_control(case, moving):
    data = np.loadtxt(case/"final.dat")
    result = json.loads((case/"result.json").read_text())
    assert not result["failed"] and result["time"] == .25
    center = data[:, 1]+data[:, 16]
    left, right = np.maximum(center-.5*data[:, 10], 0), np.minimum(center+.5*data[:, 10], 1)
    width = np.maximum(0, right-left)
    assert abs(width.sum()-1) < 1e-12
    contact = .5+(.25 if moving else 0)
    mass_exact = np.maximum(0, np.minimum(right, contact)-left)+.125*np.maximum(0, right-np.maximum(left, contact))
    density_error = float(np.sum(abs(width*data[:, 11]-mass_exact)))
    pressure_error = float(np.sum(width*abs(data[:, 12]-1)))
    velocity_error = float(np.sum(width*abs(data[:, 13]-(1 if moving else 0))))
    assert max(density_error, pressure_error, velocity_error) < 1e-9
    return {"case": case.name, "density_L1_on_original_ROI": density_error,
            "pressure_L1_on_original_ROI": pressure_error, "velocity_L1_on_original_ROI": velocity_error,
            "boundary_ledger_error": result["maximum_ledger_error"],
            "mass_weight_with_residual": result["maximum_integrated_residual_invariant_error"], "pass": True}


def padding_comparison(first, second):
    a, b = np.loadtxt(first/"final.dat"), np.loadtxt(second/"final.dat")
    left_a, right_a = a[:, 1]+a[:, 16]-.5*a[:, 10], a[:, 1]+a[:, 16]+.5*a[:, 10]
    left_b, right_b = b[:, 1]+b[:, 16]-.5*b[:, 10], b[:, 1]+b[:, 16]+.5*b[:, 10]
    difference = np.zeros(5)
    i = j = 0
    while i < len(a) and j < len(b):
        width = max(0, min(right_a[i], right_b[j], 1)-max(left_a[i], left_b[j], 0))
        if width:
            x = [a[i, 11], a[i, 12], a[i, 13], (a[i, 5]-a[i, 20])/a[i, 20], a[i, 4]/(a[i, 20]/a[i, 19])]
            y = [b[j, 11], b[j, 12], b[j, 13], (b[j, 5]-b[j, 20])/b[j, 20], b[j, 4]/(b[j, 20]/b[j, 19])]
            difference += width*abs(np.asarray(x)-np.asarray(y))
        if right_a[i] < right_b[j]:
            i += 1
        else:
            j += 1
    assert max(difference) < 1e-9
    return {"cases": [first.name, second.name],
            "ROI_density_pressure_velocity_linked_mass_weight_over_h0square_L1_differences": difference.tolist(),
            "scope": "Contact boundary control only; not proof of shock ROI padding independence", "pass": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("legacy_registry", type=Path)
    args = parser.parse_args()
    if not args.root.resolve().is_relative_to(Path("/gpfs/kjhan/LagForce")):
        raise ValueError("Use GPFS")
    args.root.mkdir(exist_ok=False)
    specs = json.loads(args.legacy_registry.read_text())
    save(args.root/"registration.json", {"scope": "Nine independent IC-only checks, four actual contact boundary controls; no full shock-suite pass",
                                       "binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                                       "legacy_registry": str(args.legacy_registry), "physical_resolution": 200,
                                       "IC_normalized_limit": 1e-11, "contact_ROI_and_padding_limit": 1e-9,
                                       "OMP_NUM_THREADS": 1, "no_parameter_tuning": True})
    records = []
    for spec in specs:
        case = run(args.binary, args.root, "IC_"+spec["problem"], spec["problem"], spec["endpoint"], 1, True)
        records.append(initial_check(case, spec))
    boundary, padding = [], []
    for problem in ["contact_open", "moving_contact_open"]:
        cases = [run(args.binary, args.root, problem+"_guard"+str(g), problem, .25, g, False) for g in [1, 2]]
        boundary.extend(contact_control(case, problem.startswith("moving")) for case in cases)
        padding.append(padding_comparison(*cases))
    save(args.root/"assessment.json", {"original_IC_checks": records, "contact_boundary_evolution": boundary,
                                      "contact_padding_controls": padding, "pass": True,
                                      "unqualified": "Strong-shock endpoint, shock guard independence, full original nine-test evolution and general multidimensional open boundaries"})
    print(json.dumps({"IC_only_passes": len(records), "actual_contact_endpoints": len(boundary), "padding_controls": len(padding), "pass": True}))


if __name__ == "__main__":
    main()
