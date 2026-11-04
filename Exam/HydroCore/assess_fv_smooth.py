"""Assess saved native smooth flows using physical node and actual-cell references."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    result = json.loads((args.case/"result.json").read_text())
    registration = {"case": str(args.case), "reference": str(args.reference),
                    "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
                    "scope": "Read-only actual-cell analytic means and physical generator velocities; no PDE evolution"}
    save(args.output/"registration.json", registration)
    records = {}
    for state in ("initial", "final"):
        checkpoint = args.case/(state+".fv")
        command = [str(args.reference), str(checkpoint), str(args.output/(state+"_reference.csv"))]
        process = subprocess.run(command, text=True, capture_output=True, check=True,
                                 env=dict(os.environ, TMPDIR=str(args.output), OMP_NUM_THREADS="1"))
        save(args.output/(state+"_process.json"),
             {"command": command, "exit_code": process.returncode,
              "checkpoint_sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
              "reference_diagnostic": json.loads(process.stdout)})
        field = np.loadtxt(args.case/(state+".dat"))
        ref = np.genfromtxt(args.output/(state+"_reference.csv"), delimiter=",", names=True)
        assert np.array_equal(field[:, 0], ref["id"])
        volume = field[:, 10]
        assert np.max(abs(ref["V_ref"]/volume-1)) < 1e-10
        mean_error = np.hypot(field[:, 13]-ref["vx_mean_ref"], field[:, 14]-ref["vy_mean_ref"])
        pressure_error = abs(field[:, 12]-ref["P_mean_ref"])
        thermal = field[:, 9]-.5*np.sum(field[:, 6:9]**2, axis=1)/field[:, 5]
        records[state] = {
            "volume": float(volume.sum()),
            "density_L1": float(np.sum(volume*abs(field[:, 11]-ref["rho_ref"]))),
            "mass_mean_velocity_L1": float(np.sum(volume*mean_error)),
            "physical_mean_pressure_L1": float(np.sum(volume*pressure_error)),
            "physical_mean_pressure_relative_L1": float(np.sum(volume*pressure_error/ref["P_mean_ref"])),
            "quadrature_velocity_L1_difference": float(np.sum(volume*ref["quadrature_dv"])),
            "quadrature_pressure_L1_difference": float(np.sum(volume*ref["quadrature_dP"])),
            "thermal_mean_state_total": float(thermal.sum()),
            "minimum_pressure": float(field[:, 12].min()),
        }
        if state == "final":
            node = np.loadtxt(args.case/"final_velocity.dat")
            assert np.array_equal(node[:, 0], field[:, 0])
            records[state]["physical_node_velocity_L1"] = float(np.sum(volume*np.hypot(
                node[:, 1]-ref["vx_node_ref"], node[:, 2]-ref["vy_node_ref"])))
    initial, final = records["initial"], records["final"]
    assessment = {"case": args.case.name, "problem": result["problem"], "time": result["time"],
                  "endpoint_pass": not result["failed"] and abs(result["time"]-result["target"]) < 1e-13,
                  "initial": initial, "final": final,
                  "thermal_mean_state_total_change": final["thermal_mean_state_total"]-initial["thermal_mean_state_total"],
                  "GCL_energy_budget_over_initial_thermal": result["cumulative_GCL_energy_budget"]/initial["thermal_mean_state_total"],
                  "warning": "The GCL energy budget is a diagnostic, not a rigorous bound on pressure error."}
    if result["problem"] == "cold":
        assessment["cold_pressure_gate_pass"] = assessment["endpoint_pass"] and final["physical_mean_pressure_relative_L1"] <= .01
    save(args.output/"assessment.json", assessment)
    print(json.dumps(assessment, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
