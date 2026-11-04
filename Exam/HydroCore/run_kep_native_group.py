"""First actual KEP coupling/vortex/cold group, with immutable settings.

Original blast/collision/open-boundary problems remain required; this
group alone cannot qualify the full candidate. Native C evolves the PDE.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
from queue import Queue

import numpy as np

from run_fv_native_suite import execute, save


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("binary", type=Path)
    args = parser.parse_args()
    if not args.root.resolve().is_relative_to(Path("/gpfs/kjhan/LagForce")):
        raise ValueError("Use the authorized GPFS project root")
    cpu_sets = [list(range(8)), list(range(8, 16))]
    if not set(range(16)) <= os.sched_getaffinity(0):
        raise ValueError("Registered disjoint CPU sets are unavailable")
    args.root.mkdir(exist_ok=False)
    cases = [("D1_contact_N200", 1, 200, "contact", .25),
             ("D1_periodic_moving_contact_N200", 1, 200, "moving_contact", .25),
             ("D2_contact_N8", 2, 8, "contact", .25),
             ("D2_periodic_moving_contact_N8", 2, 8, "moving_contact", .25),
             ("D3_contact_N4", 3, 4, "contact", .25),
             ("D3_periodic_moving_contact_N4", 3, 4, "moving_contact", .25),
             ("D1_sod_N200", 1, 200, "sod", .2),
             ("D1_sod_N400", 1, 400, "sod", .2),
             ("D2_cold_N64", 2, 64, "cold", .25),
             ("D2_cold_N64_quarterdt", 2, 64, "cold", .25, .25),
             ("D2_gresho_N64", 2, 64, "gresho", 1.)]
    save(args.root/"registration.json", {
        "scope": "First actual KEP/EC ALE native coupling group; not full candidate adoption",
        "binary": str(args.binary), "binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        "cases": cases, "cpu_sets": cpu_sets, "OMP_NUM_THREADS": 8,
        "maximum_parallel_cases": 2, "wave_CFL": .4, "cone_safety": .9,
        "mu": "immutable birth_mass/birth_volume^(2/d)",
        "weight_solver": "one scalar SPD solve, PCG cap64, relative1e-12, true residual1e-12",
        "core": "stored cell-mean KEP/EC ALE; no reconstructed face states in the core",
        "stress": "rho_bar*(c_bar*shat+2*shat^2), shat=hbar/r*max(0,-minmod(actual_delta,reconstructed_delta))",
        "stress_dt": ".4*min mass/sum(A*Pi/(-delta))",
        "step_halving_cap": 8, "minimum_dt": 1e-14,
        "conservation_limit": 1e-10, "mass_weight_residual_limit": 1e-9,
        "cold_pressure_relative_L1_limit": .01, "floor_steering_reservoir_reset": False,
        "checkpoint_version": 2,
        "still_required": "Original blast/collision and remaining 1D tests, higher-resolution/3D shock, full 2D suite, actual GPU/MPI evolution and matched GIZMO timing"})
    slots = Queue()
    for cpus in cpu_sets:
        slots.put(cpus)
    def run(specification):
        cpus = slots.get()
        try:
            record = execute(args.root, args.binary, specification, cpus)
            case = args.root/specification[0]
            if (case/"kinetic_ledger.json").exists():
                ledger = json.loads((case/"kinetic_ledger.json").read_text())
                assert ledger["core"] == "cell_mean_KEP_EC_ALE_compression_stress"
                record["kinetic_ledger"] = ledger
            if specification[3] == "cold" and (case/"final.dat").exists():
                final = np.loadtxt(case/"final.dat")
                norm = float(np.sum(final[:, 10]*abs(final[:, 12]/1e-6-1)))
                record["cold_pressure_relative_L1"] = norm
                record["cold_pressure_gate_pass"] = record["endpoint_pass"] and norm <= .01
            save(case/"group_assessment.json", record)
            return record
        finally:
            slots.put(cpus)
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(run, cases))
    save(args.root/"processes.json", results)
    print(json.dumps(results, indent=2, allow_nan=False))
    if not all(r["exit_code"] == 0 and r["endpoint_pass"] and r.get("positive_final_state", False)
               and r.get("cold_pressure_gate_pass", True) for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
