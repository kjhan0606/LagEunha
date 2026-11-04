"""Fixed grouped gates for the LSQ compression/residual stress closure."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import traceback
from pathlib import Path
from queue import Queue

import numpy as np

from run_fv_native_suite import execute, save


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--select", nargs="+", choices=("alternating", "sod", "blast", "collision", "cold", "gresho"))
    parser.add_argument("--restart-blast", type=Path)
    parser.add_argument("--core-version", type=int, choices=(3, 4), default=3)
    args = parser.parse_args()
    if not args.root.resolve().is_relative_to(Path("/gpfs/kjhan/LagForce")):
        raise ValueError("Use GPFS")
    args.root.mkdir(exist_ok=False)
    cases = [("D1_alternating_N16", 1, 16, "alternating", float(1/16/np.sqrt(1.4))),
             ("D1_sod_N200", 1, 200, "sod", .2),
             ("D1_blast_N200", 1, 200, "blast", .038),
             ("D1_collision_N200", 1, 200, "collision", .035),
             ("D2_cold_N64", 2, 64, "cold", .25),
             ("D2_gresho_N64", 2, 64, "gresho", 1.)]
    full_matrix = list(cases)
    if args.select:
        cases = [case for case in cases if case[3] in args.select]
    save(args.root/"registration.json", {
        "core": "cell-mean KEP/EC ALE; checkpoint version"+str(args.core_version),
        "binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(), "cases": cases,
        "full_matrix": full_matrix, "selected_subset": args.select,
        "blast_restart": str(args.restart_blast) if args.restart_blast else None,
        "detector": "chi=(Eperp+c*Eparallel)/(Eperp+Eparallel), c=min(1,max(0,-trG)^2/Frobenius(G)^2); original unlimited centroid LSQ weights A/r^(d+1)",
        "pair_fraction": "max(chi_i,chi_j)", "unchanged_stress_coefficients": [1, 2],
        "wave_CFL": .4, "cone_safety": .9, "PCG_cap": 64, "true_root_tolerance": 1e-12,
        "step_halving_cap": 8, "minimum_dt": 1e-14, "GFS_TEST_GUARD_MULTIPLIER": 1,
        "CPU_sets": [list(range(8)), list(range(8, 16))], "threads_each": 8,
        "cold_pressure_relative_L1_limit": .01,
        "Gresho64_physical_node_error_must_improve": .036911931828247264,
        "alternating_final_kinetic_fraction_limit": .99,
        "strong_failure_control": "If a strong shock fails, run chi=1 once for the first failed blast/collision."
        if args.core_version == 3 else "If PCG fails, test mutual-strongest-edge 2x2 block Jacobi on the exact failed trial matrix before changing the common numerical backend; keep cap/tolerance/mu.",
        "velocity_query": "Own-cell extrapolation" if args.core_version == 3 else "Actual containing power cell, native relative halfspaces/adjacency, deterministic ties, explicit failure, unchanged own-cell value on inside points",
        "not_yet_qualified": "Original nine full tests/refinement, guard-doubling equivalence, higher-resolution/3D suite, native GPU/MPI and matched GIZMO timing"})
    slots = Queue()
    for cpus in [list(range(8)), list(range(8, 16))]:
        slots.put(cpus)
    def run(specification):
        cpus = slots.get()
        try:
            record = execute(args.root, args.binary, specification, cpus,
                             {"GFS_TEST_GUARD_MULTIPLIER": "1"},
                             args.restart_blast if specification[3] == "blast" else None)
            case = args.root/specification[0]
            if (case/"kinetic_ledger.json").exists():
                record["kinetic_ledger"] = json.loads((case/"kinetic_ledger.json").read_text())
                assert record["kinetic_ledger"]["checkpoint_version"] == args.core_version
            if (case/"final.dat").exists():
                final = np.loadtxt(case/"final.dat")
                if specification[3] == "cold":
                    record["cold_pressure_relative_L1"] = float(np.sum(final[:, 10]*abs(final[:, 12]/1e-6-1)))
                    record["physical_gate_pass"] = record["endpoint_pass"] and record["cold_pressure_relative_L1"] <= .01
                if specification[3] == "alternating":
                    birth = np.loadtxt(case/"initial.dat")
                    kinetic = lambda a: float(np.sum(.5*np.sum(a[:, 6:9]**2, axis=1)/a[:, 5]))
                    record["kinetic_fraction"] = kinetic(final)/kinetic(birth)
                    record["physical_gate_pass"] = record["endpoint_pass"] and record["kinetic_fraction"] < .99
                if specification[3] == "gresho" and record["endpoint_pass"]:
                    record["physical_node_velocity_L1"] = float(record["last_history"]["physical_velocity_node_L1_from_initial_profile"])
                    record["physical_gate_pass"] = record["physical_node_velocity_L1"] < .036911931828247264
            save(case/"gate_assessment.json", record)
            return record
        except Exception:
            case = args.root/specification[0]
            record = {"case": specification[0], "endpoint_pass": False,
                      "exit_code": json.loads((case/"process.json").read_text())["exit_code"]
                      if (case/"process.json").exists() else None,
                      "controller_error": traceback.format_exc()}
            save(case/"controller_error.json", record)
            return record
        finally:
            slots.put(cpus)
    with ThreadPoolExecutor(max_workers=2) as executor:
        records = list(executor.map(run, cases))
    save(args.root/"processes.json", records)
    print(json.dumps(records, indent=2, allow_nan=False))
    if not all(r["exit_code"] == 0 and r["endpoint_pass"] and r.get("physical_gate_pass", True) for r in records):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
