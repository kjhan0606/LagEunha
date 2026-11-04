"""Schedule the new native FV branch; no Python fluid evolution.

The first suite checks coupling and complete endpoints. Periodic translating
contact is a separate boundary control; it is not relabelled as the original
open-boundary moving-contact test. The remaining original 1D suite is pending.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
from queue import Queue
import shutil
import subprocess
import time

import numpy as np


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def execute(root, binary, specification, cpus=None, extra_env=None, restart=None):
    name, dim, side, problem, target, *control = specification
    case = root/name
    case.mkdir(exist_ok=False)
    command = [str(binary), str(dim), str(side), problem, str(target)]
    if control:
        command.append(str(control[0]))
    if restart is not None:
        if not control:
            command.append("1")
        command.append(str(restart))
    if cpus is not None:
        command = [shutil.which("taskset"), "-c", ",".join(map(str, cpus))]+command
    save(case/"command.json", command)
    env = dict(os.environ, TMPDIR=str(root), OMP_NUM_THREADS="8", OMP_DYNAMIC="FALSE",
               OMP_WAIT_POLICY="ACTIVE", OMP_PLACES="cores", OMP_PROC_BIND="spread")
    if extra_env:
        env.update(extra_env)
    save(case/"resources.json", {"CPU_set": cpus, "OMP_NUM_THREADS": 8,
                                 "explicit_environment": extra_env or {}})
    start = time.monotonic()
    with (case/"stdout.log").open("x") as out, (case/"stderr.log").open("x") as err:
        process = subprocess.Popen(command, cwd=case, env=env, stdout=out, stderr=err)
        save(case/"running.json", {"pid": process.pid, "command": command})
        code = process.wait()
    record = {"case": name, "exit_code": code, "wall_seconds": time.monotonic()-start}
    save(case/"process.json", record)
    return assess(case, specification, record)


def assess(case, specification, record):
    """Read saved output; also recovers analysis after a controller error."""
    name, dim, side, problem, target, *control = specification
    code = record["exit_code"]
    if not (case/"result.json").exists():
        return dict(record, endpoint_pass=False, reason="Native result missing; inspect stderr")
    result = json.loads((case/"result.json").read_text())
    initial = np.loadtxt(case/"initial.dat")
    final = np.loadtxt(case/"final.dat")
    endpoint = bool(code == 0 and not result["failed"] and abs(result["time"]-target) < 1e-13)
    positive = np.isfinite(final).all() and bool(np.all(final[:, [5, 10, 11, 12]] > 0))
    metrics = {"endpoint_pass": endpoint, "positive_final_state": positive,
               "new_steps": result["new_steps"], "max_ledger": result["maximum_ledger_error"],
               "max_mass_weight_with_residual": result["maximum_integrated_residual_invariant_error"],
               "retries": result["retries"], "maximum_Newton": result["maximum_Newton"],
               "maximum_CG_total_per_rhs": result["maximum_CG_total_per_rhs"]}
    if problem in ("contact", "moving_contact"):
        metrics["contact_pressure_max_error"] = float(np.max(abs(final[:, 12]-1)))
        exact_velocity = 1 if problem == "moving_contact" else 0
        metrics["contact_velocity_max_error"] = float(np.max(abs(final[:, 13]-exact_velocity)))
        metrics["contact_transverse_velocity_max"] = float(np.max(abs(final[:, 14:16])))
        metrics["contact_density_by_ID_max_change"] = float(np.max(abs(final[:, 11]-initial[:, 11])))
        metrics["contact_noise_pass"] = max(metrics["contact_pressure_max_error"],
                                              metrics["contact_velocity_max_error"],
                                              metrics["contact_transverse_velocity_max"]) < 1e-9
    history = list(csv.DictReader((case/"history.csv").open()))
    if history:
        metrics["last_history"] = history[-1]
    save(case/"assessment.json", dict(record, **metrics))
    return dict(record, **metrics)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--phase", choices=("coupling", "sod_large", "smooth"), default="coupling")
    parser.add_argument("--cpu-sets", nargs=2, metavar=("SLOT1", "SLOT2"),
                        help="Two disjoint comma-separated sets of eight allowed CPUs")
    args = parser.parse_args()
    cpu_sets = [None, None]
    if args.cpu_sets:
        cpu_sets = [sorted({int(x) for x in entry.split(",")}) for entry in args.cpu_sets]
        allowed = os.sched_getaffinity(0)
        if not shutil.which("taskset") or any(len(x) != 8 or not set(x) <= allowed for x in cpu_sets):
            parser.error("taskset and two sets of eight available CPUs are required")
        if set(cpu_sets[0]) & set(cpu_sets[1]):
            parser.error("The two CPU sets must be disjoint")
    if not args.root.resolve().is_relative_to(Path("/gpfs/kjhan/LagForce")):
        raise ValueError("New calculation directories must be under the authorized GPFS project root")
    args.root.mkdir(exist_ok=False)
    cases = [("D1_contact_N200", 1, 200, "contact", .25),
             ("D1_periodic_moving_contact_N200", 1, 200, "moving_contact", .25),
             ("D1_sod_N200", 1, 200, "sod", .2),
             ("D1_sod_N400", 1, 400, "sod", .2),
             ("D2_contact_N8", 2, 8, "contact", .25),
             ("D2_periodic_moving_contact_N8", 2, 8, "moving_contact", .25),
             ("D3_contact_N4", 3, 4, "contact", .25),
             ("D3_periodic_moving_contact_N4", 3, 4, "moving_contact", .25),
             ("D2_sod_N8", 2, 8, "sod", .2),
             ("D3_sod_N8", 3, 8, "sod", .2)]
    if args.phase == "sod_large":
        cases = [("D2_sod_N64", 2, 64, "sod", .2),
                 ("D2_sod_N128", 2, 128, "sod", .2),
                 ("D3_sod_N16", 3, 16, "sod", .2)]
    if args.phase == "smooth":
        cases = [("D2_gresho_N64", 2, 64, "gresho", 1.),
                 ("D2_gresho_N128", 2, 128, "gresho", 1.),
                 ("D2_gresho_N64_halfdt", 2, 64, "gresho", 1., .5),
                 ("D2_cold_N64", 2, 64, "cold", .25),
                 ("D2_cold_N128", 2, 128, "cold", .25)]
    registration = {"scope": "Actual original geometry and native FV coupling/evolution endpoints, not full candidate adoption",
                    "binary": str(args.binary), "binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                    "cases": cases, "phase": args.phase, "CPU_threads_per_case": 8, "maximum_parallel_cases": 2,
                    "cpu_sets": cpu_sets,
                    "wave_CFL": .4, "cone_safety": .9,
                    "reconstruction": "Centroid face-angular LSQ; scalar midpoint-extrema limiter, actual-point positivity slopes (0.9 safety), covariant midpoint velocity-ball bound with one common alpha",
                    "weight_mu": "birth_mass / birth_volume^(2/d), immutable",
                    "root_relative_tolerance": 1e-12, "Newton_cap": 8, "PCG_cap_per_Newton": 64,
                    "PCG_relative_tolerance": 1e-12, "Newton_backtrack_cap": 16,
                    "SSPRK_step_halving_cap": 8, "minimum_dt": 1e-14,
                    "conservation_gate": 1e-10, "residual_corrected_mass_weight_gate": 1e-9,
                    "contact_P_velocity_gate": 1e-9,
                    "floor_steering_energy_repair": False,
                    "missing_scope": "Original open-boundary 1D problems, full 2D/3D resolutions, temporal/GCL convergence, physical reference norms and matched GIZMO speed remain unqualified."}
    save(args.root/"registration.json", registration)
    slots = Queue()
    for cpus in cpu_sets:
        slots.put(cpus)
    def execute_on_slot(case):
        cpus = slots.get()
        try:
            return execute(args.root, args.binary, case, cpus)
        finally:
            slots.put(cpus)
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(execute_on_slot, cases))
    save(args.root/"processes.json", results)
    print(json.dumps(results, indent=2))
    if not all(r["exit_code"] == 0 and r["endpoint_pass"] and r.get("positive_final_state", False)
               for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
