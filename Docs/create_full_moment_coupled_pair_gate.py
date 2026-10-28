"""Build the paired fine/coarse gate for one accepted native coupled endpoint.

This is read-only with respect to simulation state. It checks native acceptance,
the existing physical monitor, the approved temporal-difference budget, source
identity, and the predeclared compute-wall forecast before emitting a restart
pair gate. It never advances or modifies a physical state.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

EPS = sys.float_info.epsilon
INTEGRATOR_LABEL = "self-consistent midpoint + Sylvester/exponential/positive QR"
MAX_COMBINED_WALL_SECONDS = 24 * 60 * 60
MAX_TOTAL_PROVIDER_CALLS_TO_005 = 150
EXPECTED_CELLS = 19945


class GateFailure(RuntimeError):
    pass


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        while block := f.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def close_time(a, b):
    return abs(float(a) - float(b)) <= 32 * EPS * max(1.0, abs(float(a)), abs(float(b)))


def parse_params(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("define ") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip()
    return result


def slurm_accounting(job_id):
    if not job_id.isdigit():
        raise GateFailure(f"Invalid Slurm job id: {job_id!r}")
    command = [
        "sacct", "-n", "-P", "-j", job_id,
        "--format=JobID,State,ElapsedRaw,ExitCode",
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise GateFailure(f"sacct failed for {job_id}: {result.stderr[-1200:]}")
    rows = [line.split("|") for line in result.stdout.splitlines() if line.strip()]
    parent = [row for row in rows if len(row) >= 4 and row[0].strip() == job_id]
    if len(parent) != 1:
        raise GateFailure(f"Expected one parent accounting row for {job_id}; got {parent!r}")
    state = parent[0][1].strip().split("+", 1)[0]
    if state != "COMPLETED":
        raise GateFailure(f"Slurm job {job_id} ended in {state}")
    if parent[0][3].strip() != "0:0":
        raise GateFailure(f"Slurm job {job_id} has nonzero exit code {parent[0][3]!r}")
    try:
        elapsed = int(parent[0][2].strip())
    except ValueError as exc:
        raise GateFailure(f"Missing ElapsedRaw for {job_id}: {parent[0]!r}") from exc
    if elapsed <= 0:
        raise GateFailure(f"Nonpositive Slurm elapsed time for {job_id}")
    return dict(job_id=job_id, state=state, elapsed_seconds=elapsed)


def validate_branch(runtime, job_id, step, start_step, expected_h, target_time,
                    restart_provenance, state_digest):
    runtime = Path(runtime).resolve(strict=True)
    if not close_time(start_step * expected_h, 0.025):
        raise GateFailure(f"The branch does not start from the approved t=.025 step: {start_step}")
    if not close_time(target_time, 0.03):
        raise GateFailure("This producer implements only the initial paired t=.03 gate")
    checkpoint = runtime / f"accepted_{step:04d}"
    accepted_path = runtime / f"accepted_{step:04d}.json"
    completed_path = runtime / f"completed_{step:04d}.json"
    accepted = read_json(accepted_path)
    completed = read_json(completed_path)
    if accepted["status"] != "accepted_native_coupled_midpoint_step":
        raise GateFailure(f"{accepted_path} is not an accepted coupled step")
    if completed["status"] != "completed_native_evolution_batch":
        raise GateFailure(f"{completed_path} is not a completed batch")
    for report in (accepted, completed):
        if not report["accepted_commit"] or not report["fixed_mass"] or not report["Re_infinity"]:
            raise GateFailure(f"Native acceptance invariants failed in {report}")
        if report["ranks"] != 2 or report["cells"] != EXPECTED_CELLS:
            raise GateFailure(f"Unexpected rank/cell count in {report}")
        if report["step"] != step or report["h"] != expected_h:
            raise GateFailure(f"Unexpected step or h in {report}")
        if not close_time(report["time"], target_time) or not close_time(report["time"], step * expected_h):
            raise GateFailure(f"Unexpected endpoint time in {report}")
        if report["integrator"] != INTEGRATOR_LABEL:
            raise GateFailure(f"Unexpected integrator in {report}")
    if completed["step"] != accepted["step"] or completed["h"] != accepted["h"]:
        raise GateFailure("Completed and accepted reports disagree")
    if not checkpoint.is_dir():
        raise GateFailure(f"Missing accepted checkpoint: {checkpoint}")

    steps = []
    for number in range(start_step + 1, step + 1):
        report_path = runtime / f"accepted_{number:04d}.json"
        control_path = runtime / f"coupled_{number:04d}" / "control.json"
        report = read_json(report_path)
        control = read_json(control_path)
        if report["status"] != "accepted_native_coupled_midpoint_step" or not report["accepted_commit"]:
            raise GateFailure(f"Unaccepted prefix step {number}")
        if report["h"] != expected_h or not close_time(report["time"], number * expected_h):
            raise GateFailure(f"Unexpected prefix endpoint at step {number}")
        if not math.isfinite(report["step_wall_seconds"]) or report["step_wall_seconds"] <= 0:
            raise GateFailure(f"Missing accepted-step wall time at step {number}")
        if (control["status_code"] != 0 or control["converged"] is not True
                or control["final_checks_passed"] is not True
                or control["outer_native_commit_required"] is not True
                or control["accepted"] is not False):
            raise GateFailure(f"Coupled control failed at step {number}")
        if not 0 < control["provider_calls"] <= 10:
            raise GateFailure(f"Provider cap failed at step {number}")
        audit = runtime / f"candidate_audit_{number:04d}"
        admissibility = read_json(audit / "admissibility.json")
        physical = read_json(audit / "physical_gate.json")
        density = read_json(audit / "density" / "report.json")
        agreement = read_json(runtime / f"coupled_{number:04d}" / "agreement" / "report.json")
        journal = read_json(runtime / f"coupled_{number:04d}" / "accepted_journal" / "report.json")
        if (admissibility["status"] != "passed_native_tentative_admissibility"
                or not admissibility["mass_bitwise_initial"]
                or admissibility["minimum_thermal"] <= 0):
            raise GateFailure(f"Admissibility gate failed at step {number}")
        ledger_residual = admissibility["P_E_L_ledger_residual"]
        ledger_bound = admissibility["P_E_L_ledger_bound"]
        if (not isinstance(ledger_residual, list) or not isinstance(ledger_bound, list)
                or len(ledger_residual) != 4 or len(ledger_bound) != 4):
            raise GateFailure(f"Malformed P/E/L ledger vector at step {number}")
        ledger_residual = [float(value) for value in ledger_residual]
        ledger_bound = [float(value) for value in ledger_bound]
        if (not all(math.isfinite(value) for value in ledger_residual + ledger_bound)
                or any(bound <= 0 or abs(residual) > bound
                       for residual, bound in zip(ledger_residual, ledger_bound))):
            raise GateFailure(f"P/E/L ledger failed at step {number}")
        combined_energy_change = float(admissibility["combined_energy_change"])
        combined_energy_bound = float(admissibility["combined_energy_roundoff_bound"])
        if (not math.isfinite(combined_energy_change)
                or not math.isfinite(combined_energy_bound)
                or combined_energy_bound <= 0
                or abs(combined_energy_change) > combined_energy_bound):
            raise GateFailure(f"Combined energy gate failed at step {number}")
        if (physical["status"] != "passed_checkpoint_physical_monitor"
                or not physical["regions"]
                or not all(row["passed"] for row in physical["regions"].values())):
            raise GateFailure(f"Physical monitor failed at step {number}")
        if (density["status"] != "passed_full_stage_Psi_density_preparation"
                or density["source"] != "Native-owned fits with unchanged independent q48"
                or density["cells"] != EXPECTED_CELLS
                or not close_time(density["time"], report["time"])):
            raise GateFailure(f"Independent q48/density audit failed at step {number}")
        candidate = audit / "candidate.npz"
        if sha(candidate) != density["state_sha256"]:
            raise GateFailure(f"Candidate density-audit state hash mismatch at step {number}")
        if (agreement["status"] != "passed_actual_midpoint_agreement"
                or journal["status"] != "captured_verified_map_cell_journal"
                or not close_time(journal["time"], report["time"])
                or not close_time(journal["h"], expected_h)):
            raise GateFailure(f"Midpoint/journal audit failed at step {number}")
        steps.append(dict(step=number, time=report["time"], h=expected_h,
                          step_wall_seconds=report["step_wall_seconds"],
                          provider_calls=control["provider_calls"],
                          checkpoint_path=str(runtime / f"accepted_{number:04d}")))

    identity_path, identity, operator_sha, sources, params_path, provenance_path, build_manifest_path, build = \
        restart_provenance.source_identity(checkpoint, 2)
    params = parse_params(params_path)
    if not checkpoint.exists():
        raise GateFailure(f"Missing final checkpoint: {checkpoint}")
    return dict(
        accepted=True,
        physical_gates_passed=True,
        integrator=INTEGRATOR_LABEL,
        h=expected_h,
        step=step,
        time=float(accepted["time"]),
        checkpoint_path=str(checkpoint),
        checkpoint_sha256=state_digest(checkpoint),
        job_id=job_id,
        runtime_identity_path=str(identity_path),
        runtime_identity_sha256=sha(identity_path),
        binary_sha256=identity["binary_sha256"],
        operator_sha256=operator_sha,
        provider_sources_sha256=sha(provenance_path),
        build_manifest_sha256=sha(build_manifest_path),
        user_main_sha256=build["user_main_sha256"],
        params=params,
        accepted_steps=steps,
    )


def run(args, output):
    global restart_provenance, state_digest
    sys.path.insert(0, str(args.restart_provider.resolve(strict=True)))
    import full_moment_coupled_restart_provenance as restart_provenance
    from full_moment_state_io import digest as state_digest

    fine = validate_branch(args.fine_runtime, args.fine_job, args.fine_step,
                           args.fine_start_step, 0.0025, args.time,
                           restart_provenance, state_digest)
    coarse = validate_branch(args.coarse_runtime, args.coarse_job, args.coarse_step,
                             args.coarse_start_step, 0.005, args.time,
                             restart_provenance, state_digest)
    shared_identity = (
        "binary_sha256", "operator_sha256", "provider_sources_sha256",
        "build_manifest_sha256", "user_main_sha256",
    )
    for key in shared_identity:
        if fine[key] != coarse[key]:
            raise GateFailure(f"Fine/coarse runtime identity mismatch: {key}")
    if args.expected_user_source_sha and fine["user_main_sha256"] != args.expected_user_source_sha:
        raise GateFailure("The user main-source SHA does not match the requested source")
    shared_params = (
        "moment_state", "moment_initial", "moment_physics_file", "moment_physics_sha",
        "moment_fork_manifest", "moment_integrator", "moment_workers", "moment_threads",
    )
    for key in shared_params:
        if fine["params"].get(key) != coarse["params"].get(key):
            raise GateFailure(f"Fine/coarse input parameter mismatch: {key}")
    if fine["params"].get("moment_integrator") != "self_consistent_midpoint":
        raise GateFailure("The pair did not use self-consistent midpoint")
    for label, branch in (("fine", fine), ("coarse", coarse)):
        for key in ("moment_state", "moment_initial", "moment_physics_file", "moment_fork_manifest"):
            path = Path(branch["params"][key]).resolve(strict=True)
            branch.setdefault("input_sha256", {})[key] = sha(path)
    if fine["input_sha256"] != coarse["input_sha256"]:
        raise GateFailure("Fine/coarse branch inputs are not byte-identical")

    comparison_script = Path(__file__).with_name("compare_full_moment_checkpoint_errors.py").resolve()
    budget_script = Path(__file__).with_name("evaluate_full_moment_checkpoint_difference_budget.py").resolve()
    comparison_dir = output / "comparison"
    budget_dir = output / "temporal_budget"
    fine_physics = Path(args.fine_runtime).resolve() / f"candidate_audit_{args.fine_step:04d}" / "physics"
    coarse_physics = Path(args.coarse_runtime).resolve() / f"candidate_audit_{args.coarse_step:04d}" / "physics"
    compare = subprocess.run(
        [sys.executable, str(comparison_script), "--a", str(fine_physics),
         "--b", str(coarse_physics), "--output", str(comparison_dir)],
        text=True, capture_output=True,
    )
    if compare.returncode:
        raise GateFailure(f"Checkpoint comparison failed: {compare.stderr[-1600:]}")
    comparison_path = comparison_dir / "report.json"
    comparison = read_json(comparison_path)
    if (comparison["status"] != "completed_checkpoint_actual_cell_error_difference_diagnostics"
            or not comparison["existing_physical_monitors_both_passed"]
            or not close_time(comparison["time"], args.time)):
        raise GateFailure("Checkpoint comparison did not establish the required paired monitor")
    budget = subprocess.run(
        [sys.executable, str(budget_script), "--comparison", str(comparison_path),
         "--kind", "temporal", "--output", str(budget_dir)],
        text=True, capture_output=True,
    )
    if budget.returncode:
        raise GateFailure(f"Temporal-budget analysis failed: {budget.stderr[-1600:]}")
    budget_path = budget_dir / "report.json"
    budget_report = read_json(budget_path)
    consistency = budget_report["status"] == "passed_temporal_physical_difference_budget"

    fine_job = slurm_accounting(args.fine_job)
    coarse_job = slurm_accounting(args.coarse_job)
    fine_wall = sum(row["step_wall_seconds"] for row in fine["accepted_steps"])
    coarse_wall = sum(row["step_wall_seconds"] for row in coarse["accepted_steps"])
    fine_mean = fine_wall / len(fine["accepted_steps"])
    coarse_mean = coarse_wall / len(coarse["accepted_steps"])
    fine_remaining = round((0.05 - args.time) / fine["h"])
    coarse_remaining = round((0.05 - args.time) / coarse["h"])
    if fine_remaining < 0 or coarse_remaining < 0:
        raise GateFailure("The .05 forecast is incompatible with the current endpoint")
    if (fine_remaining != 8 or coarse_remaining != 4
            or len(fine["accepted_steps"]) + len(coarse["accepted_steps"])
            + fine_remaining + coarse_remaining != 15):
        raise GateFailure("The continuation forecast does not match the approved 15-step .025-to-.05 bundle")
    actual_providers = sum(row["provider_calls"] for row in fine["accepted_steps"])
    actual_providers += sum(row["provider_calls"] for row in coarse["accepted_steps"])
    provider_forecast = actual_providers + 10 * (fine_remaining + coarse_remaining)
    if provider_forecast > MAX_TOTAL_PROVIDER_CALLS_TO_005:
        raise GateFailure("The provider-call forecast exceeds the approved .025-to-.05 budget")
    observed_job_wall = fine_job["elapsed_seconds"] + coarse_job["elapsed_seconds"]
    forecast = observed_job_wall + fine_mean * fine_remaining + coarse_mean * coarse_remaining
    compute_budget = dict(
        cap_seconds=MAX_COMBINED_WALL_SECONDS,
        observed_job_wall_seconds=observed_job_wall,
        accepted_step_wall_seconds=dict(fine=fine_wall, coarse=coarse_wall),
        mean_accepted_step_wall_seconds=dict(fine=fine_mean, coarse=coarse_mean),
        remaining_steps_to_005=dict(fine=fine_remaining, coarse=coarse_remaining),
        provider_calls_to_005=dict(
            cap=MAX_TOTAL_PROVIDER_CALLS_TO_005,
            already_used=actual_providers,
            forecast_maximum=provider_forecast,
            remaining_step_cap=10 * (fine_remaining + coarse_remaining),
            passed=provider_forecast <= MAX_TOTAL_PROVIDER_CALLS_TO_005,
        ),
        forecast_combined_wall_seconds=forecast,
        passed=forecast <= MAX_COMBINED_WALL_SECONDS,
        slurm_jobs=dict(fine=fine_job, coarse=coarse_job),
    )
    all_physical = fine["physical_gates_passed"] and coarse["physical_gates_passed"]
    status = "passed_coupled_pair_gate" if all_physical and consistency and compute_budget["passed"] else "failed_coupled_pair_gate"
    return dict(
        status=status,
        all_physical_gates_passed=all_physical,
        fine_coarse_consistency_passed=consistency,
        time=float(args.time),
        branches=dict(fine=fine, coarse=coarse),
        comparison_artifact_path=str(comparison_path.resolve()),
        comparison_artifact_sha256=sha(comparison_path),
        analysis_script_path=str(comparison_script),
        analysis_script_sha256=sha(comparison_script),
        temporal_budget_artifact_path=str(budget_path.resolve()),
        temporal_budget_artifact_sha256=sha(budget_path),
        temporal_budget_script_path=str(budget_script),
        temporal_budget_script_sha256=sha(budget_script),
        pair_gate_script_path=str(Path(__file__).resolve()),
        pair_gate_script_sha256=sha(__file__),
        temporal_budget=budget_report.get("budget"),
        compute_budget=compute_budget,
        continuation_target_time=0.05,
        automatic_extension_authorized=False,
        physical_providers=0,
        accepted_time_steps=0,
        limits=("Pair gate only. Requires both native branches, the existing per-region temporal "
                "budget, and the predeclared combined 24-hour compute forecast. It does not establish "
                "t>6 stability or authorize continuation beyond the .05 evaluation."),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fine-runtime", type=Path, required=True)
    parser.add_argument("--coarse-runtime", type=Path, required=True)
    parser.add_argument("--fine-job", required=True)
    parser.add_argument("--coarse-job", required=True)
    parser.add_argument("--fine-start-step", type=int, required=True)
    parser.add_argument("--coarse-start-step", type=int, required=True)
    parser.add_argument("--fine-step", type=int, required=True)
    parser.add_argument("--coarse-step", type=int, required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--restart-provider", type=Path, required=True)
    parser.add_argument("--expected-user-source-sha", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    try:
        report = run(args, output)
    except Exception as exc:
        report = dict(
            status="failed_coupled_pair_gate",
            all_physical_gates_passed=False,
            fine_coarse_consistency_passed=False,
            time=args.time,
            failure=repr(exc),
            physical_providers=0,
            accepted_time_steps=0,
        )
    report_path = output / "pair_gate.json"
    with report_path.open("x") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps(dict(status=report["status"], time=report["time"],
                         pair_gate=str(report_path)), flush=True))
    return 0 if report["status"] == "passed_coupled_pair_gate" else 1


if __name__ == "__main__":
    raise SystemExit(main())
