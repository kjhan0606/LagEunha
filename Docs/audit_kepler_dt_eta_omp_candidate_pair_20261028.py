"""Compare matched eta=0.06 RK4 runs before and after the OMP/cache candidate."""

import argparse
import json
import re
import struct
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot


FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den",
          "pressure", "w2")
GPU_STEP = re.compile(r"\[GPU step (\d+)\].*total=([0-9.]+)ms")
CPU_STEP = re.compile(r"\[CPU step (\d+)\].*total=([0-9.]+)s")
RUN_WALL = re.compile(r"\[RUN_TIME\] wall=([0-9.]+)")


def endpoint(run: Path) -> tuple[Path, dict]:
    receipt = json.loads((run / "audit_t1p1.json").read_text())
    path = Path(receipt["endpoint"])
    if not path.is_absolute():
        path = run / path
    if path.parent.resolve() != run.resolve():
        raise ValueError(f"Endpoint is outside its run directory: {path}")
    return path, receipt


def ordered_snapshot(path: Path) -> tuple[np.ndarray, dict]:
    with path.open("rb") as stream:
        magic, count, stride, stress = struct.unpack("<4i", stream.read(16))
    if (magic, count, stride, stress) != (-1, 147456, 576, 1):
        raise ValueError(f"Unexpected native snapshot header: {path}")
    if path.stat().st_size != 16 + count * stride + 52:
        raise ValueError(f"Incomplete native snapshot: {path}")
    particles, meta = read_snapshot(path)
    ids = (particles["u4if"].astype(np.uint64) & np.uint64((1 << 56) - 1))
    order = np.argsort(ids)
    if not np.array_equal(ids[order], np.arange(count, dtype=np.uint64)):
        raise ValueError(f"Invalid cell ID set: {path}")
    return particles[order], meta


def exact_field(a: np.ndarray, b: np.ndarray) -> dict:
    aa, bb = np.asarray(a), np.asarray(b)
    if aa.shape != bb.shape:
        raise ValueError("Matched endpoint field shapes differ")
    delta = np.abs(aa.astype(np.longdouble) - bb.astype(np.longdouble))
    scale = np.abs(aa.astype(np.longdouble))
    relative = np.zeros_like(delta)
    nonzero = scale > 0
    relative[nonzero] = delta[nonzero] / scale[nonzero]
    relative[~nonzero] = np.where(delta[~nonzero] > 0, np.inf, 0)
    max_relative = np.max(relative)
    return {
        "same_bits": aa.dtype == bb.dtype and aa.tobytes() == bb.tobytes(),
        "different_values": int(np.count_nonzero(delta)),
        "max_abs_difference": float(np.max(delta)),
        "max_relative_difference": (float(max_relative)
                                     if np.isfinite(max_relative) else None),
        "nonzero_difference_from_zero_reference": bool(
            np.any((scale == 0) & (delta > 0))),
    }


def timing(run: Path) -> dict:
    log = (run / "log").read_text(errors="replace")
    gpu = [(int(step), float(ms)) for step, ms in GPU_STEP.findall(log)]
    cpu = [(int(step), float(sec)) for step, sec in CPU_STEP.findall(log)]
    if not gpu or not cpu:
        raise ValueError(f"Missing per-step timings in {run / 'log'}")
    slurm_logs = sorted(run.glob("slurm-*.out"))
    if len(slurm_logs) != 1:
        raise ValueError(f"Expected one Slurm stdout file in {run}, found {len(slurm_logs)}")
    wall_values = RUN_WALL.findall(slurm_logs[0].read_text(errors="replace"))
    if len(wall_values) != 1:
        raise ValueError(f"Expected one timed mpirun record in {slurm_logs[0]}")
    return {
        "exit_zero": "EXIT:0" in log,
        "mpirun_wall_seconds": float(wall_values[0]),
        "gpu_samples": len(gpu),
        "gpu_stage_total_median_ms": float(np.median([row[1] for row in gpu])),
        "cpu_samples": len(cpu),
        "cpu_step_total_median_s": float(np.median([row[1] for row in cpu])),
    }


def passed(receipt: dict) -> bool:
    checks = receipt.get("checks", {})
    return bool(checks) and all(value is True for value in checks.values())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    base_path, base_receipt = endpoint(args.baseline)
    cand_path, cand_receipt = endpoint(args.candidate)
    base, base_meta = ordered_snapshot(base_path)
    cand, cand_meta = ordered_snapshot(cand_path)
    if base_meta["t"] != cand_meta["t"]:
        raise ValueError(f"Endpoint times differ: {base_meta['t']} vs {cand_meta['t']}")
    if base_meta["t"] != 1.1:
        raise ValueError(f"Unexpected endpoint time: {base_meta['t']}")

    fields = {name: exact_field(base[name], cand[name]) for name in FIELDS}
    base_timing = timing(args.baseline)
    cand_timing = timing(args.candidate)
    base_radial = json.loads((args.baseline / "audit_radial_t1p1.json").read_text())
    cand_radial = json.loads((args.candidate / "audit_radial_t1p1.json").read_text())
    base_radial_ok = passed(base_radial)
    cand_radial_ok = passed(cand_radial)
    exact_fields = all(item["same_bits"] for item in fields.values())
    same_gpu_samples = base_timing["gpu_samples"] == cand_timing["gpu_samples"]
    same_cpu_samples = base_timing["cpu_samples"] == cand_timing["cpu_samples"]
    mpirun_speedup = (base_timing["mpirun_wall_seconds"] /
                      cand_timing["mpirun_wall_seconds"])
    instrumented_cpu_ratio = (base_timing["cpu_step_total_median_s"] /
                              cand_timing["cpu_step_total_median_s"])
    result = {
        "baseline": str(args.baseline.resolve()),
        "candidate": str(args.candidate.resolve()),
        "endpoint_time": float(base_meta["t"]),
        "endpoint_files": {"baseline": str(base_path), "candidate": str(cand_path)},
        "per_field": fields,
        "timing": {
            "baseline": base_timing,
            "candidate": cand_timing,
            "mpirun_speedup_baseline_over_candidate": mpirun_speedup,
            "instrumented_cpu_phase_sum_ratio_baseline_over_candidate": instrumented_cpu_ratio,
            "gpu_stage_total_ratio_baseline_over_candidate": (
                base_timing["gpu_stage_total_median_ms"] /
                cand_timing["gpu_stage_total_median_ms"]),
        },
        "checks": {
            "both_global_physics_gates_pass": passed(base_receipt) and passed(cand_receipt),
            "both_radial_physics_gates_pass": base_radial_ok and cand_radial_ok,
            "both_radial_times_match_endpoint": (
                base_radial.get("time") == base_meta["t"] and
                cand_radial.get("time") == cand_meta["t"]),
            "both_runs_exit_zero": base_timing["exit_zero"] and cand_timing["exit_zero"],
            "same_gpu_step_sample_count": same_gpu_samples,
            "same_cpu_step_sample_count": same_cpu_samples,
            "all_endpoint_fields_bitwise_equal": exact_fields,
            "candidate_reduces_full_mpirun_wall_time": mpirun_speedup > 1.0,
        },
        "scope": "Matched 384^2 q1 restart and eta=0.06 t=1.1 endpoint; speed gate uses complete mpirun wall time, not the partial CPU phase timer; this does not establish long-time stability.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    if not all(result["checks"].values()):
        raise SystemExit("Matched eta=0.06 OMP/cache comparison failed")


if __name__ == "__main__":
    main()
