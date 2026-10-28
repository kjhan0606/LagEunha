"""Compare short native q1 polar variants against their bitwise baseline."""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

SNAPSHOT_MODULE = Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
sys.path.insert(0, str(SNAPSHOT_MODULE))
from audit_snapshot import read_snapshot

FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den",
          "pressure", "w2")
GPU_PATTERN = re.compile(r"\[GPU step \d+\].*SoA=([0-9.]+)ms.*total=([0-9.]+)ms")


def field_equal_exact(a: np.ndarray, b: np.ndarray) -> bool:
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def summary(log_path: Path) -> dict:
    content = log_path.read_text(errors="replace")
    pairs = np.array([(float(a), float(b)) for a, b in GPU_PATTERN.findall(content)])
    if len(pairs) == 0:
        raise ValueError(f"No GPU step timing in {log_path}")
    return {"n_steps": int(len(pairs)),
            "soa_median_ms": float(np.median(pairs[:, 0])),
            "gpu_wrapper_total_median_ms": float(np.median(pairs[:, 1])),
            "local_periodic_sync_used": "[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh" in content,
            "host_halo_sync_used": "[HALO_SYNC] rank=0 refreshed" in content,
            "run_success": "EXIT:0" in content,
            "log": str(log_path)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline", type=Path)
    ap.add_argument("candidate", type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--expected-refresh", choices=("local", "host"),
                    default="local")
    ap.add_argument("--require-wall-profile", action="store_true")
    args = ap.parse_args()
    result = {"baseline": str(args.baseline), "candidate": str(args.candidate),
              "expected_refresh": args.expected_refresh,
              "snapshots": {}, "timing": {}}
    for suffix in ("000000", "000014"):
        name = f"khout.{suffix}.dat"
        a, ma = read_snapshot(args.baseline / name)
        b, mb = read_snapshot(args.candidate / name)
        ids_a = (a["u4if"].astype(np.uint64) & ((1 << 56) - 1)).astype(np.int64)
        ids_b = (b["u4if"].astype(np.uint64) & ((1 << 56) - 1)).astype(np.int64)
        ia, ib = np.argsort(ids_a), np.argsort(ids_b)
        fields = {}
        for field in FIELDS:
            aa, bb = np.asarray(a[field][ia]), np.asarray(b[field][ib])
            fields[field] = {"same_bits": field_equal_exact(aa, bb),
                             "max_abs_difference": float(np.max(abs(aa.astype(float)-bb.astype(float))))}
        result["snapshots"][name] = {
            "same_time": ma["t"] == mb["t"],
            "same_ids": bool(np.array_equal(ids_a[ia], ids_b[ib])),
            "fields": fields,
        }
    base_timing = summary(args.baseline / "log")
    cand_timing = summary(args.candidate / "log")
    result["timing"] = {"baseline": base_timing, "candidate": cand_timing,
                        "soa_median_ratio_baseline_over_candidate":
                        base_timing["soa_median_ms"] / cand_timing["soa_median_ms"]}
    profile_log = (args.candidate / "log").read_text(errors="replace")
    wall_pattern = re.compile(
        r"\[KH_WALL\] step=\d+ t=[0-9.eE+-]+ evol=([0-9.]+)s "
        r"outside=([0-9.]+)s total=([0-9.]+)s")
    wall_rows = np.array([(float(a), float(b), float(c))
                          for a, b, c in wall_pattern.findall(profile_log)])
    if len(wall_rows):
        result["timing"]["kh_wall_profile"] = {
            "n_steps": int(len(wall_rows)),
            "evolution_median_s": float(np.median(wall_rows[:, 0])),
            "outside_median_s": float(np.median(wall_rows[:, 1])),
            "whole_iteration_median_s": float(np.median(wall_rows[:, 2]))}
    refresh_ok = (cand_timing["local_periodic_sync_used"] and
                  not cand_timing["host_halo_sync_used"])
    if args.expected_refresh == "host":
        refresh_ok = (cand_timing["host_halo_sync_used"] and
                      not cand_timing["local_periodic_sync_used"])
    result["checks"] = {
        "candidate_run_success": cand_timing["run_success"],
        "candidate_used_expected_refresh": bool(refresh_ok),
        "wall_profile_present_if_required":
        (not args.require_wall_profile or len(wall_rows) == cand_timing["n_steps"]),
        "same_step_count": base_timing["n_steps"] == cand_timing["n_steps"],
        "all_snapshot_fields_bitwise_equal": all(
            entry["same_time"] and entry["same_ids"] and
            all(f["same_bits"] for f in entry["fields"].values())
            for entry in result["snapshots"].values()),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not all(result["checks"].values()):
        raise SystemExit("Local periodic refresh A/B did not pass")


if __name__ == "__main__":
    main()
