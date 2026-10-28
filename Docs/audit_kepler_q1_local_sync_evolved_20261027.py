"""Compare matched q1 halo12 host/device refresh from accepted t=.618764."""

import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, "/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020")
from audit_snapshot import read_snapshot

from prepare_kepler_polar_long_restart_20261027 import digest, snapshot_time


FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")
GPU_PATTERN = re.compile(r"\[GPU step \d+\].*SoA=([0-9.]+)ms.*total=([0-9.]+)ms")
CPU_PATTERN = re.compile(r"\[CPU step \d+\].*total=([0-9.]+)s")
N = 147456
SOURCE_T = .6187644470884642
TARGET_T = .65


def ordered(path):
    p, meta = read_snapshot(path)
    ids = (p["u4if"].astype(np.uint64) & np.uint64((1 << 56)-1)).astype(np.int64)
    if len(p) != N or not np.array_equal(np.sort(ids), np.arange(N)):
        raise ValueError(f"Incomplete q1 native snapshot IDs: {path}")
    return p[np.argsort(ids)], meta


def one_run(run, local):
    plan = json.loads((run / "recovery_plan.json").read_text())
    audit = json.loads((run / "audit_long_t065.json").read_text())
    source = Path(plan["source"])
    if (abs(plan["source_time"]-SOURCE_T) > 1e-12 or
            abs(plan["target_time"]-TARGET_T) > 1e-12 or
            plan["source_step"] != 8119 or plan["new_halo_cells"] != 12 or
            abs(snapshot_time(source)-SOURCE_T) > 1e-12 or
            plan["inputs_sha256"].get(source.name) != digest(run / source.name) or
            not all(audit["checks"].values()) or
            abs(audit["time"]-TARGET_T) > 1e-12):
        raise ValueError(f"Invalid matched evolved checkpoint lineage: {run}")
    endpoint = Path(audit["endpoint"])
    if endpoint.parent.resolve() != run.resolve() or abs(snapshot_time(endpoint)-TARGET_T) > 1e-12:
        raise ValueError(f"Invalid exact t=.65 endpoint: {run}")
    text = (run / "log").read_text(errors="replace")
    pairs = np.asarray([(float(a), float(b)) for a, b in GPU_PATTERN.findall(text)])
    cpu = np.asarray([float(x) for x in CPU_PATTERN.findall(text)])
    if len(pairs) < 100 or len(cpu) < 100 or "EXIT:0" not in text:
        raise ValueError(f"Incomplete measured evolution: {run}")
    local_seen = "[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh" in text
    host_seen = "[HALO_SYNC] rank=0 refreshed" in text
    if (local_seen, host_seen) != ((True, False) if local else (False, True)):
        raise ValueError(f"Wrong GPU periodic refresh path: {run}")
    return {"run": str(run), "source": str(source),
            "source_sha256": plan["inputs_sha256"][source.name],
            "endpoint": str(endpoint), "endpoint_sha256": digest(endpoint),
            "local_refresh": local_seen, "host_refresh": host_seen,
            "gpu_samples": int(len(pairs)), "cpu_samples": int(len(cpu)),
            "gpu_soa_median_ms": float(np.median(pairs[5:, 0])),
            "gpu_wrapper_median_ms": float(np.median(pairs[5:, 1])),
            "cpu_reported_step_median_s": float(np.median(cpu[5:])),
            "audit_checks": audit["checks"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("host", type=Path)
    ap.add_argument("device", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    host = one_run(args.host, False)
    device = one_run(args.device, True)
    a, ma = ordered(Path(host["endpoint"]))
    b, mb = ordered(Path(device["endpoint"]))
    fields = {}
    for field in FIELDS:
        aa = np.asarray(a[field])
        bb = np.asarray(b[field])
        fields[field] = {"same_bits": aa.dtype == bb.dtype and aa.tobytes() == bb.tobytes(),
                         "max_abs_difference": float(np.max(np.abs(aa.astype(float)-bb.astype(float))))}
    checks = {
        "same_accepted_source": host["source_sha256"] == device["source_sha256"],
        "same_exact_endpoint_time": abs(ma["t"]-TARGET_T) < 1e-12 and ma["t"] == mb["t"],
        "same_step_count": host["gpu_samples"] == device["gpu_samples"],
        "all_fields_bitwise_equal": all(row["same_bits"] for row in fields.values()),
    }
    result = {"host": host, "device": device, "fields": fields,
              "checks": checks,
              "gpu_wrapper_median_ratio_host_over_device":
                  host["gpu_wrapper_median_ms"]/device["gpu_wrapper_median_ms"],
              "scope": "Same q1 accepted t=.618764 source, halo12, exact t=.65. GPU wrapper timing is not full Slurm wall time."}
    if args.output.exists():
        raise FileExistsError(args.output)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("Evolved q1 GPU refresh A/B failed")


if __name__ == "__main__":
    main()
