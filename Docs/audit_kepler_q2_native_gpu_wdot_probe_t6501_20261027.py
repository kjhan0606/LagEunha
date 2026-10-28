"""Check native CUDA FaceCSR Wdot against CPU geometry and W0 replay."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
ORIGINAL = ROOT / "faceforce_pair_t6501_a10_v1"


def one_case(root, original_root, label, reference):
    run = root / label
    old = original_root / label
    new_receipt = json.loads((run / "audit_t6501.json").read_text())
    old_receipt = json.loads((old / "audit_t6501.json").read_text())
    if (len(new_receipt["checks"]) != 11 or not all(new_receipt["checks"].values())
            or len(old_receipt["checks"]) != 11 or not all(old_receipt["checks"].values())
            or new_receipt["restart_sha256"] != old_receipt["restart_sha256"]
            or abs(snapshot_time(Path(new_receipt["endpoint"]), N)-6.501) > 1e-10
            or digest(Path(new_receipt["endpoint"])) != new_receipt["endpoint_sha256"]
            or digest(Path(old_receipt["endpoint"])) != old_receipt["endpoint_sha256"]):
        raise ValueError(f"{label}: replay or original comparison gate failed")
    native, _ = ordered(Path(new_receipt["endpoint"]), N)
    control, _ = ordered(Path(old_receipt["endpoint"]), N)
    fields = {}
    for name in FIELDS:
        a, b = native[name], control[name]
        fields[name] = {"bitwise_equal": bool(np.array_equal(a, b)),
                        "max_absolute_difference": float(np.max(np.abs(
                            a.astype(float)-b.astype(float))))}
    path = run / "gpu_wdot_probe_r0_call1.tsv"
    if not path.is_file():
        raise ValueError(f"{label}: missing native GPU Wdot dump")
    rows = np.loadtxt(path, skiprows=1, dtype=[("id", "<u8"), ("rate", "<f8")])
    ids = rows["id"].astype(int)
    if len(ids) != N or not np.array_equal(np.sort(ids), np.arange(N)):
        raise ValueError(f"{label}: incomplete native Wdot IDs")
    rate = rows["rate"][np.argsort(ids)]
    if not np.all(np.isfinite(rate)):
        raise ValueError(f"{label}: nonfinite native Wdot")
    mask = reference["active"] > 0
    expected = reference[label]
    residual = rate[mask]-expected[mask]
    relative_l2 = float(np.linalg.norm(residual)/np.linalg.norm(expected[mask]))
    log = (run / "log").read_text(errors="replace")
    lines = [line for line in log.splitlines() if "[GPU_WDOT_PROBE] active=" in line]
    if len(lines) != 1 or "active=49038" not in lines[0] or "jacobi=16 omega=0.8" not in lines[0]:
        raise ValueError(f"{label}: missing one native probe timing row")
    return {
        "new_endpoint_sha256": new_receipt["endpoint_sha256"],
        "original_endpoint_sha256": old_receipt["endpoint_sha256"],
        "new_short_gate": new_receipt["checks"],
        "physical_fields": fields,
        "gpu_probe_dump_sha256": digest(path),
        "gpu_probe_log": lines[0],
        "active_rate_max_abs_difference": float(np.max(np.abs(residual))),
        "active_rate_relative_l2_difference": relative_l2,
        "inactive_rate_max_abs": float(np.max(np.abs(rate[~mask]))),
        "native_rate_abs_max": float(np.max(np.abs(rate[mask]))),
        "reference_rate_abs_max": float(np.max(np.abs(expected[mask]))),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("staged_pair", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = args.staged_pair.resolve()
    manifest = json.loads((root / "stage_manifest.json").read_text())
    if manifest["binary_sha256"] != digest(root / "aligned/eunha2.exe"):
        raise ValueError("Staged native binary changed")
    if digest(root / "disk_active_mask.bin") != manifest["mask_sha256"]:
        raise ValueError("Active mask changed")
    if digest(root / "cpu_reference.npz") != manifest["reference_sha256"]:
        raise ValueError("CPU reference changed")
    with np.load(root / "cpu_reference.npz") as reference:
        cases = {label: one_case(root, ORIGINAL, label, reference)
                 for label in ("aligned", "jitter")}
    checks = {"accepted_parent_replay": all(
                  all(case["new_short_gate"].values()) for case in cases.values()),
              "no_physical_change": all(
                  all(row["bitwise_equal"] for row in case["physical_fields"].values())
                  for case in cases.values()),
              "finite_native_wdot": all(
                  np.isfinite(case["active_rate_relative_l2_difference"])
                  for case in cases.values()),
              "inactive_rate_zero": all(case["inactive_rate_max_abs"] == 0.
                                        for case in cases.values())}
    result = {"run": str(root), "binary_sha256": manifest["binary_sha256"],
              "mask_sha256": manifest["mask_sha256"],
              "reference_sha256": manifest["reference_sha256"],
              "cases": cases, "checks": checks,
              "scope": "One native W0 t=6.5->6.501 RK4 replay per IC with opt-in read-only GPU Jacobi on the actual tessellator FaceCSR. CPU comparison uses independent SciPy power geometry; no weighted hydro evolution."}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"checks": checks,
                      "rates": {label: {k: cases[label][k] for k in
                          ("active_rate_max_abs_difference", "active_rate_relative_l2_difference")}
                                for label in cases}}, indent=2))
    if not all(checks.values()):
        raise SystemExit("Native GPU Wdot replay gate failed")


if __name__ == "__main__":
    main()
