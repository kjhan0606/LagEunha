"""Verify the A10 synthetic and full power-face Jacobi16 GPU probes."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


DOCS = Path(__file__).resolve().parent
SOURCES = (
    "kepler_gpu_jacobi_probe_20261027.cu",
    "kepler_gpu_jacobi_probe_smoke_20261027.cu",
    "kepler_gpu_jacobi_full_check_20261027.cu",
    "export_kepler_q2_weight_rate_csr_t65_20261027.py",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    run = args.run.resolve()
    subprocess.run(["sha256sum", "-c", "inputs.sha256"], cwd=run,
                   check=True, capture_output=True, text=True)
    smoke = (run / "smoke.log").read_text()
    if not ("both_active status=0 unsupported=0" in smoke
            and "one_active status=0 unsupported=0" in smoke
            and "unsupported_wall status=-2" in smoke):
        raise ValueError("Synthetic solve/boundary checks failed")
    rows = {}
    pattern = re.compile(r"case=(\S+) active=(\d+) faces=(\d+) status=(-?\d+) "
                         r"unsupported=(\d+) nonfinite=(\d+) elapsed_seconds=([\deE+.-]+) "
                         r"max_abs=([\deE+.-]+) relative_l2=([\deE+.-]+)")
    for label in ("aligned", "jitter"):
        content = (run / f"{label}.log").read_text().strip()
        match = pattern.fullmatch(content)
        if not match or match.group(1) != f"{label}.bin":
            raise ValueError(f"Unparseable {label} full-face result")
        values = match.groups()
        row = {"active": int(values[1]), "directed_faces": int(values[2]),
               "status": int(values[3]), "unsupported": int(values[4]),
               "nonfinite": int(values[5]), "elapsed_seconds": float(values[6]),
               "max_abs_rate_error": float(values[7]),
               "relative_l2_rate_error": float(values[8])}
        if (row["active"] != 49038 or row["directed_faces"] != 297114
                or row["status"] or row["unsupported"] or row["nonfinite"]
                or row["max_abs_rate_error"] >= 1e-10
                or row["relative_l2_rate_error"] >= 1e-9):
            raise ValueError(f"{label}: full-face CUDA/CPU gate failed")
        manifest = json.loads((run.parent / "gpu_jacobi_full_h200_v1" /
                               f"{label}.json").read_text())
        if digest(run / f"{label}.bin") != manifest["sha256"]:
            raise ValueError(f"{label}: exported face input changed")
        row["input_sha256"] = manifest["sha256"]
        row["log_sha256"] = digest(run / f"{label}.log")
        rows[label] = row
    result = {
        "job_id": 437363,
        "run": str(run),
        "smoke_log_sha256": digest(run / "smoke.log"),
        "input_sha256": {name: digest(run / name) for name in
                         ("gfs_volume_rate_probe_smoke_a10",
                          "gfs_volume_rate_probe_full_check_a10",
                          "aligned.bin", "jitter.bin", "job.slurm")},
        "source_sha256": {name: digest(DOCS / name) for name in SOURCES},
        "checks": {"staged_inputs_verified": True,
                   "synthetic_two_cell": True,
                   "dirichlet_inactive_neighbor": True,
                   "unsupported_wall_rejected": True,
                   "aligned_full_faces": True,
                   "jitter_full_faces": True},
        "cases": rows,
        "scope": "Independent H200/A10-compatible diagnostic CUDA kernel on CPU-exported t=6.5 power FaceCSR; no native GPU tessellation attachment, hydro evolution, conservation or end-to-end RK4 speed yet.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"cases": rows, "checks": result["checks"]}, indent=2))


if __name__ == "__main__":
    main()
