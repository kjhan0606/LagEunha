"""Verify projected-rate CUDA/CPU comparisons on accepted W0 t=6.5 faces."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


DOCS = Path(__file__).resolve().parent
PATTERN = re.compile(
    r"case=(\S+) active=(\d+) faces=(\d+) status=(-?\d+) "
    r"unsupported=(\d+) nonfinite=(\d+) elapsed_seconds=([\deE+.-]+) "
    r"max_rate=([\deE+.-]+) max_gradient=([\deE+.-]+) "
    r"max_projected=([\deE+.-]+) relative_orthogonality=([\deE+.-]+) "
    r"reference_relative_orthogonality=([\deE+.-]+) alpha=([\deE+.-]+)"
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve()
    subprocess.run(["sha256sum", "-c", "inputs.sha256"], cwd=run,
                   check=True, capture_output=True, text=True)
    rows = {}
    for label in ("aligned", "jitter"):
        manifest = json.loads((run.parent/"gpu_energy_neutral_w0_export_v1"/
                               f"{label}.json").read_text())
        match = PATTERN.fullmatch((run/f"{label}.log").read_text().strip())
        if not match or Path(match.group(1)) != Path(manifest["binary"]):
            raise ValueError(f"Unrecognized {label} GPU output")
        fields = match.groups()
        row = {
            "active": int(fields[1]), "directed_faces": int(fields[2]),
            "status": int(fields[3]), "unsupported": int(fields[4]),
            "nonfinite": int(fields[5]), "elapsed_seconds": float(fields[6]),
            "max_rate_difference": float(fields[7]),
            "max_pressure_gradient_difference": float(fields[8]),
            "max_projected_rate_difference": float(fields[9]),
            "relative_work_orthogonality": float(fields[10]),
            "reference_relative_work_orthogonality": float(fields[11]),
            "projection_alpha": float(fields[12]),
            "input_sha256": digest(Path(manifest["binary"])),
            "log_sha256": digest(run/f"{label}.log"),
        }
        if (row["input_sha256"] != manifest["sha256"] or
                row["active"] != 49038 or row["directed_faces"] != 297114 or
                row["status"] or row["unsupported"] or row["nonfinite"] or
                row["max_rate_difference"] >= 1e-10 or
                row["max_pressure_gradient_difference"] >= 1e-10 or
                row["max_projected_rate_difference"] >= 1e-10 or
                row["relative_work_orthogonality"] >= 1e-12):
            raise ValueError(f"{label} GPU projection failed")
        rows[label] = row
    result = {
        "job_id": 437374, "run": str(run),
        "binary_sha256": digest(run/"gfs_energy_neutral_full_check"),
        "source_sha256": {
            name: digest(DOCS/name) for name in (
                "kepler_gpu_jacobi_probe_20261027.cu",
                "kepler_gpu_energy_neutral_probe_validated_v1_20261027.cu",
                "kepler_gpu_energy_neutral_full_check_20261027.cu",
                "export_kepler_q2_energy_neutral_csr_t65_20261027.py",
            )},
        "cases": rows,
        "checks": {"input_hashes": True, "aligned_full_w0_faces": True,
                   "jitter_full_w0_faces": True, "work_neutral": True},
        "scope": "A10 CUDA Jacobi16 plus pressure-gradient projection on CPU-exported accepted t=6.5 W0 FaceCSR. No native weighted-force attachment or hydro evolution yet. H200 submission 437373 was cancelled while QOSGrpGRES pending.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"cases": rows, "checks": result["checks"]},
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
