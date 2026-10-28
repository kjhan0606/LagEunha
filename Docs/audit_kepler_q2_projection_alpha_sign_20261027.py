"""Summarize sign changes of the native global pressure-work projection."""

import json
from pathlib import Path
import re

import numpy as np

from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
RUNS = {
    "early": (ROOT / "gpu_energy_neutral_projected_v11b_t0825_to1_ada_v1", "audit_t1.json", 1.),
    "late": (ROOT / "gpu_energy_neutral_projected_v11b_t65_to66_ada_v1", "audit_t66.json", 6.6),
}
OUTPUT = Path(__file__).resolve().parent / "KEPLER_Q2_PROJECTION_ALPHA_SIGN_EARLY_LATE_2026-10-27.json"
ALPHA = re.compile(r"\[GPU_VOLW2\] stage=(\d+) active=49038 faces=\d+ "
                   r"alpha=([\deE+.-]+) work=([\deE+.-]+) elapsed=[\deE+.-]+ s")
TIME = re.compile(r"Time=\s*([\deE+.-]+)\s+step=\s*\d+\s+dt=\s*([\deE+.-]+)")


def summarize(run, receipt_name, end_time):
    receipt_path = run / receipt_name
    receipt = json.loads(receipt_path.read_text())
    if (receipt["time"] != end_time or len(receipt["checks"]) != 21
            or not all(receipt["checks"].values())
            or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]):
        raise ValueError(f"Unaccepted native projection run: {run}")
    stages = []
    effective = []
    times = []
    for line in (run / "log").read_text(errors="replace").splitlines():
        match = ALPHA.fullmatch(line)
        if match:
            stages.append((int(match.group(1)), float(match.group(2))))
        clock = TIME.fullmatch(line)
        if clock:
            if len(stages) != 4 or [row[0] for row in stages] != list(range(4)):
                raise ValueError(f"Missing four stages before {line}")
            effective.append(float(np.dot([1., 2., 2., 1.],
                                          [row[1] for row in stages])/6.))
            times.append(float(clock.group(1)))
            stages.clear()
    if stages or not times or abs(times[-1]-end_time) > 1e-10:
        raise ValueError("Incomplete native projection sign ledger")
    alpha = np.asarray(effective)
    return {"run": str(run), "receipt_sha256": digest(receipt_path),
            "log_sha256": digest(run / "log"), "steps": len(alpha),
            "positive_effective_steps": int(np.count_nonzero(alpha > 0)),
            "negative_effective_steps": int(np.count_nonzero(alpha < 0)),
            "effective_alpha_min": float(alpha.min()),
            "effective_alpha_median": float(np.median(alpha)),
            "effective_alpha_max": float(alpha.max()),
            "first_time": times[0], "last_time": times[-1]}


def main():
    report = {name: summarize(*arguments) for name, arguments in RUNS.items()}
    report["interpretation_scope"] = (
        "With linearized local volume Jacobian J and pressure gradient h=J^T P, "
        "the projection changes volume rate by -alpha J h. For a constant-entropy "
        "cell this produces a pressure-rate term proportional to alpha J J^T P; "
        "negative alpha damps a positive pressure mode and positive alpha can "
        "amplify it. This is a local linear inference, not a measured modal "
        "growth rate or a whole-disk stability proof.")
    with OUTPUT.open("x") as out:
        json.dump(report, out, indent=2, allow_nan=False)
        out.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
