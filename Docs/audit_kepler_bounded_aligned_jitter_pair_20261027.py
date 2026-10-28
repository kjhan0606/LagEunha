"""Compare two fully accepted bounded Kepler native states at one time."""

import argparse
import json
from pathlib import Path

from prepare_kepler_sparse_area_restart_20261027 import digest


METRICS = (
    "disk_density_abs_mass_weighted",
    "disk_ballistic_density_abs_mass_weighted",
    "disk_native_minus_ballistic_abs_mass_weighted",
    "disk_entropy_log_abs_mass_weighted",
    "disk_radial_speed_rms_over_circular",
    "disk_angular_momentum_relative",
    "total_energy_relative",
)


def load_case(path, tails_path, expected_failed=()):
    receipt = json.loads(path.read_text())
    failed = {key for key, passed in receipt["checks"].items() if not passed}
    if failed != set(expected_failed):
        raise ValueError(f"Native failed checks {sorted(failed)} differ from expected {sorted(set(expected_failed))}: {path}")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(receipt[key])) != receipt[sha_key]:
            raise ValueError(f"Changed {key}: {path}")
    tails = json.loads(tails_path.read_text())
    if (tails["endpoint_sha256"] != receipt["endpoint_sha256"] or
            tails["initial_sha256"] != receipt["initial_sha256"] or
            abs(tails["time"] - receipt["time"]) > 1e-10 or
            set(tails.get("native_failed_checks", [])) != failed):
        raise ValueError(f"Radial-tail receipt does not match native state: {tails_path}")
    return receipt, tails


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("aligned_receipt", type=Path)
    parser.add_argument("jitter_receipt", type=Path)
    parser.add_argument("aligned_tails", type=Path)
    parser.add_argument("jitter_tails", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jitter-expected-failed-check", action="append", default=[])
    args = parser.parse_args()
    aligned, aligned_tails = load_case(args.aligned_receipt, args.aligned_tails)
    jitter, jitter_tails = load_case(args.jitter_receipt, args.jitter_tails,
                                     args.jitter_expected_failed_check)
    if abs(aligned["time"] - jitter["time"]) > 1e-10:
        raise ValueError("Different native physical times")
    metric_rows = {}
    for key in METRICS:
        a, j = float(aligned[key]), float(jitter[key])
        metric_rows[key] = {
            "aligned": a, "jitter": j,
            "jitter_over_aligned_magnitude": abs(j)/abs(a) if a else None,
        }
    inner = {}
    for label, tails in (("aligned", aligned_tails), ("jitter", jitter_tails)):
        row = next((r for r in tails["radial_bands"]
                    if r["r_min"] == .5 and r["r_max"] == .75), None)
        if row is None:
            raise ValueError(f"Missing 0.5<=r<0.75 band: {label}")
        inner[label] = row
    result = {
        "time": aligned["time"],
        "aligned_failed_checks": [],
        "jitter_failed_checks": sorted(args.jitter_expected_failed_check),
        "aligned_receipt": str(args.aligned_receipt.resolve()),
        "jitter_receipt": str(args.jitter_receipt.resolve()),
        "aligned_tails": str(args.aligned_tails.resolve()),
        "jitter_tails": str(args.jitter_tails.resolve()),
        "metrics": metric_rows,
        "inner_disk_r05_r075": inner,
        "inner_density_residual_jitter_over_aligned": (
            inner["jitter"]["density_native_minus_ballistic_abs_mass_weighted"]
            / inner["aligned"]["density_native_minus_ballistic_abs_mass_weighted"]),
        "scope": ("Equal-time native closed-wall W=0 comparison. Explicitly listed failed jitter gates "
                  "remain failures; this exploratory comparison does not turn them into accepted physics. "
                  "One jitter seed and one resolution cannot establish general solver accuracy."),
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
