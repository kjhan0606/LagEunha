"""Compare equal-time accepted aligned and site-jitter bounded Kepler runs."""

import argparse
import json
from pathlib import Path

from prepare_kepler_sparse_area_restart_20261027 import digest


def load(native_path, tails_path):
    native = json.loads(native_path.read_text())
    tails = json.loads(tails_path.read_text())
    if not all(native["checks"].values()):
        raise ValueError(f"Failed native physics gate: {native_path}")
    if tails["endpoint_sha256"] != native["endpoint_sha256"] or tails["initial_sha256"] != native["initial_sha256"]:
        raise ValueError(f"Radial-tail receipt uses different native state: {tails_path}")
    if abs(float(tails["time"]) - float(native["time"])) > 1e-10:
        raise ValueError(f"Radial-tail receipt uses different physical time: {tails_path}")
    if digest(Path(native["initial"])) != native["initial_sha256"] or digest(Path(native["endpoint"])) != native["endpoint_sha256"]:
        raise ValueError(f"Native state changed after audit: {native_path}")
    inner = next((row for row in tails["radial_bands"]
                  if row["r_min"] == .5 and row["r_max"] == .75), None)
    if inner is None:
        raise ValueError("Missing inner disk band")
    return native, tails, inner


def ratio(a, b):
    if b == 0:
        raise ValueError("Undefined ratio")
    return float(a / b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("aligned_native", type=Path)
    parser.add_argument("aligned_tails", type=Path)
    parser.add_argument("jitter_native", type=Path)
    parser.add_argument("jitter_tails", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    a, at, ai = load(args.aligned_native, args.aligned_tails)
    j, jt, ji = load(args.jitter_native, args.jitter_tails)
    if abs(float(a["time"]) - float(j["time"])) > 1e-10:
        raise ValueError("The native endpoints are not at the same physical time")
    if a["n"] != j["n"] or a["disk_sites"] != j["disk_sites"]:
        raise ValueError("The two discretizations have different cell counts")
    keys = (
        "disk_density_abs_mass_weighted",
        "disk_ballistic_density_abs_mass_weighted",
        "disk_native_minus_ballistic_abs_mass_weighted",
        "disk_entropy_log_abs_mass_weighted",
        "disk_position_error_abs_mass_weighted",
        "disk_radial_speed_rms_over_circular",
        "disk_angular_momentum_relative",
        "total_energy_relative",
    )
    comparisons = {key: {"aligned": a[key], "jitter": j[key],
                         "jitter_over_aligned": ratio(j[key], a[key])
                         if a[key] != 0 else None} for key in keys}
    comparisons["disk_entropy_log_signed_mass_weighted"] = {
        "aligned": at["whole_disk"]["entropy_log_signed_mass_weighted"],
        "jitter": jt["whole_disk"]["entropy_log_signed_mass_weighted"],
    }
    inner_keys = (
        "entropy_log_signed_mass_weighted",
        "entropy_log_abs_mass_weighted",
        "density_native_minus_ballistic_abs_mass_weighted",
        "density_native_minus_ballistic_abs_p99",
        "density_native_minus_ballistic_abs_max",
        "velocity_error_over_circular_rms_mass_weighted",
    )
    inner = {key: {"aligned": ai[key], "jitter": ji[key],
                   "jitter_over_aligned": ratio(ji[key], ai[key])
                   if ai[key] != 0 else None} for key in inner_keys}
    result = {
        "time": float(a["time"]),
        "aligned_native_receipt": str(args.aligned_native.resolve()),
        "aligned_tails_receipt": str(args.aligned_tails.resolve()),
        "jitter_native_receipt": str(args.jitter_native.resolve()),
        "jitter_tails_receipt": str(args.jitter_tails.resolve()),
        "aligned_endpoint_sha256": a["endpoint_sha256"],
        "jitter_endpoint_sha256": j["endpoint_sha256"],
        "native_density_change_residual_fraction": {
            "aligned": ratio(a["disk_native_minus_ballistic_abs_mass_weighted"], a["disk_density_abs_mass_weighted"]),
            "jitter": ratio(j["disk_native_minus_ballistic_abs_mass_weighted"], j["disk_density_abs_mass_weighted"]),
        },
        "whole_disk": comparisons,
        "inner_disk_0p5_to_0p75": inner,
        "scope": "Equal-time two-IC comparison of accepted fixed-mass W0 native states. Different site locations and one-time cell masses; ratios do not isolate a force mechanism or prove late-time stability.",
    }
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
