"""Compare same-start native q1 RK4 energy arithmetic with direct-U opt-in."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_polar_long_gate_20261027 import ordered, total_energy
from audit_kepler_q1_entropy_k_pairwise_20261027 import REGIONS, region_audit, weighted_mean


def run_wall_seconds(run):
    logs = list(run.glob("slurm-*.out"))
    if len(logs) != 1:
        return None
    found = re.findall(r"\[RUN_TIME\] wall=([0-9.]+)",
                       logs[0].read_text(errors="replace"))
    return float(found[-1]) if len(found) == 1 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline_bundle", type=Path)
    ap.add_argument("direct_bundle", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    base_manifest = json.loads((args.baseline_bundle / "manifest.json").read_text())
    direct_manifest = json.loads((args.direct_bundle / "manifest.json").read_text())
    if (base_manifest["parent_checkpoint_sha256"] != direct_manifest["start_sha256"]
            or base_manifest["old_binary_sha256"] != direct_manifest["binary_sha256"]
            or base_manifest["runs"]["energy"]["entropy_mode"] != 0):
        raise ValueError("Not the same accepted q1 binary/start or baseline mode")
    start_path = Path(base_manifest["parent_checkpoint"])
    if start_path.resolve() != Path(direct_manifest["accepted_start"]).resolve():
        raise ValueError("Reference and direct-U runs do not share the exact checkpoint")
    start, ts = ordered(start_path)
    reference, t0 = ordered(Path(base_manifest["initial"]))
    if abs(ts - .825) > 1e-10 or t0 != 0:
        raise ValueError("Unexpected start/reference times")
    mass = start["mass"].astype(float)
    radius = np.hypot(reference["x"].astype(float) - 3.,
                      reference["y"].astype(float) - 3.)
    disk = (radius >= .5) & (radius < 2.)
    cases = {}
    endpoints = {}
    runs = {"energy": args.baseline_bundle / "energy", "direct_u": args.direct_bundle}
    for name, run in runs.items():
        receipt = json.loads((run / "audit_t0845.json").read_text())
        endpoint = Path(receipt["endpoint"])
        if endpoint.parent.resolve() != run.resolve():
            raise ValueError(f"Endpoint outside expected run: {name}")
        end, te = ordered(endpoint)
        if (abs(te - .845) > 1e-10 or
                not np.array_equal(start["mass"], end["mass"]) or
                np.any(end["w2"] != 0) or np.any(end["ie"] <= 0) or
                np.any(end["volume"] <= 0)):
            raise ValueError(f"Native state incompatibility: {name}")
        if "EXIT:0" not in (run / "log").read_text(errors="replace"):
            raise ValueError(f"Native executable failed: {name}")
        radial_path = run / "audit_radial_t0845.json"
        radial = json.loads(radial_path.read_text()) if radial_path.exists() else None
        cases[name] = {
            "endpoint": str(endpoint), "time": te,
            "native_checks": receipt["checks"],
            "native_gate_pass": bool(all(receipt["checks"].values())),
            "radial_checks": radial["checks"] if radial else None,
            "wall_seconds": run_wall_seconds(run),
            "total_energy_relative_from_start":
                (total_energy(end) - total_energy(start)) / abs(total_energy(start)),
            "regions": {label: region_audit(start, end, radius, mass, lo, hi)
                        for label, lo, hi in REGIONS},
        }
        endpoints[name] = end
    anchor = json.loads((args.direct_bundle / "audit_anchor_t0845.json").read_text())
    if not all(anchor["checks"].values()):
        raise ValueError("Direct-U branch-use certificate failed")
    base, direct = endpoints["energy"], endpoints["direct_u"]
    report = {
        "accepted_start": str(start_path), "time_start": ts, "time_end": .845,
        "cases": cases, "direct_u_anchor_receipt": anchor,
        "same_time_disk_difference": {
            "position_abs_mass_weighted": weighted_mean(
                np.hypot(direct["x"][disk] - base["x"][disk],
                         direct["y"][disk] - base["y"][disk]), mass[disk]),
            "density_relative_abs_mass_weighted": weighted_mean(
                abs(direct["den"][disk] / base["den"][disk] - 1.), mass[disk]),
        },
        "scope": "Same-checkpoint same-equation native RK4 arithmetic attribution only; no t>6 stability or general speedup claim",
    }
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"energy_disk": cases["energy"]["regions"]["disk"],
                      "direct_disk": cases["direct_u"]["regions"]["disk"],
                      "same_time": report["same_time_disk_difference"]}, indent=2))


if __name__ == "__main__":
    main()
