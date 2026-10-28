"""Track along-ring pressure roughness that directly enters face torque."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipts", nargs="+", type=Path)
    parser.add_argument("--reference-receipt", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipts = [json.loads(path.read_text()) for path in args.receipts]
    if not receipts:
        raise ValueError("Need native receipts")
    sha0 = receipts[0]["initial_sha256"]
    n = receipts[0]["n"]
    initial_path = Path(receipts[0]["initial"])
    if digest(initial_path) != sha0:
        raise ValueError("Changed initial snapshot")
    initial, _ = ordered(initial_path, n)
    x0 = initial["x"].astype(float)-3.
    y0 = initial["y"].astype(float)-3.
    radius = np.hypot(x0, y0)
    ids = np.flatnonzero((radius >= .5) & (radius < .52))
    if len(ids) != 294:
        raise ValueError("Unexpected q2 first-ring population")
    angle = np.arctan2(y0[ids], x0[ids])
    order = np.argsort(angle)
    ids, angle = ids[order], angle[order]
    gaps = np.diff(np.r_[angle, angle[0]+2*np.pi])
    initial_pressure = initial["pressure"][ids].astype(float)
    if np.max(initial_pressure)-np.min(initial_pressure) != 0.:
        raise ValueError("Initial same-ring pressures differ")
    reference_error = None
    reference = None
    if args.reference_receipt:
        reference = json.loads(args.reference_receipt.read_text())
        if (reference["initial_sha256"] != sha0 or reference["n"] != n
                or digest(Path(reference["endpoint"])) != reference["endpoint_sha256"]):
            raise ValueError("Changed or mismatched reference endpoint")
        state, _ = ordered(Path(reference["endpoint"]), n)
        reference_error = state["pressure"][ids].astype(float)/initial_pressure-1.
    rows = []
    for path, receipt in zip(args.receipts, receipts):
        if (receipt["initial_sha256"] != sha0 or receipt["n"] != n
                or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]):
            raise ValueError(f"Mismatched or changed endpoint: {path}")
        end, _ = ordered(Path(receipt["endpoint"]), n)
        error = end["pressure"][ids].astype(float)/initial_pressure-1.
        jump = np.roll(error, -1)-error
        angle_gradient = jump/gaps
        spectrum = np.abs(np.fft.rfft(error-np.mean(error)))**2
        high_fraction = float(spectrum[30:].sum()/spectrum[1:].sum())
        rows.append({"receipt": str(path), "time": receipt["time"],
                     "endpoint_sha256": receipt["endpoint_sha256"],
                     "gates_passed": sum(receipt["checks"].values()),
                     "gates_total": len(receipt["checks"]),
                     "pressure_relative_mean_abs": float(np.mean(np.abs(error))),
                     "pressure_relative_rms": float(np.sqrt(np.mean(error**2))),
                     "angular_neighbor_jump_rms": float(np.sqrt(np.mean(jump**2))),
                     "angular_gradient_rms": float(np.sqrt(np.mean(angle_gradient**2))),
                     "fft_mode_ge30_power_fraction_approx": high_fraction})
        if reference_error is not None:
            increment = error-reference_error
            neighbor_increment = np.roll(increment, -1)-increment
            rows[-1]["pressure_error_increment_rms_from_reference"] = float(
                np.sqrt(np.mean(increment**2)))
            rows[-1]["angular_neighbor_increment_rms_from_reference"] = float(
                np.sqrt(np.mean(neighbor_increment**2)))
    result = {"initial_sha256": sha0, "first_ring_sites": len(ids),
              "initial_pressure_jump_max": 0.,
              "initial_angle_gap_max_over_min": float(gaps.max()/gaps.min()),
              "reference_receipt": str(args.reference_receipt)
                  if args.reference_receipt else None,
              "reference_endpoint_sha256": reference["endpoint_sha256"]
                  if reference else None,
              "rows": rows,
              "scope": "Same initial q2 first-ring angular order for all endpoints. "
                       "Neighbor jumps are a pressure-roughness proxy, not "
                       "native face traction; the FFT band is approximate "
                       "because angular spacing is not perfectly uniform."}
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    for row in rows:
        print(row["time"], Path(row["receipt"]).parent.name,
              row["pressure_relative_mean_abs"], row["angular_neighbor_jump_rms"],
              row["fft_mode_ge30_power_fraction_approx"])


if __name__ == "__main__":
    main()
