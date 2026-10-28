"""Verify that the native opt-in direct-U RK4 branch used its stored anchor."""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_polar_long_gate_20261027 import ordered


N = 384 * 384
STRIDE = 576
ANCHOR_OFFSET = 472  # C offsetof(treevorostressrk4particletype, stress.E_inv_xx)


def anchors_by_id(path):
    dtype = np.dtype({"names": ("id_word", "anchor"),
                      "formats": ("<u8", "<f8"),
                      "offsets": (8, ANCHOR_OFFSET), "itemsize": STRIDE})
    data = np.memmap(path, mode="r", dtype=dtype, offset=16, shape=(N,))
    ids = (data["id_word"] & np.uint64((1 << 56) - 1)).astype(np.int64)
    order = np.argsort(ids)
    if not np.array_equal(ids[order], np.arange(N)):
        raise ValueError("Incomplete native ID set in direct-U anchor payload")
    return np.asarray(data["anchor"][order], dtype=float)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accepted_start", type=Path)
    parser.add_argument("native_receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = json.loads(args.native_receipt.read_text())
    endpoint = Path(receipt["endpoint"])
    start, ts = ordered(args.accepted_start)
    end, te = ordered(endpoint)
    a0 = anchors_by_id(args.accepted_start)
    a1 = anchors_by_id(endpoint)
    u = end["ie"].astype(float)
    checks = {
        "accepted_start_t0825": bool(abs(ts - .825) < 1e-10),
        "endpoint_t0845": bool(abs(te - .845) < 1e-10),
        "accepted_start_anchor_zero": bool(np.all(a0 == 0)),
        "all_endpoint_anchors_positive_finite": bool(np.all(np.isfinite(a1) & (a1 > 0))),
        "anchor_tracks_final_U": bool(np.max(abs(a1 / u - 1.)) < .1),
        "native_physics_gate": bool(all(receipt["checks"].values())),
        "mass_exact": bool(np.array_equal(start["mass"], end["mass"])),
    }
    report = {
        "accepted_start": str(args.accepted_start), "endpoint": str(endpoint),
        "time_start": ts, "time_end": te, "anchor_offset": ANCHOR_OFFSET,
        "endpoint_anchor_min": float(a1.min()),
        "endpoint_anchor_max": float(a1.max()),
        "max_relative_anchor_vs_final_U": float(np.max(abs(a1 / u - 1.))),
        "checks": checks,
        "scope": "Branch-use certificate for opt-in same-equation RK4 direct-U arithmetic; not a physical or t>6 stability result",
    }
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))
    if not all(checks.values()):
        raise SystemExit("Native direct-U anchor or physical gate failed")


if __name__ == "__main__":
    main()
