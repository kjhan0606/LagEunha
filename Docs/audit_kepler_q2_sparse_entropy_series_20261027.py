"""Compare saved full-q2 and corrected-sparse entropy histories without rerunning gas."""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


GM1 = float(np.float32(1.666667)) - 1.
BINS = ((.5, .75), (.75, 1.), (1., 2.), (.5, 2.))


def weighted_mean(value, weights):
    return float(np.sum(value * weights, dtype=np.longdouble) /
                 np.sum(weights, dtype=np.longdouble))


def step(path):
    match = re.fullmatch(r"khout\.(\d+)\.dat", path.name)
    if not match:
        raise ValueError(f"Not a native numbered checkpoint: {path}")
    return int(match.group(1))


def analyze(initial, checkpoints, count):
    p0, m0 = ordered(initial, count)
    if abs(snapshot_time(initial, count)) > 1e-12 or abs(float(m0["t"])) > 1e-12:
        raise ValueError("Initial native state is not t=0")
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    mass = p0["mass"].astype(float)
    if not np.all(np.isfinite(mass)) or not np.all(mass > 0):
        raise ValueError("Invalid initial mass")
    masks = [(r0 >= lo) & (r0 < hi) for lo, hi in BINS]
    rows = []
    last_time = -1.
    for path in checkpoints:
        t = snapshot_time(path, count)
        p, meta = ordered(path, count)
        if (t <= last_time or abs(float(meta["t"])-t) > 1e-12 or
                not np.array_equal(p0["mass"], p["mass"]) or
                not np.all(p0["w2"] == 0) or not np.all(p["w2"] == 0)):
            raise ValueError(f"Time, fixed-mass, or W0 gate failed: {path}")
        if (not all(np.isfinite(p[key]).all() for key in
                    ("volume", "ie", "den", "pressure")) or
                not np.all(p["volume"] > 0) or not np.all(p["ie"] > 0)):
            raise ValueError(f"Invalid native thermodynamic state: {path}")
        rho_consistency = float(np.max(np.abs(
            p["den"].astype(float)/(mass/p["volume"].astype(float))-1.)))
        eos_consistency = float(np.max(np.abs(
            p["pressure"].astype(float)/(GM1*p["ie"].astype(float)/p["volume"].astype(float))-1.)))
        if (rho_consistency >= 1e-6 or eos_consistency >= 1e-6 or
                abs(float(np.sum(p["volume"].astype(float), dtype=np.longdouble))-36.) >= 1e-5):
            raise ValueError(f"Native density/EOS/periodic-area gate failed: {path}")
        density = p["den"].astype(float)/p0["den"].astype(float)-1.
        entropy = (np.log(p["ie"].astype(float)/p0["ie"].astype(float)) +
                   GM1*np.log(p["volume"].astype(float)/p0["volume"].astype(float)))
        bins = []
        for (lo, hi), take in zip(BINS, masks):
            w = mass[take]
            bins.append({"rlo":lo,"rhi":hi,"cells":int(np.sum(take)),
                         "density_abs":weighted_mean(np.abs(density[take]), w),
                         "entropy_abs":weighted_mean(np.abs(entropy[take]), w),
                         "entropy_signed":weighted_mean(entropy[take], w)})
        rows.append({"checkpoint":str(path.resolve()),"sha256":digest(path),
                     "step":step(path),"time":t,"bins":bins,
                     "rho_mass_volume_max_relative":rho_consistency,
                     "eos_pressure_max_relative":eos_consistency})
        last_time = t
    return {"initial":str(initial.resolve()),"initial_sha256":digest(initial),
            "count":count,"rows":rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full-initial", type=Path, required=True)
    ap.add_argument("--full", type=Path, action="append", required=True)
    ap.add_argument("--sparse-initial", type=Path, required=True)
    ap.add_argument("--sparse", type=Path, action="append", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    if len(args.full) != len(args.sparse) or len(args.full) < 2:
        raise ValueError("Need corresponding full and sparse checkpoint series")
    full = analyze(args.full_initial, args.full, 384*384)
    sparse = analyze(args.sparse_initial, args.sparse, 54312)
    for a, b in zip(full["rows"], sparse["rows"]):
        if abs(a["time"]-b["time"]) > .002:
            raise ValueError("Checkpoint times differ too much for early comparison")
    result = {"full_q2":full,"sparse_narrow":sparse,
              "scope":"Early native stored-state trend only; times differ by <.002 and changing ICs/timesteps cannot be causally separated"}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(f"Compared {len(args.full)} complete native checkpoint pairs")


if __name__ == "__main__":
    main()
