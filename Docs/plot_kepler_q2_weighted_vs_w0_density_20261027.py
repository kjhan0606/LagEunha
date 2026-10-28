"""Render matched W0 and weighted native density on actual power-cell pixels."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.path import Path as MplPath
import numpy as np
from scipy.spatial import cKDTree

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


N = 55351
DOMAIN = (1., 5.)


def accepted(receipt_path, expected_failed):
    receipt = json.loads(receipt_path.read_text())
    failed = {key for key, passed in receipt["checks"].items() if not passed}
    if failed != set(expected_failed):
        raise ValueError(f"{receipt_path}: failed checks {sorted(failed)}")
    start_key = "initial" if "initial" in receipt else "restart"
    start = Path(receipt[start_key])
    end = Path(receipt["endpoint"])
    if (digest(start) != receipt[start_key + "_sha256"]
            or digest(end) != receipt["endpoint_sha256"]):
        raise ValueError("Changed native input or endpoint")
    time = snapshot_time(end, N)
    if abs(time-float(receipt.get("time", receipt.get("endpoint_time", time))))>1e-10:
        raise ValueError("Wrong native endpoint time")
    return receipt, start, end, time, failed


def certified_power_winners(tree, points, xy, weight):
    """Resolve exact raster power winners using geometric kNN lower bounds."""
    answer = np.empty(len(points), dtype=np.int32)
    pending = np.arange(len(points))
    maximum = float(np.max(weight))
    for k in (8, 16, 32, 64, 128):
        if not len(pending):
            break
        distance, nearest = tree.query(points[pending], k=k, workers=1)
        score = distance*distance-weight[nearest]
        local = np.argmin(score, axis=1)
        row = np.arange(len(pending))
        best = score[row, local]
        # Every unqueried site's Euclidean distance is at least d_k.
        certified = distance[:, -1]**2-maximum > best+1e-12
        answer[pending[certified]] = nearest[row[certified], local[certified]]
        pending = pending[~certified]
    if len(pending):
        raise ValueError(f"Power raster did not certify {len(pending)} pixels with k<=128")
    return answer


def case_fields(start, endpoint, points, weighted):
    p0, _ = ordered(start, N)
    p1, _ = ordered(endpoint, N)
    if not np.array_equal(p0["mass"], p1["mass"]):
        raise ValueError("Changed cell masses")
    if np.any(p0["w2"] != 0):
        raise ValueError("Reference requires W0")
    if not weighted and np.any(p1["w2"] != 0):
        raise ValueError("W0 endpoint has weights")
    xy = np.column_stack((p1["x"].astype(float), p1["y"].astype(float)))
    tree = cKDTree(xy)
    if weighted:
        weights = p1["w2"].astype(float)
        if np.any(~np.isfinite(weights)):
            raise ValueError("Nonfinite power weights")
        ids = certified_power_winners(tree, points, xy, weights)
    else:
        ids = tree.query(points, k=1, workers=1)[1]
    rho = p1["den"].astype(float)
    rho0 = p0["den"].astype(float)
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    disk = (r0>=.5)&(r0<2.)
    relative = rho/rho0-1.
    relative[~disk] = np.nan
    mass = p0["mass"][disk].astype(float)
    disk_abs = float(np.average(np.abs(rho[disk]/rho0[disk]-1.), weights=mass))
    return rho[ids], relative[ids], disk_abs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("w0_receipt", type=Path)
    ap.add_argument("weighted_receipt", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--pixels", type=int, default=512)
    ap.add_argument("--relative-max", type=float, default=.005)
    ap.add_argument("--domain", type=float, nargs=2, default=DOMAIN,
                    metavar=("LOW", "HIGH"),
                    help="Square x/y plotting window; default is the full 1..5 field")
    ap.add_argument("--max-time-difference", type=float, default=1e-10,
                    help="Maximum checkpoint time difference, at most 1e-6")
    ap.add_argument("--w0-expected-failed-check", action="append", default=[])
    ap.add_argument("--weighted-expected-failed-check", action="append", default=[])
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    w0, start0, end0, t0, failed0 = accepted(
        args.w0_receipt, args.w0_expected_failed_check)
    weighted, start1, end1, t1, failed1 = accepted(
        args.weighted_receipt, args.weighted_expected_failed_check)
    if not np.isfinite(args.max_time_difference) or not 0 <= args.max_time_difference <= 1e-6:
        raise ValueError("Time tolerance must be finite and between 0 and 1e-6")
    if (not np.all(np.isfinite(args.domain))
            or not 1. <= args.domain[0] < args.domain[1] <= 5.):
        raise ValueError("Plot domain must be a finite ordered subset of 1..5")
    if digest(start0) != digest(start1) or abs(t0-t1)>args.max_time_difference:
        raise ValueError("Unmatched reference snapshots or physical times")
    axis = np.linspace(*args.domain, args.pixels)
    xx, yy = np.meshgrid(axis, axis)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    theta = np.arange(1024)*2*np.pi/1024
    polygon = np.column_stack((3+2.75*np.cos(theta), 3+2.75*np.sin(theta)))
    inside = MplPath(polygon).contains_points(points).reshape(xx.shape)
    fig, panels = plt.subplots(2, 2, figsize=(10.5, 9), constrained_layout=True)
    for column, (label, start, endpoint, mode, failed) in enumerate((
            ("Original W0", start0, end0, False, failed0),
            ("Experimental weighted", start1, end1, True, failed1))):
        rho, relative, disk_abs = case_fields(start, endpoint, points, mode)
        rho = np.ma.masked_where(~inside, rho.reshape(xx.shape))
        relative = np.ma.masked_invalid(np.where(inside, relative.reshape(xx.shape), np.nan))
        density_artist = panels[0, column].imshow(
            rho, origin="lower", extent=(*args.domain, *args.domain),
            norm=LogNorm(vmin=.01, vmax=1.1, clip=True), cmap="viridis",
            interpolation="nearest")
        residual_artist = panels[1, column].imshow(
            relative, origin="lower", extent=(*args.domain, *args.domain),
            norm=SymLogNorm(linthresh=args.relative_max/500,
                            vmin=-args.relative_max, vmax=args.relative_max,
                            base=10), cmap="coolwarm", interpolation="nearest")
        status = "pass" if not failed else "failed: " + ", ".join(sorted(failed))
        panels[0, column].set_title(f"{label} ({status})")
        panels[1, column].set_title(f"disk mean |ΔΣ/Σ₀| = {disk_abs:.4g}")
    fig.colorbar(density_artist, ax=panels[0, :].tolist(), label="surface density Σ")
    fig.colorbar(residual_artist, ax=panels[1, :].tolist(), label="Σ/Σ₀ − 1")
    for panel in panels.flat:
        panel.set_aspect("equal")
        panel.set_xlabel("x")
        panel.set_ylabel("y")
        ticks = np.linspace(*args.domain, 5)
        panel.set_xticks(ticks)
        panel.set_yticks(ticks)
    fig.suptitle(f"Native bounded Kepler density at t≈{t0:.6g}; Δt={abs(t0-t1):.2g}; "
                 f"window {args.domain[0]:g}..{args.domain[1]:g}; common colors")
    fig.savefig(args.output, dpi=150)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
