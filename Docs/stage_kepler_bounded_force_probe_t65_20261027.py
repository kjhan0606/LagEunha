"""Stage matched one-step native force-input captures from accepted t=6.5 states."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from audit_kepler_sparse_native_ballistic_20261027 import disk_geometry
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DEST = ROOT / "force_inputs_t65_a10_v1"
PARENTS = {"aligned": ROOT / "to65_near32_a10_v1",
           "jitter": ROOT / "jitter005_to65_a10_v1"}
RUNNER = "kepler_q2_bounded_force_probe_t65_a10_20261027.slurm"
AUDITOR = "audit_kepler_bounded_force_probe_t65_20261027.py"
HELPERS = ("audit_kepler_bounded_short_replay_20261027.py",
           "audit_kepler_sparse_area_t0825_20261027.py",
           "audit_kepler_sparse_state_consistency_20261027.py",
           "prepare_kepler_sparse_area_restart_20261027.py")
BINARY_SHA = "3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
N = 55351
EDGES = (.5, .75, 1., 1.25, 1.5, 1.75, 2.)


def accepted(name):
    parent = PARENTS[name]
    receipt = json.loads((parent / "audit_t65.json").read_text())
    if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
            or not all(receipt["checks"].values())):
        raise ValueError(f"Unaccepted {name} t=6.5 state")
    if (digest(Path(receipt["initial"])) != receipt["initial_sha256"]
            or digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]
            or abs(snapshot_time(Path(receipt["endpoint"]), N) - 6.5) > 1e-10):
        raise ValueError(f"Changed {name} accepted snapshot")
    if digest(parent / "eunha2.exe") != BINARY_SHA:
        raise ValueError(f"Changed {name} production binary")
    return receipt


def select_ids(jitter):
    initial, _ = ordered(Path(jitter["initial"]), N)
    endpoint, _ = ordered(Path(jitter["endpoint"]), N)
    r0 = np.hypot(initial["x"].astype(float)-3., initial["y"].astype(float)-3.)
    disk = np.flatnonzero((r0 >= .5) & (r0 < 2.))
    _, ballistic_area, safety = disk_geometry(initial, 6.5, disk)
    if safety < 1:
        raise ValueError("Exact-orbit area depended on omitted images")
    residual = np.full(N, np.nan)
    residual[disk] = ballistic_area / endpoint["volume"][disk].astype(float) - 1.
    rng = np.random.default_rng(20261027)
    chosen = {}
    for low, high in zip(EDGES[:-1], EDGES[1:]):
        eligible = np.flatnonzero((r0 >= low) & (r0 < high))
        for idx in rng.choice(eligible, size=8, replace=False):
            chosen[int(idx)] = "stratified_random"
    for idx in disk[np.argsort(np.abs(residual[disk]))[::-1]]:
        if len(chosen) == 64:
            break
        chosen.setdefault(int(idx), "large_jitter_residual")
    if len(chosen) != 64:
        raise ValueError("Did not select exactly 64 distinct disk cells")
    return [{"id": i, "initial_radius": float(r0[i]),
             "jitter_density_residual": float(residual[i]), "selection": chosen[i]}
            for i in sorted(chosen)]


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    receipts = {name: accepted(name) for name in PARENTS}
    selected = select_ids(receipts["jitter"])
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    DEST.mkdir()
    shutil.copy2(DOCS / RUNNER, DEST / RUNNER)
    shutil.copy2(DOCS / AUDITOR, DEST / AUDITOR)
    (DEST / "selected_ids.json").write_text(json.dumps(selected, indent=2) + "\n")
    with (DEST / "runner.sha256").open("x") as stream:
        for target in (RUNNER, AUDITOR, "selected_ids.json"):
            stream.write(f"{digest(DEST / target)}  {target}\n")
    sources = {}
    for name, parent in PARENTS.items():
        out = DEST / name
        out.mkdir()
        files = {"eunha2.exe": parent / "eunha2.exe",
                 "params.dat": parent / "params.dat",
                 "glass.bin": parent / "glass.bin",
                 "areas.bin": parent / "areas.bin",
                 "khout.000001.dat": Path(receipts[name]["endpoint"]),
                 **{helper: DOCS / helper for helper in HELPERS}}
        for target, source in files.items():
            shutil.copy2(source, out / target)
        (out / "selected_ids.txt").write_text(
            ",".join(str(row["id"]) for row in selected) + "\n")
        with (out / "inputs.sha256").open("x") as stream:
            for target in sorted((*files, "selected_ids.txt")):
                stream.write(f"{digest(out / target)}  {target}\n")
        sources[name] = {"accepted_receipt": str(parent / "audit_t65.json"),
                         "accepted_endpoint_sha256": receipts[name]["endpoint_sha256"],
                         "staged_restart_sha256": digest(out / "khout.000001.dat"),
                         "binary_sha256": digest(out / "eunha2.exe")}
    manifest = {"run": str(DEST), "start_time": 6.5, "stop_time": 6.5005,
                "selected_cells": len(selected), "selection_seed": 20261027,
                "selection": "8 uniform random IDs in each of six initial-radius disk bands, plus 16 largest jitter native-minus-ballistic density residuals",
                "runner_sha256": digest(DEST / RUNNER),
                "auditor_sha256": digest(DEST / AUDITOR), "sources": sources,
                "scope": "One accepted-state GPU face-force input capture per IC; observational, no physics change or long-time claim"}
    (DEST / "staging.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
