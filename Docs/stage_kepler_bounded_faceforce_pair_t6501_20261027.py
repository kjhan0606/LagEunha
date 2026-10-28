"""Stage a one-step, selected-face W0 force probe from two accepted states."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
RUNNER = "kepler_q2_bounded_faceforce_pair_t6501_a10_20261027.slurm"
BINARY_SHA = "3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
CASES = {
    "aligned": ("to65_near32_a10_v1", "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json"),
    "jitter": ("jitter005_to65_a10_v1", "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json"),
}
HELPERS = (
    "audit_kepler_bounded_short_replay_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        raise FileExistsError(dest)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)

    sources = {}
    for label, (subdir, receipt_name) in CASES.items():
        receipt = json.loads((DOCS / receipt_name).read_text())
        if (len(receipt["checks"]) != 17 or not all(receipt["checks"].values())
                or abs(receipt["time"] - 6.5) > 1e-10):
            raise ValueError(f"{label}: native parent gate incomplete")
        if digest(Path(receipt["endpoint"])) != receipt["endpoint_sha256"]:
            raise ValueError(f"{label}: accepted checkpoint changed")
        if abs(snapshot_time(Path(receipt["endpoint"]), 55351)-6.5) > 1e-10:
            raise ValueError(f"{label}: accepted checkpoint time changed")
        source = ROOT / subdir
        if Path(receipt["endpoint"]).parent != source:
            raise ValueError(f"{label}: unexpected native source directory")
        if digest(source / "eunha2.exe") != BINARY_SHA:
            raise ValueError(f"{label}: binary differs from accepted near32")
        if digest(source / "params.dat") != PARAMS_SHA:
            raise ValueError(f"{label}: physical parameter file differs")
        sources[label] = (source, receipt)

    initial = Path(sources["jitter"][1]["initial"])
    endpoint = Path(sources["jitter"][1]["endpoint"])
    p0, _ = ordered(initial, 55351)
    p, _ = ordered(endpoint, 55351)
    radius = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    inner = np.flatnonzero((radius >= .5) & (radius < .75))
    rel_pressure = np.abs(p["pressure"].astype(float)/p0["pressure"].astype(float)-1.)
    highest = inner[np.argsort(rel_pressure[inner], kind="stable")[-32:]]
    remaining = np.setdiff1d(inner, highest)
    sampled = np.random.default_rng(20261027).choice(remaining, size=32, replace=False)
    selected = np.sort(np.concatenate((highest, sampled)))
    if len(selected) != 64 or len(np.unique(selected)) != 64:
        raise ValueError("Selected-face ID list is not 64 unique cells")

    dest.mkdir(parents=True)
    shutil.copy2(DOCS / RUNNER, dest / RUNNER)
    with (dest / "selected_ids.txt").open("x") as stream:
        stream.write(",".join(str(int(i)) for i in selected) + "\n")
    with (dest / "selected_ids.sha256").open("x") as stream:
        stream.write(f"{digest(dest / 'selected_ids.txt')}  selected_ids.txt\n")
    staged = {}
    for label, (source, receipt) in sources.items():
        directory = dest / label
        directory.mkdir()
        inputs = {
            "eunha2.exe": source / "eunha2.exe",
            "params.dat": source / "params.dat",
            "khout.000001.dat": Path(receipt["endpoint"]),
            **{name: DOCS / name for name in HELPERS},
        }
        staged[label] = {}
        for name, path in inputs.items():
            shutil.copy2(path, directory / name)
            staged[label][name] = digest(directory / name)
        with (directory / "inputs.sha256").open("x") as stream:
            for name, sha in sorted(staged[label].items()):
                stream.write(f"{sha}  {name}\n")
    manifest = {
        "destination": str(dest), "restart_time": 6.5, "target_time": 6.501,
        "selected_ids": [int(i) for i in selected],
        "highest_pressure_change_ids": [int(i) for i in highest],
        "random_inner_band_ids": [int(i) for i in sampled],
        "selected_ids_sha256": digest(dest / "selected_ids.txt"),
        "runner_sha256": digest(dest / RUNNER),
        "binary_sha256": BINARY_SHA, "params_sha256": PARAMS_SHA,
        "accepted_receipts": {label: sources[label][1]["endpoint_sha256"] for label in CASES},
        "inputs": staged,
        "scope": ("One-step probe using original accepted near32 binary and identical physics. "
                  "64 selected inner-disk cells; face dump is read-only diagnostics."),
    }
    with (dest / "stage_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"destination": str(dest), "ids": len(selected),
                      "runner_sha256": manifest["runner_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
