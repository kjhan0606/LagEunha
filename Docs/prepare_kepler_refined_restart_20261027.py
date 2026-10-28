"""Stage a refined W=0 continuation only from its passing t=.001 gate."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import struct


COUNT = 420 * 420
STRIDE = 576
BIN_SHA = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"
REQUIRED_GATES = (
    "periodic_area_sidecar_tiles_box", "initial_native_areas_match_sidecar",
    "initial_native_areas_tile_box", "evolved_native_areas_tile_box",
    "initial_mass_matches_target_profile_and_sidecar",
    "initial_total_mass_matches_profile", "mass_exact_across_step", "W0",
    "positive_state", "no_floor_injection", "glass_initial_positions",
    "disk_pressure", "disk_orbit",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def snapshot_time(path):
    with path.open("rb") as stream:
        header = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        time = struct.unpack("<d", stream.read(8))[0]
    if (header != (-1, COUNT, STRIDE, 1) or
            path.stat().st_size != 16 + COUNT * STRIDE + 52):
        raise ValueError(f"Incomplete or incompatible checkpoint: {path}")
    return time


def require_zero_injection(log_path):
    lines = [line for line in log_path.read_text(errors="replace").splitlines()
             if "[RK4E]" in line]
    if not lines:
        raise ValueError("Preflight RK4 energy ledger is missing")
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        match = re.search(rf"\b{key}=([0-9.eE+-]+)", lines[-1])
        if not match or float(match.group(1)) != 0.:
            raise ValueError(f"Preflight {key} is missing or nonzero")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("preflight", type=Path)
    ap.add_argument("destination", type=Path)
    ap.add_argument("--target-time", type=float, required=True)
    args = ap.parse_args()
    parent = args.preflight.resolve()
    dest = args.destination.resolve()
    if dest == parent or args.target_time <= .001:
        raise ValueError("Expected an independent later-time destination")

    receipt_path = parent / "audit_t001.json"
    receipt = json.loads(receipt_path.read_text())
    if (receipt.get("n") != COUNT or abs(receipt.get("t0", float("inf"))) > 1e-12 or
            abs(receipt.get("t1", -1.) - .001) > 1e-12 or
            not all(receipt["checks"].values()) or
            not all(receipt["checks"].get(key) is True for key in REQUIRED_GATES)):
        raise ValueError("Refined W=0 startup gate did not pass")
    if (abs(receipt["initial_volume_sum"] - 36.) > 1e-6 or
            abs(receipt["evolved_volume_sum"] - 36.) > 1e-6):
        raise ValueError("Preflight native areas do not tile the box")

    inputs = json.loads((parent / "prepare_manifest.json").read_text())
    areas = json.loads((parent / "periodic_areas.json").read_text())
    expected = {"glass.bin": inputs["glass_sha256"],
                "periodic_areas.bin": areas["area_sha256"],
                "eunha2.exe": BIN_SHA}
    for name, sha in expected.items():
        if digest(parent / name) != sha or receipt["sha256"].get(name) != sha:
            raise ValueError(f"Preflight input changed: {name}")
    if digest(parent / "params.dat") != inputs["staged_params_sha256"]:
        raise ValueError("Preflight parameters changed")

    initial = Path(receipt["initial"]).resolve()
    checkpoint = Path(receipt["evolved"]).resolve()
    if (initial.parent != parent or checkpoint.parent != parent or
            abs(snapshot_time(initial)) > 1e-12 or
            abs(snapshot_time(checkpoint) - .001) > 1e-12):
        raise ValueError("Preflight checkpoint lineage or time is invalid")
    match = re.fullmatch(r"khout\.(\d+)\.dat", checkpoint.name)
    if (not match or int(match.group(1)) <= 0 or
            digest(initial) != receipt["sha256"].get(initial.name) or
            digest(checkpoint) != receipt["sha256"].get(checkpoint.name)):
        raise ValueError("Preflight checkpoint identity changed")
    require_zero_injection(parent / "log")

    if dest.exists() and any((dest / name).exists() for name in
                              (checkpoint.name, "eunha2.exe", "params.dat",
                               "restart_plan.json")):
        raise ValueError("Destination already contains restart state")
    dest.mkdir(parents=True, exist_ok=True)
    copied = {}
    for source in (checkpoint, parent / "eunha2.exe", parent / "params.dat"):
        target = dest / source.name
        shutil.copy2(source, target)
        copied[target.name] = digest(target)
    plan = {
        "preflight_receipt": str(receipt_path),
        "initial": str(initial), "initial_sha256": digest(initial),
        "preflight_checkpoint": str(checkpoint),
        "preflight_checkpoint_sha256": digest(checkpoint),
        "restart_step": int(match.group(1)), "particle_count": COUNT,
        "predecessor_time": .001, "target_time": args.target_time,
        "glass_sha256": expected["glass.bin"],
        "periodic_area_sha256": expected["periodic_areas.bin"],
        "binary_sha256": BIN_SHA, "inputs_sha256": copied,
    }
    with (dest / "restart_plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2)
        stream.write("\n")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
