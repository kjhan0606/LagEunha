"""Start a corrected sparse Kepler continuation only from a passing t=.001 gate."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import struct


BIN_SHA = "eff735174de1b487759f168d75e1e317a201692ef8d6cadba03cd03e49afb4b5"
IC = {
    54312: ("03448f6151f85195b49e8ea7bd7543611094b8b0a957c3891134d1d83f0156a7",
            "828c3431f27e0e0ae8d1652ebc4bba1c9687e96e3c922c45f20a18e2b69b9417"),
    55931: ("de50b97639a243f9384fc7b4bfc402f9eccb52386f30ecdd3450f7b0acaa9e8f",
            "899ef274a318ce2677b716b9a2050c87802d2d110ac4aaf2d6c78d14fe7ff5b5"),
}
REQUIRED_GATES = (
    "periodic_area_sidecar_tiles_box", "initial_native_areas_match_sidecar",
    "initial_native_areas_tile_box", "evolved_native_areas_tile_box",
    "initial_mass_matches_target_profile_and_sidecar",
    "initial_total_mass_matches_profile", "mass_exact_across_step",
    "W0", "positive_state", "no_floor_injection", "disk_orbit",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def snapshot_time(path, count):
    with path.open("rb") as stream:
        header = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        time = struct.unpack("<d", stream.read(8))[0]
    if header != (-1, count, 576, 1) or path.stat().st_size != 16 + 576 * count + 52:
        raise ValueError(f"Incomplete or incompatible native snapshot: {path}")
    return time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("preflight", type=Path)
    ap.add_argument("destination", type=Path)
    ap.add_argument("--expected-count", type=int, choices=tuple(IC), required=True)
    ap.add_argument("--target-time", type=float, required=True)
    args = ap.parse_args()
    if args.target_time <= .001:
        raise ValueError("Target must follow the preflight")
    parent = args.preflight.resolve()
    dest = args.destination.resolve()
    receipt = json.loads((parent / "audit_t001.json").read_text())
    if (receipt["n"] != args.expected_count or abs(receipt["t0"]) > 1e-12 or
            abs(receipt["t1"] - .001) > 1e-12 or
            not all(receipt["checks"].values()) or
            not all(receipt["checks"].get(key) is True for key in REQUIRED_GATES) or
            abs(receipt["initial_volume_sum"] - 36) > 1e-6 or
            abs(receipt["evolved_volume_sum"] - 36) > 1e-6):
        raise ValueError("Corrected periodic-mass preflight did not pass")
    expected_glass, expected_area = IC[args.expected_count]
    for name, sha in (("glass.bin", expected_glass),
                      ("periodic_areas.bin", expected_area),
                      ("eunha2.exe", BIN_SHA)):
        if digest(parent / name) != sha or receipt["sha256"].get(name) != sha:
            raise ValueError(f"Preflight input is not the approved {name}")
    initial = Path(receipt["initial"])
    source = Path(receipt["evolved"])
    if initial.parent.resolve() != parent or source.parent.resolve() != parent:
        raise ValueError("Preflight receipt points outside its run")
    if abs(snapshot_time(initial, args.expected_count)) > 1e-12 or abs(
            snapshot_time(source, args.expected_count) - .001) > 1e-12:
        raise ValueError("Preflight native snapshot time changed")
    match = re.fullmatch(r"khout\.(\d+)\.dat", source.name)
    if (not match or int(match.group(1)) <= 0 or
            digest(initial) != receipt["sha256"].get(initial.name) or
            digest(source) != receipt["sha256"].get(source.name)):
        raise ValueError("Preflight checkpoint identity changed")
    lines = [line for line in (parent / "log").read_text(errors="replace").splitlines()
             if "[RK4E]" in line]
    if not lines:
        raise ValueError("Preflight energy ledger is missing")
    for key in ("floor_cum", "de_cum", "sfl_cum"):
        ledger_match = re.search(rf"\b{key}=([0-9.eE+-]+)", lines[-1])
        if not ledger_match or float(ledger_match.group(1)) != 0.:
            raise ValueError(f"Preflight {key} ledger did not pass")
    if dest.exists() and (any(dest.glob("khout.*.dat")) or
                          any((dest / name).exists() for name in
                              ("eunha2.exe", "params.dat", "restart_plan.json"))):
        raise ValueError("Destination already has native state inputs")
    dest.mkdir(parents=True, exist_ok=True)
    inputs = {}
    for original in (source, parent / "eunha2.exe", parent / "params.dat"):
        target = dest / original.name
        shutil.copy2(original, target)
        inputs[target.name] = digest(target)
    plan = {"preflight_receipt": str(parent / "audit_t001.json"),
            "initial": str(initial), "preflight_checkpoint": str(source),
            "step": int(match.group(1)), "expected_time": .001,
            "target_time": args.target_time, "n": args.expected_count,
            "glass_sha256": expected_glass, "periodic_area_sha256": expected_area,
            "binary_sha256": BIN_SHA, "inputs_sha256": inputs}
    with (dest / "restart_plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2)
        stream.write("\n")
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
