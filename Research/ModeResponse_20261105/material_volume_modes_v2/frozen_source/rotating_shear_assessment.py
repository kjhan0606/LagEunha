"""Analyse native checkpoint fields; no fluid evolution is performed here."""
from pathlib import Path
import argparse
import csv
import json
import re
import struct

import numpy as np


def read_native(path):
    path = Path(path)
    with path.open("rb") as stream:
        header = struct.unpack("<iIQIIQQQIIQ", stream.read(64))
        stream.seek(-52, 2)
        time = struct.unpack("<d", stream.read(8))[0]
    magic, version, count, record, prefix, _, tag, payload, size, _, _ = header
    assert (magic, version, record, prefix, tag, size) == (
        -2147483607, 2, 696, 576, 0x4746535736393632, 64
    )
    assert payload == count * record
    assert path.stat().st_size == size + payload + 52
    raw = np.memmap(path, dtype="u1", mode="r", offset=size, shape=(count, record))

    def field(offset, dtype="f8"):
        return np.ndarray((count,), dtype=dtype, buffer=raw,
                          offset=offset, strides=(record,)).copy()

    ids = field(8, "u8") & np.uint64(2**56 - 1)
    order = np.argsort(ids)
    assert np.unique(ids).size == count
    fields = {key: field(offset)[order] for key, offset in (
        ("x", 24), ("y", 32), ("m", 48), ("vx", 56), ("vy", 64),
        ("W", 576), ("pi", 584), ("U", 672), ("V", 688))}
    fields["id"] = ids[order]
    assert all(np.isfinite(v).all() for v in fields.values())
    return fields, time


def energy(state):
    return float(np.sum(state["U"], dtype=np.longdouble) + np.sum(
        0.5 * state["m"] * (state["vx"]**2 + state["vy"]**2),
        dtype=np.longdouble))


def assess(directory):
    cfg = json.loads((directory / "registration.json").read_text())
    birth, birth_time = read_native(cfg["birth_initial"])
    assert birth_time == 0
    gamma = float(np.float32(5 / 3))
    domain_volume = 1 if cfg["problem"] == "gresho" else 2
    initial_energy = energy(birth)
    entropy0 = np.log((gamma - 1) * birth["U"] / birth["V"] /
                      (birth["m"] / birth["V"])**gamma)
    rows = []
    for path in sorted(directory.glob("khout.*.dat")):
        if path.is_symlink():
            continue
        s, time = read_native(path)
        assert np.array_equal(s["id"], birth["id"])
        rho = s["m"] / s["V"]
        pressure = (gamma - 1) * s["U"] / s["V"]
        entropy = np.log(pressure / rho**gamma)
        row = {
            "time": time, "path": str(path),
            "mass_identical_per_id": bool(np.array_equal(s["m"], birth["m"])),
            "W_and_pi_identically_zero": bool(np.all(s["W"] == 0) and np.all(s["pi"] == 0)),
            "minimum_U": float(s["U"].min()), "minimum_V": float(s["V"].min()),
            "relative_volume_sum_error": float(s["V"].sum() / domain_volume - 1),
            "relative_energy_change_from_birth": energy(s) / initial_energy - 1,
            "momentum_change": [float(np.sum(s["m"] * (s[v] - birth[v]), dtype=np.longdouble))
                                for v in ("vx", "vy")],
            "minimum_material_entropy_change": float((entropy - entropy0).min()),
            "mass_weighted_entropy_change": float(np.sum(s["m"] * (entropy - entropy0), dtype=np.longdouble)),
        }
        if cfg["problem"] == "gresho":
            dx, dy = s["x"] - 0.5, s["y"] - 0.5
            radius = np.hypot(dx, dy)
            speed = np.where(radius <= 0.2, 5 * radius,
                             np.where(radius <= 0.4, 2 - 5 * radius, 0))
            ex = np.divide(dx, radius, out=np.zeros_like(dx), where=radius > 0)
            ey = np.divide(dy, radius, out=np.zeros_like(dy), where=radius > 0)
            azimuthal = -ey * s["vx"] + ex * s["vy"]
            vector_error = np.hypot(s["vx"] + ey * speed, s["vy"] - ex * speed)
            row.update({
                "azimuthal_velocity_L1_unweighted": float(np.mean(abs(azimuthal - speed))),
                "azimuthal_velocity_L1_volume_weighted": float(np.sum(s["V"] * abs(azimuthal - speed))),
                "vector_velocity_L1_volume_weighted": float(np.sum(s["V"] * vector_error)),
                "radial_velocity_L1_volume_weighted": float(np.sum(s["V"] * abs(ex * s["vx"] + ey * s["vy"]))),
                "kinetic_energy_fraction_of_birth": float(np.sum(s["m"] * (s["vx"]**2 + s["vy"]**2)) /
                                                         np.sum(birth["m"] * (birth["vx"]**2 + birth["vy"]**2))),
            })
        rows.append(row)
    result_file = directory / "result.json"
    text = (directory / "log").read_text()
    heat = re.findall(r"\[UNIFIED_RADIAL_HEAT\].*segment_cumulative=\(([^)]*)\)", text)
    radial = [float(x) for x in heat[-1].split(",")] if heat else None
    if radial is not None and cfg.get("retains_preflight"):
        preflight = directory.parent.parent / "smooth_preflight_v1" / cfg["problem"] / "log"
        older = re.findall(r"\[UNIFIED_RADIAL_HEAT\].*segment_cumulative=\(([^)]*)\)", preflight.read_text())
        assert older
        radial = [a + float(b) for a, b in zip(radial, older[-1].split(","))]
    flow = list(csv.DictReader((directory / "unified_flow.csv").open()))
    return {"case": cfg["name"], "gamma": gamma,
            "native_exit_result": json.loads(result_file.read_text()) if result_file.exists() else None,
            "checkpoints": rows, "birth_to_latest_radial_heat": radial,
            "flow_first": flow[0], "flow_latest": flow[-1]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    cases = [assess(d) for d in sorted(args.suite.iterdir())
             if d.is_dir() and (d / "registration.json").is_file()]
    with args.output.open("x") as output:
        json.dump({"scope": "Native field and conservation audit; KH profile changes are not exact-solution error norms.",
                   "cases": cases}, output, indent=2)
        output.write("\n")
    for case in cases:
        last = case["checkpoints"][-1]
        print(case["case"], json.dumps(last), flush=True)
