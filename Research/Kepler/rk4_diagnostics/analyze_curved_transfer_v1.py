#!/usr/bin/env python3
"""Read native transport/affinity artifacts; no hydro or time stepping."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import struct


FIELDS = ("mass", "px", "py", "angular", "K", "U", "volume", "qx", "qy", "polar", "thermal_rate", "potential")


def moments(path):
    rows = list(csv.DictReader(path.open()))
    source, target = {}, {}
    for row in rows:
        i, j = int(row["source"]), int(row["target"])
        source.setdefault(i, {k: [] for k in FIELDS})
        target.setdefault(j, {k: [] for k in FIELDS})
        for k in FIELDS:
            source[i][k].append(float(row[k]))
            target[j][k].append(float(row[k]))
    return ({i: {k: math.fsum(v) for k, v in d.items()} for i, d in group.items()} for group in (source, target))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("directory", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    root = args.directory
    assert (root / "status.txt").read_text().strip() == "AUDIT_COMPLETE"
    native = {mode: json.loads((root / mode / "contract.json").read_text()) for mode in ("coarse_mpi1", "coarse_mpi4", "fine_mpi4", "weighted_mpi4")}
    assert all(d["budget_failures"] == 0 for d in native.values())
    assert (root / "coarse_mpi1/transfer_moments.csv").read_bytes() == (root / "coarse_mpi4/transfer_moments.csv").read_bytes()
    assert (root / "coarse_mpi1/angular_uncertainty.csv").read_bytes() == (root / "coarse_mpi4/angular_uncertainty.csv").read_bytes()
    reference = json.loads((root.parent / "material_p3_finite_v2/mpi1_cpu/contract.json").read_text())["cases"][1]["cells"]
    reference = {c["id"]: c for c in reference}
    parsed = {mode: tuple(moments(root / mode / "transfer_moments.csv")) for mode in native}
    report = {"scope": "native curved quadrature estimates; target radial convergence fails; no certified remap, weight Newton, compression or evolution", "native": native, "MPI_coarse_bitmatch": True, "modes": {}}
    for mode, (source, target) in parsed.items():
        report["modes"][mode] = {
            "source_mass_relative_error_max": max(abs(source[i]["mass"] / c["mass"] - 1) for i, c in reference.items()),
            "source_K_error_over_U_max": max(abs(source[i]["K"] - c["K"]) / c["U"] for i, c in reference.items()),
            "source_U_relative_error_max": max(abs(source[i]["U"] / c["U"] - 1) for i, c in reference.items()),
            "target_W0_mass_change_or_pattern_mass_change_max": max(abs(target[i]["mass"] / c["mass"] - 1) for i, c in reference.items()),
            "global_source_target_conservation": {k: abs(math.fsum(source[i][k] for i in source) - math.fsum(target[i][k] for i in target)) for k in FIELDS},
            "target_moments": target,
        }
    _, coarse = parsed["coarse_mpi4"]
    _, fine = parsed["fine_mpi4"]
    report["radial_quadrature_sensitivity"] = {
        "target_mass_difference_over_original_mass_max": max(abs(coarse[i]["mass"] - fine[i]["mass"]) / c["mass"] for i, c in reference.items()),
        "target_K_difference_over_original_U_max": max(abs(coarse[i]["K"] - fine[i]["K"]) / c["U"] for i, c in reference.items()),
        "target_U_difference_over_original_U_max": max(abs(coarse[i]["U"] - fine[i]["U"]) / c["U"] for i, c in reference.items()),
        "cold_transfer_gate": "FAIL; source/global conservation does not certify target allocation",
    }
    ids = sorted(reference)
    mass = {i: [] for i in ids}
    point_count = 0
    packet_integrity = []
    for rank in range(4):
        path = root / "fine_mpi4" / f"positive_quadrature_rank{rank}.bin"
        raw = path.read_bytes()
        h = struct.unpack_from("=6Q", raw)
        assert h[:4] == (0x474650334F564C31, 1, rank, 4) and h[5] == 56
        assert len(raw) == 48 + h[4] * 56
        for tri, target, flags, x, y, w, rho in struct.iter_unpack("=3Q4d", raw[48:]):
            assert tri < 384 and target < 24 and flags in (0, 1, 2)
            assert all(math.isfinite(v) for v in (x, y, w, rho)) and w > 0 and rho > 0
            mass[ids[target]].append(w * rho)
            point_count += 1
        packet_integrity.append({"rank": rank, "points": h[4], "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()})
    report["positive_quadrature_packets"] = packet_integrity
    report["saved_points"] = point_count
    report["saved_packet_target_mass_roundtrip_relative_max"] = max(abs(math.fsum(mass[i]) - fine[i]["mass"]) / reference[i]["mass"] for i in ids)
    trace = json.loads((root / "fine_mpi4/whole_trace_bounds.json").read_text())
    report["whole_edge_trace_bounds"] = trace
    affinity = root.parent / "curved_affinity_v1"
    old = native["coarse_mpi4"]["seconds"]
    new = json.loads((affinity / "transfer_mpi4/contract.json").read_text())["seconds"]
    assert (affinity / "status.txt").read_text().strip() == "AUDIT_COMPLETE"
    assert (affinity / "transfer_mpi4/transfer_moments.csv").read_bytes() == (root / "coarse_mpi4/transfer_moments.csv").read_bytes()
    report["affinity_fix"] = {"old_seconds": old, "new_seconds": new, "speedup": old / new, "transport_CSV_bitmatch": True, "old_mapping": (affinity / "inherited.txt").read_text(), "new_mapping": (affinity / "unbound.txt").read_text()}
    with args.output.open("x") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps({"source_fine": {k: v for k, v in report["modes"]["fine_mpi4"].items() if k not in ("target_moments", "global_source_target_conservation")}, "radial_sensitivity": report["radial_quadrature_sensitivity"], "speedup": report["affinity_fix"]["speedup"], "positive_points": point_count, "packet_mass_roundtrip_relative": report["saved_packet_target_mass_roundtrip_relative_max"]}, indent=2))


if __name__ == "__main__":
    main()
