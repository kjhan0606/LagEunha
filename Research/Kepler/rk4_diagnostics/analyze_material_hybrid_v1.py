#!/usr/bin/env python3
"""Small read-only native-result comparison; never a hydro evolution path."""
import argparse
import csv
import json
import math
from pathlib import Path
import struct


def state(path):
    data = path.read_bytes()
    magic, nb, nr = struct.unpack_from("=QQQ", data)
    assert magic == 0x4746485942523031
    count = 6 * nb + 2 * nr
    assert len(data) == 24 + 8 * count
    values = struct.unpack_from(f"={count}d", data, 24)
    out, offset = {}, 0
    for name, size in [(k, nb) for k in ("x", "y", "vx", "vy", "ax", "ay")] + [("lambda_x", nr), ("lambda_y", nr)]:
        out[name] = values[offset:offset + size]
        offset += size
    return out


def traction(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def comparisons(a, b, root_a, root_b):
    out = {}
    for ca, cb in zip(a["cases"], b["cases"]):
        assert ca["case"] == cb["case"]
        name = ca["case"]
        sa, sb = state(root_a / f"{name}_state.bin"), state(root_b / f"{name}_state.bin")
        fields = {k: max(abs(x - y) for x, y in zip(sa[k], sb[k])) for k in sa}
        ta, tb = traction(root_a / f"{name}_traction.csv"), traction(root_b / f"{name}_traction.csv")
        tdiff = None
        if len(ta) == len(tb):
            assert all((x["face"], x["q"]) == (y["face"], y["q"]) for x, y in zip(ta, tb))
            tdiff = max(abs(float(x[k]) - float(y[k])) for x, y in zip(ta, tb) for k in ("tx_reference", "ty_reference"))
        metrics = ["max_numerical_energy_residual_over_U", "max_numerical_angular_residual_over_2K", "max_physical_energy_residual_over_U", "max_physical_angular_residual_over_2K", "traction_difference_L2_over_P"]
        local_metrics = ["acceleration_error_RMS", "gravity_RMS", "acceleration_error_max", "independent_source_projection_work_over_U", "independent_source_projection_torque_over_2K", "pressure_material_curvature_error_L1_volume", "pressure_material_curvature_error_max", "pressure_material_curvature_exact_max", "pressure_material_curvature_numerical_max"]
        local = {}
        assert [c["id"] for c in ca["cells"]] == [c["id"] for c in cb["cells"]]
        for key in local_metrics:
            if key in ca["cells"][0] and key in cb["cells"][0]:
                local[key] = max(abs(x[key] - y[key]) for x, y in zip(ca["cells"], cb["cells"]))
        out[name] = {"state_max_absolute_differences": fields, "sampled_traction_max_absolute_difference": tdiff, "metric_absolute_differences": {k: abs(ca[k] - cb[k]) for k in metrics}, "cell_metric_max_absolute_differences": local}
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("versions", nargs="+", help="material_hybrid_v1 etc., absolute paths")
    parser.add_argument("--output", required=True)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    output = Path(args.output)
    assert not output.exists(), "preserve previous analysis"
    report = {"scope": "initial material-hybrid operator reference only; numerical-traction closure is not physical-traction accuracy or long-time stability", "versions": {}, "quadrature_comparison": {}}
    docs = []
    for arg in args.versions:
        root = Path(arg)
        assert (root / "status.txt").read_text().strip() == "AUDIT_COMPLETE"
        a = json.loads((root / "mpi1/contract.json").read_text())
        b = json.loads((root / "mpi4/contract.json").read_text())
        assert len(a["cases"]) == len(b["cases"]) == 4
        for c in a["cases"] + b["cases"]:
            assert c["old_space_not_in_kernel_rows"] == 0
            assert c["disconnected_couplings"] == 0
            assert c["conforming_nullity"] == c["broken_scalar_dofs"] - c["trace_rank"]
            assert c["Schur_min_eigenvalue"] > 0
            assert all(math.isfinite(x) for x in c.values() if isinstance(x, (float, int)))
        report["versions"][root.name] = {"mpi1": a, "mpi4": b, "mpi_comparison": comparisons(a, b, root / "mpi1", root / "mpi4")}
        docs.append((root, a))
    for (ra, a), (rb, b) in zip(docs, docs[1:]):
        report["quadrature_comparison"][f"{ra.name}_to_{rb.name}"] = comparisons(a, b, ra / "mpi1", rb / "mpi1")
    output.write_text(json.dumps(report, indent=2) + "\n")
    for version, doc in ([] if args.quiet else report["versions"].items()):
        for c in doc["mpi1"]["cases"]:
            print(version, c["case"], json.dumps({k: c[k] for k in ("conforming_nullity", "extra_conforming_modes", "Schur_condition", "trace_Gram_condition_max", "momentum_equation_relative", "direct_traction_load_relative", "old_residual_kernel_projection_relative", "acceleration_vs_old_max", "max_numerical_energy_residual_over_U", "max_numerical_angular_residual_over_2K", "max_physical_energy_residual_over_U", "traction_difference_L2_over_P")}))
    print("MPI and quadrature differences saved:", output)


if __name__ == "__main__":
    main()
