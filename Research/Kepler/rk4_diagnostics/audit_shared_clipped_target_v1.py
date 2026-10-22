"""Read-only audit of named frozen runs. Does not run field/quadrature code."""
import collections
import csv
import hashlib
import json
from pathlib import Path
import struct
import sys

ROOT = Path(__file__).resolve().parent
EXPECTED_STATE = "74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(name):
    path = ROOT / name
    assert path.parent == ROOT and path.is_dir()
    native = int((path / "native_exit_code.txt").read_text())
    contract = json.loads((path / "contract.json").read_text())
    report = {"run": str(path), "native_exit_code": native, "contract": contract}
    checkpoints = {}
    for filename in ("caller_unchanged_state.bin", "immutable_heat_generating_state.bin"):
        if (path / filename).is_file():
            checkpoints[filename] = sha(path / filename)
            assert checkpoints[filename] == EXPECTED_STATE
    report["checkpoint_SHA256"] = checkpoints
    integrity = {}
    for filename in ("final_source_integrity.txt", "final_input_integrity.txt"):
        lines = (path / filename).read_text().splitlines()
        assert lines and all(line.endswith(": OK") for line in lines)
        integrity[filename] = len(lines)
    report["runner_final_integrity_checks"] = integrity
    manifest = {}
    for line in (path / "output.sha256").read_text().splitlines():
        digest, filename = line.split("  ", 1)
        manifest[Path(filename).name] = digest
    # Verify small scientific ledgers/checkpoints now; large query packets
    # were hashed by the terminal runner, and are counted here by byte size.
    for filename, digest in manifest.items():
        if filename.endswith((".json", ".csv")) or filename in checkpoints or filename == "caller_target_packet.bin":
            assert sha(path / filename) == digest
    report["small_output_manifest_rechecked"] = True
    source_path = path / "all384_source_completion.csv"
    if source_path.is_file():
        source = rows(source_path)
        assert len(source) == 384 and {int(r["source"]) for r in source} == set(range(384))
        assert all(r["coverage_complete"] == "1" for r in source)
        report["source_reasons"] = dict(collections.Counter(r["reason"] for r in source))
        report["plan_unproved_source_IDs"] = [int(r["source"]) for r in source if r["graph_atlas_complete"] != "1"]
        report["new_planned_panels"] = sum(int(r["planned_panels"]) for r in source)
        report["new_certified_panels"] = sum(int(r["certified_panels"]) for r in source)
        report["new_positive_queries"] = sum(int(r["new_positive_points"]) for r in source)
        report["retained_whole_source_IDs"] = [int(r["source"]) for r in source if r["retained_group"] == "1"]
        assert report["retained_whole_source_IDs"] == [240, 244, 248, 252]
        chart = []
        for rank in range(contract["MPI_ranks"]):
            chart.extend(rows(path / f"actual_clipped_charts_rank{rank}.csv"))
        assert len(chart) == 156
        report["chart_reasons"] = dict(collections.Counter(r["reason"] for r in chart))
        report["chart_unproved_source_chart_IDs"] = [[int(r["source"]), int(r["chart"])] for r in chart if r["complete"] != "1"]
        certificate_paths = [path / f"actual_clipped_panel_certificates_rank{rank}.csv" for rank in range(contract["MPI_ranks"])]
        if all(p.is_file() for p in certificate_paths):
            certificate = [r for p in certificate_paths for r in rows(p)]
            report["analytic_panels_evaluated"] = len(certificate)
            report["analytic_panels_passed"] = sum(r["pass"] == "1" for r in certificate)
            report["analytic_panel_failure_reasons"] = dict(collections.Counter(r["reason"] for r in certificate if r["pass"] != "1"))
            report["analytic_unproved_source_panel_IDs"] = [[int(r["source"]), int(r["panel"])] for r in certificate if r["pass"] != "1"]
            assert report["analytic_panels_passed"] == report["new_certified_panels"]
            for r in certificate:
                if r["pass"] == "1" and r["implicit_observed"] == "1":
                    assert float(r["contraction"]) < 1 and float(r["invariance"]) < 1
                    if "prior_radius" in r:
                        assert float(r["real_root_radius"]) <= float(r["correction_radius"]) <= float(r["prior_radius"])
        if "pullback_AD2_comparisons" in source[0]:
            report["pullback_AD2_comparisons"] = sum(int(r["pullback_AD2_comparisons"]) for r in source)
            assert report["pullback_AD2_comparisons"] == 21528
        report["caller_epoch_target"] = list(struct.unpack("<QQ", (path / "caller_target_packet.bin").read_bytes()[:16]))
        if not contract["target_published"]:
            assert report["caller_epoch_target"] == [9, 1234]
        if any(r["all12_complete"] != "1" for r in source):
            assert all(r["total_reported_error"] == "unproved" for r in rows(path / "all12_ORIGINAL_whole_target_budget.csv"))
        record_bytes = contract.get("reference80_record_bytes", 448)
        sizes = [(path / f"NEW_clipped_positive_queries_rank{rank}.bin").stat().st_size for rank in range(contract["MPI_ranks"])]
        assert all(size % record_bytes == 0 for size in sizes)
        assert sum(sizes) // record_bytes == report["new_positive_queries"]
        report["new_query_binary_bytes"] = sum(sizes)
        report["new_query_record_bytes"] = record_bytes
    return report


def main():
    coverage = rows(ROOT / "shared_boundary_witness_v5_mpi1/source_coverage.csv")
    assert len(coverage) == 384 and all(r["coverage_complete"] == "1" for r in coverage)
    classes = collections.Counter(int(r["classification"]) for r in coverage)
    assert classes == {0: 344, 1: 4, 2: 36}
    empty = sum(r["classification"] == "2" and r["empty"] == "1" for r in coverage)
    assert empty == 10
    whole = rows(ROOT / "shared_atlas_v4_mpi4/whole_inside_source_integrals.csv")
    assert {int(r["source"]) for r in whole} == {240, 244, 248, 252}
    assert all(r["all12_complete"] == "1" for r in whole)
    assert sum(int(r["panels"]) for r in whole) == sum(int(r["certified_panels"]) for r in whole) == 176
    assert sum(int(r["positive_points"]) for r in whole) == 45056
    whole_bytes = sum((ROOT / "shared_atlas_v4_mpi4" / f"PARTIAL_whole_inside_queries_rank{rank}.bin").stat().st_size for rank in range(4))
    assert whole_bytes == 45056 * 448
    report = {
        "scope": "Read-only frozen coverage, retained whole4 and named actual clipped transaction audit; no trajectory validation",
        "state_SHA256": EXPECTED_STATE,
        "coverage_all384_complete": True,
        "classification_counts": dict(classes),
        "candidate_empty_complete": empty,
        "actual_clipped_sources": 26,
        "retained_whole_source_IDs": [240, 244, 248, 252],
        "retained_whole_panels": 176,
        "retained_whole_positive_queries": 45056,
        "runs": [audit(name) for name in sys.argv[1:]],
        "accepted_RK4_evolution": False,
        "CUDA_production_transfer_certified": False,
        "long_Kepler_validated": False,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
