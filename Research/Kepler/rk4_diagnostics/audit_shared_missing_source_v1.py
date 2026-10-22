"""Read-only exact accounting of reused old25 plus actual missing source."""
from fractions import Fraction
import json
from pathlib import Path
import struct
import sys
from audit_shared_regular_q_transfer_v1 import Q, quad_bits, decode_quad, binary_round, extended, rows, sha

ROOT = Path(__file__).resolve().parent
STATE = "74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "shared_missing_source_v1_mpi1"
    run = ROOT / name
    original = ROOT / "shared_clipped_target_v8_mpi4"
    proof = ROOT / "shared_potential_p2_v1_mpi4"
    assert run.parent == ROOT
    contract = json.loads((run / "contract.json").read_text())
    native_exit = int((run / "native_exit_code.txt").read_text())
    for filename in ("final_source_integrity.txt", "final_input_integrity.txt", "originating_source_integrity.txt",
                     "originating_output_integrity.txt", "potential_proof_source_integrity.txt", "potential_proof_output_integrity.txt"):
        lines = (run / filename).read_text().splitlines()
        assert lines and all(line.endswith(": OK") for line in lines)
    for filename in ("caller_unchanged_state.bin", "immutable_heat_generating_state.bin"):
        assert sha(run / filename) == STATE
    for line in (run / "output.sha256").read_text().splitlines():
        digest, filename = line.split("  ", 1)
        assert sha(Path(filename)) == digest
    completion = rows(run / "all384_source_completion.csv")
    prior_completion = {int(r["source"]): r for r in rows(original / "all384_source_completion.csv")}
    assert len(completion) == 384 and {int(r["source"]) for r in completion} == set(range(384))
    for row in completion:
        source = int(row["source"])
        for field in ("cell", "classification", "coverage_complete", "graph_atlas_complete", "analytic_complete", "retained_group", "planned_panels", "certified_panels", "pullback_AD2_comparisons"):
            assert row[field] == prior_completion[source][field]
        if source != 301:
            assert row == prior_completion[source]
    assert [int(r["source"]) for r in completion if r["retained_group"] == "1"] == [240, 244, 248, 252]
    report = {"scope": "Independent exact saved-record/source/global/serialization accounting; not an independent PDE or derivative evaluation",
              "run": str(run), "native_exit_code": native_exit, "contract": contract,
              "original_current_and_physical_heat_SHA256": STATE,
              "all384_original_coverage_and_old383_source_rows_retained": True}
    panel_log = rows(run / "missing_source_panels.csv")
    assert len(panel_log) == 11 and {int(r["panel"]) for r in panel_log} == set(range(11))
    if not contract["all384_complete"]:
        assert not contract["target_published"] and native_exit == 2
        report["incomplete_actual_panels"] = [r for r in panel_log if r["complete"] != "1"]
        report["target_verdict"] = "unproved; actual missing-source quadrature remains incomplete"
        print(json.dumps(report, indent=2))
        return
    assert all(r["complete"] == "1" and r["points"] == "256" for r in panel_log)
    assert contract["parallel_panel_reduction_matches_serial_saved_moment_tree"]
    assert contract["new_actual_queries"] == 2816 and contract["total_actual_queries"] == 268032
    assert all(r["all12_complete"] == "1" for r in completion)
    current = {(int(r["source"]), r["moment"]): r for r in rows(run / "actual_clipped_source_integrals.csv")}
    previous = {(int(r["source"]), r["moment"]): r for r in rows(original / "actual_clipped_source_integrals.csv")}
    assert len(current) == len(previous) == 312 and set(current) == set(previous)
    patched = {}
    for rank in range(4):
        for row in rows(proof / f"potential_P2_certificates_rank{rank}.csv"):
            source = int(row["source"])
            patched[source] = patched.get(source, Fraction()) + Q(row["chosen_R"])
    for (source, moment), row in current.items():
        assert row["all12_complete"] == "1"
        old = previous[source, moment]
        if source != 301:
            assert row["positive_points"] == old["positive_points"]
            assert Q(row["value"]) == Q(old["value"])
            assert Q(row["native_geometry_pairwise_sum_error"]) == Q(old["native_geometry_pairwise_sum_error"])
        else:
            assert row["positive_points"] == "2816"
        if moment != "potential":
            assert Q(row["analytic_remainder"]) == Q(old["analytic_remainder"])
        else:
            assert patched[source] <= Q(row["analytic_remainder"]) <= Q(old["analytic_remainder"])
    raw = (run / "NEW_missing_source301_queries.bin").read_bytes()
    assert len(raw) == 2816 * 1008
    assert raw == b"".join((run / f"source301_panel{j}.bin").read_bytes() for j in range(11))
    witness = (ROOT / "shared_actual_query_v3_source301_mpi1/diagnostic_positive_query.bin").read_bytes()
    assert raw[16:816] == witness[16:816]
    moment_sum = [Fraction() for _ in range(12)]
    radius_sum = [Fraction() for _ in range(12)]
    for index in range(2816):
        offset = index * 1008
        assert struct.unpack_from("<Q", raw, offset)[0] == 301
        def field(k):
            return decode_quad(int.from_bytes(raw[offset+16+16*k:offset+32+16*k], "little"))
        assert field(2) > 0 and field(26) > 0 and field(31) > 0
        for k in range(12):
            moment_sum[k] += field(26+k)
            radius = field(38+k)
            assert radius >= 0
            radius_sum[k] += radius
    names = ("mass", "px", "py", "angular", "kinetic", "physical_U", "volume", "qx", "qy", "polar", "physical_Udot", "potential")
    for k, moment in enumerate(names):
        row = current[301, moment]
        error = Q(row["native_geometry_pairwise_sum_error"])
        assert error >= radius_sum[k]
        assert abs(Q(row["value"]) - moment_sum[k]) <= error - radius_sum[k]
    report["all2816_positive_records_source_sums_and_forward_radii_accounted"] = True
    report["actual_first_query_fields_bitwise_match_repaired_witness"] = True
    report["old25_native_values_errors_and_other11_remainders_unchanged"] = True
    whole = {r["moment"]: r for r in rows(ROOT / "shared_atlas_v4_mpi4/PARTIAL_four_source_all12_errors.csv")}
    original_tau = {r["moment"]: r["unchanged_physical_tolerance"] for r in rows(ROOT / "curved_whole_target_v3_mpi4/whole_target_ledger.csv") if r["target"] == "68526"}
    aliases = {"kinetic": "K", "physical_U": "U", "physical_Udot": "thermal_rate"}
    caller = (run / "caller_target_packet.bin").read_bytes()
    assert len(caller) == 400
    published = contract["target_published"]
    assert struct.unpack("<QQ", caller[:16]) == ((1, 68526) if published else (9, 1234))
    failures, budgets = [], []
    positive = True
    for k, row in enumerate(rows(run / "all12_ORIGINAL_whole_target_budget.csv")):
        moment = row["moment"]
        assert moment == names[k] and row["all384_complete"] == "1"
        assert bool(int(row["published"])) == published
        selected = [r for (source, name), r in current.items() if name == moment]
        expected = Q(whole[moment]["partial_value"]) + sum((Q(r["value"]) for r in selected), Fraction())
        native_error = Q(whole[moment]["partial_native_geometry_sum_error"]) + sum((Q(r["native_geometry_pairwise_sum_error"]) for r in selected), Fraction())
        remainder_sum = Q(whole[moment]["partial_analytic_remainder"]) + sum((Q(r["analytic_remainder"]) for r in selected), Fraction())
        candidate = Q(row["candidate"]);error = Q(row["completed_native_geometry_sum_error"])
        remainder = Q(row["completed_remainder"]);conversion = Q(row["final_R80_conversion_error"])
        total = Q(row["total_reported_error"])
        tau = decode_quad(quad_bits(original_tau[aliases.get(moment, moment)]) - 1)
        assert Q(row["ORIGINAL_tau"]) == tau
        assert error >= native_error and abs(candidate - expected) <= error - native_error
        assert remainder >= remainder_sum
        rounded = binary_round(candidate, 64, -16382)[0]
        assert conversion == abs(rounded - candidate)
        assert total >= error + remainder + conversion
        if total >= tau:
            failures.append(moment)
        if k in (0, 5):
            positive = positive and rounded > total
        if published:
            assert extended(caller[16+16*k:32+16*k]) == rounded
            assert extended(caller[208+16*k:224+16*k]) == total
        budgets.append({"moment": moment, "reported_error_over_ORIGINAL_tau": float(total/tau),
                        "forward_over_tau": float(error/tau), "analytic_remainder_over_tau": float(remainder/tau),
                        "final_R80_conversion_over_tau": float(conversion/tau)})
    assert len(budgets) == 12 and published == (not failures and positive)
    assert native_exit == (0 if published else 2)
    if not published:
        assert all(extended(caller[16+16*k:32+16*k]) == 0 for k in range(24))
    report["complete_all12_ORIGINAL_post_serialization_budget"] = budgets
    report["ORIGINAL_tau_failures"] = failures
    report["serialized_mass_and_physical_heat_positive"] = positive
    report["target_published"] = published
    report["accepted_RK4_evolution"] = False
    report["production_CUDA_transfer_certified"] = False
    report["long_Kepler_validated"] = False
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
