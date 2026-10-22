"""Exact saved-certificate arithmetic; no physical query evaluation."""
from fractions import Fraction
import json
from math import factorial
from pathlib import Path
import struct
import sys
from audit_shared_regular_q_transfer_v1 import Q, quad_bits, decode_quad, rows, sha

ROOT = Path(__file__).resolve().parent
STATE = "74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005"


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "shared_potential_p2_v1_mpi4"
    run = ROOT / name
    original_run = ROOT / "shared_clipped_target_v8_mpi4"
    assert run.parent == ROOT
    contract = json.loads((run / "contract.json").read_text())
    assert contract["original_panels"] == 871 and not contract["target_published"]
    for filename in ("final_source_integrity.txt", "final_input_integrity.txt",
                     "originating_source_integrity.txt", "originating_output_integrity.txt",
                     "originating_input_integrity.txt"):
        lines = (run / filename).read_text().splitlines()
        assert lines and all(line.endswith(": OK") for line in lines)
    for filename in ("caller_unchanged_state.bin", "immutable_heat_generating_state.bin"):
        assert sha(run / filename) == STATE
    for line in (run / "output.sha256").read_text().splitlines():
        digest, filename = line.split("  ", 1)
        assert sha(Path(filename)) == digest
    old, geometry = {}, {}
    for owner in range(4):
        for row in rows(original_run / f"actual_clipped_panel_certificates_rank{owner}.csv"):
            key = int(row["source"]), int(row["panel"])
            assert key not in old and row["pass"] == "1"
            old[key] = row
        raw = (original_run / f"actual_reference_panel_packets_rank{owner}.bin").read_bytes()
        magic, version, epoch, width, count, rank, z1, z2 = struct.unpack("<8Q", raw[:64])
        assert (magic, version, epoch, rank, z1, z2) == (0x47514650414e4c31, 1, 1, owner, 0, 0)
        assert len(raw) == 64 + width * count
        for index in range(count):
            key_s, key_p, policy, dependent, lower, upper = struct.unpack_from("<QQiiii", raw, 64 + width * index)
            geometry[key_s, key_p] = owner, index, policy, dependent, lower, upper
    exact_C = Fraction(2, 33) * Fraction(factorial(16)**2, factorial(32))**2 * Fraction(4, 3)**32
    chosen, proved, improved = {}, 0, 0
    attempted = 0
    for rank in range(contract["MPI_ranks"]):
        for row in rows(run / f"potential_P2_certificates_rank{rank}.csv"):
            key = int(row["source"]), int(row["panel"])
            assert key not in chosen and key[0] % contract["MPI_ranks"] == rank
            recorded = tuple(int(row[field]) for field in ("origin_rank", "origin_record", "coordinate_policy", "dependent", "lower", "upper"))
            assert recorded == geometry[key]
            old_R = Q(row["old_R"])
            assert old_R == Q(old[key]["potential_remainder"])
            selected = Q(row["chosen_R"])
            valid = row["new_proof"] == "1"
            attempted += int(row["AD2_evaluations"])
            if valid:
                assert row["center_complete"] == row["axis0_complete"] == row["axis1_complete"] == "1"
                majorants = []
                for axis in range(2):
                    radius = (Fraction(5, 4), Fraction(1, 2)) if axis == 0 else (Fraction(1, 2), Fraction(5, 4))
                    bound = Q(row["center_value_error"])
                    bound += sum(Q(row[f"center_grad{i}_error"]) * radius[i] for i in range(2))
                    bound += sum(Fraction(1, 2) * Q(row[f"Hdiff{axis}_{i}{j}"]) * radius[i] * radius[j] for i in range(2) for j in range(2))
                    actual = Q(row[f"axis{axis}_majorant"])
                    assert actual >= bound >= 0
                    majorants.append(actual)
                new_R = Q(row["new_R"])
                assert new_R >= exact_C * max(majorants)
                assert selected == min(old_R, new_R)
                proved += 1
            else:
                assert row["new_R"] == "unproved" and selected == old_R
            improved += selected < old_R
            chosen[key] = selected
    assert len(chosen) == len(old) == len(geometry) == 871
    assert proved == contract["new_valid_P2_proofs"] and improved == contract["strictly_improved_potential_remainders"]
    assert attempted == contract["full_potential_AD2_evaluations_attempted"]
    whole = {r["moment"]: r for r in rows(ROOT / "shared_atlas_v4_mpi4/PARTIAL_four_source_all12_errors.csv")}
    original_tau = {r["moment"]: r["unchanged_physical_tolerance"] for r in rows(ROOT / "curved_whole_target_v3_mpi4/whole_target_ledger.csv") if r["target"] == "68526"}
    aliases = {"kinetic": "K", "physical_U": "U", "physical_Udot": "thermal_rate"}
    native = rows(original_run / "actual_clipped_source_integrals.csv")
    report = []
    possible = True
    for row in rows(run / "known_all12_budget.csv"):
        moment = row["moment"]
        selected = [x for x in native if x["moment"] == moment]
        assert len(selected) == 26
        forward_min = Q(whole[moment]["partial_native_geometry_sum_error"])
        forward_min += sum((Q(x["native_geometry_pairwise_sum_error"]) for x in selected if x["all12_complete"] == "1"), Fraction())
        remainder_min = Q(whole[moment]["partial_analytic_remainder"])
        remainder_min += sum(chosen.values()) if moment == "potential" else sum((Q(x["analytic_remainder"]) for x in selected), Fraction())
        forward = Q(row["completed_source_forward_and_Qsum_error"])
        remainder = Q(row["all_source_analytic_remainder"])
        known = Q(row["known_total_without_missing_source_or_final_conversion"])
        tau = decode_quad(quad_bits(original_tau[aliases.get(moment, moment)]) - 1)
        assert Q(row["ORIGINAL_tau"]) == tau
        assert forward >= forward_min and remainder >= remainder_min and known >= forward + remainder
        ok = known < tau
        assert bool(int(row["known_terms_strictly_below_tau"])) == ok
        possible = possible and ok
        report.append({"moment": moment, "known_forward_over_original_tau": float(forward/tau),
                       "all_source_analytic_remainder_over_original_tau": float(remainder/tau),
                       "known_total_over_original_tau": float(known/tau)})
    assert len(report) == 12 and possible == contract["known_budget_compatible"]
    assert int((run / "native_exit_code.txt").read_text()) == (0 if possible else 2)
    print(json.dumps({"scope": "Exact arithmetic/provenance audit of saved fixed-P2 certificates, not an independent derivative evaluation or completed target publication",
                      "run": str(run), "contract": contract, "exact_Cauchy_coefficient_and_Hessian_integral_accounting": True,
                      "all871_typed_domain_branch_bindings_match": True, "original_state_and_heat_SHA256": STATE,
                      "known_budget_compatible": possible, "known_all12": report,
                      "missing_source301_forward_and_final_serialization_unproved": True,
                      "target_published": False}, indent=2))


if __name__ == "__main__":
    main()
