"""Named frozen-run audit, exact binary128/binary80 accounting; no PDE work."""
import csv
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import struct
import sys

from audit_shared_clipped_target_v1 import ROOT, EXPECTED_STATE, audit, rows, sha


def scale_two(value, exponent):
    return value * (1 << exponent) if exponent >= 0 else value / (1 << -exponent)


def binary_round(value, precision, minimum_exponent):
    """Exact nearest/even binary rounding from an exact rational."""
    if not value:
        return Fraction(), 0, 0
    sign = -1 if value < 0 else 1
    value = abs(value)
    p, q = value.numerator, value.denominator
    exponent = p.bit_length() - q.bit_length()
    if value < scale_two(Fraction(1), exponent):
        exponent -= 1
    unit = max(exponent, minimum_exponent) - precision + 1
    scaled = scale_two(value, -unit)
    mantissa, remainder = divmod(scaled.numerator, scaled.denominator)
    if 2 * remainder > scaled.denominator or (2 * remainder == scaled.denominator and mantissa % 2):
        mantissa += 1
    if mantissa == 1 << precision:
        mantissa >>= 1
        unit += 1
    exact = scale_two(Fraction(sign * mantissa), unit)
    return exact, mantissa, unit


def quad_bits(text):
    number = Decimal(text)
    assert number.is_finite()
    exact, mantissa, unit = binary_round(Fraction(number), 113, -16382)
    if not exact:
        return 0
    exponent = unit + 112
    encoded_exponent = exponent + 16383 if mantissa >= 1 << 112 else 0
    assert 0 <= encoded_exponent < 0x7fff
    fraction = mantissa - (1 << 112) if encoded_exponent else mantissa
    return ((1 if exact < 0 else 0) << 127) | (encoded_exponent << 112) | fraction


def decode_quad(bits):
    exponent = (bits >> 112) & 0x7fff
    mantissa = bits & ((1 << 112) - 1)
    assert exponent < 0x7fff
    if exponent:
        mantissa |= 1 << 112
    return scale_two(Fraction(-mantissa if bits >> 127 else mantissa),
                     (exponent - 16383 if exponent else -16382) - 112)


def Q(text):
    return decode_quad(quad_bits(text))


def extended(raw):
    significand, sign_exponent = struct.unpack("<QH", raw[:10])
    exponent = sign_exponent & 0x7fff
    assert exponent < 0x7fff
    if exponent:
        assert significand >> 63 == 1
    return scale_two(Fraction(-significand if sign_exponent >> 15 else significand),
                     (exponent - 16383 if exponent else -16382) - 63)


def diagnostic_integrity(name):
    path = ROOT / name
    assert int((path / "native_exit_code.txt").read_text()) == 0
    for filename in ("caller_unchanged_state.bin", "immutable_heat_generating_state.bin"):
        assert sha(path / filename) == EXPECTED_STATE
    for filename in ("final_source_integrity.txt", "final_input_integrity.txt"):
        lines = (path / filename).read_text().splitlines()
        assert lines and all(line.endswith(": OK") for line in lines)
    return json.loads((path / "contract.json").read_text())


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "shared_clipped_target_v8_mpi4"
    path = ROOT / name
    report = {
        "scope": "Independent exact arithmetic accounting of one actual mixed-reference target; no quadrature replay or trajectory validation",
        "original_state_SHA256": EXPECTED_STATE,
        "witness29": diagnostic_integrity("shared_clip_witness_v5_mpi4"),
        "reference80_pairwise_diagnostic": diagnostic_integrity("shared_reference80_v2_mpi4"),
        "whole_call": audit(name),
    }
    if not report["whole_call"]["contract"]["target_published"]:
        rejected_packet = (path / "caller_target_packet.bin").read_bytes()
        assert len(rejected_packet) == 400
        assert struct.unpack("<QQ", rejected_packet[:16]) == (9, 1234)
        assert all(extended(rejected_packet[16+16*k:32+16*k]) == 0 for k in range(24))
        report["rejected_caller_header_and_all24_numeric_fields_unchanged"] = True
    certificate, packet_policies = {}, {}
    output_manifest = {}
    for line in (path / "output.sha256").read_text().splitlines():
        digest, filename = line.split("  ", 1)
        output_manifest[Path(filename).name] = digest
    previous_policy = {}
    for rank in range(4):
        for row in rows(ROOT / "shared_clipped_target_v7_mpi4" / f"actual_clipped_panel_certificates_rank{rank}.csv"):
            key = int(row["source"]), int(row["panel"])
            assert key not in previous_policy
            previous_policy[key] = row
    for rank in range(4):
        for row in rows(path / f"actual_clipped_panel_certificates_rank{rank}.csv"):
            key = int(row["source"]), int(row["panel"])
            assert key not in certificate and row["pass"] == "1"
            certificate[key] = row
            original_panel = previous_policy[key]
            for field in ("chart", "whole", "free_lo", "free_hi", "eta_lo", "eta_hi"):
                assert Q(row[field]) == Q(original_panel[field])
            assert (int(row["coordinate_policy"]) == 0) == (original_panel["pass"] == "1")
        filename = path / f"actual_reference_panel_packets_rank{rank}.bin"
        assert sha(filename) == output_manifest[filename.name]
        raw = filename.read_bytes()
        magic, version, epoch, width, count, owner, zero1, zero2 = struct.unpack("<8Q", raw[:64])
        assert (magic, version, epoch, owner, zero1, zero2) == (0x47514650414e4c31, 1, 1, rank, 0, 0)
        assert width % 16 == 0 and len(raw) == 64 + width * count
        for index in range(count):
            start = 64 + index * width
            source, panel, policy, dependent, lower, upper = struct.unpack("<QQiiii", raw[start:start+32])
            assert source % 4 == rank and (source, panel) not in packet_policies
            packet_policies[source, panel] = policy
            assert int(certificate[source, panel]["coordinate_policy"]) == policy
            assert policy in (0, 1, 2)
            if policy == 1:
                assert dependent == 1
            if policy == 2:
                assert dependent == 0
    assert len(certificate) == len(packet_policies) == 871
    assert [sum(p == mode for p in packet_policies.values()) for mode in (0, 1, 2)] == [842, 27, 2]
    report["exact_typed_panel_packets_and_fixed_coordinate_policy"] = {"panels": 871, "original_Duffy": 842, "regular_q": 27, "radial_prior": 2}

    contributions = {}
    for row in rows(path / "actual_clipped_source_integrals.csv"):
        key = int(row["source"]), row["moment"]
        assert key not in contributions and row["all12_complete"] in ("0", "1")
        contributions[key] = row
    sources = {s for s, _ in contributions}
    assert len(contributions) == 312 and len(sources) == 26
    for (source, moment), row in contributions.items():
        actual = Q(row["analytic_remainder"])
        raw_sum = sum((Q(c[f"{moment}_remainder"]) for (s, _), c in certificate.items() if s == source), Fraction())
        assert actual >= raw_sum
    incomplete = sorted({s for (s, _), row in contributions.items() if row["all12_complete"] != "1"})
    if incomplete:
        assert not report["whole_call"]["contract"]["target_published"]
        assert report["whole_call"]["native_exit_code"] == 2
        report["incomplete_integral_sources"] = [row for row in rows(path / "all384_source_completion.csv") if int(row["source"]) in incomplete]
        report["all871_analytic_remainders_retained_and_accounted"] = True
        report["whole_target_published"] = False
        report["original_twelve_tau_status"] = "unproved because one or more actual integrals are incomplete"
        whole = {r["moment"]: r for r in rows(ROOT / "shared_atlas_v4_mpi4/PARTIAL_four_source_all12_errors.csv")}
        original = {r["moment"]: r["unchanged_physical_tolerance"] for r in rows(ROOT / "curved_whole_target_v3_mpi4/whole_target_ledger.csv") if r["target"] == "68526"}
        aliases = {"kinetic": "K", "physical_U": "U", "physical_Udot": "thermal_rate"}
        partial = []
        for moment in whole:
            tau = decode_quad(quad_bits(original[aliases.get(moment, moment)]) - 1)
            complete = [row for (s, k), row in contributions.items() if k == moment and row["all12_complete"] == "1"]
            all_sources = [row for (s, k), row in contributions.items() if k == moment]
            retained_forward = Q(whole[moment]["partial_native_geometry_sum_error"])
            new_forward = sum((Q(row["native_geometry_pairwise_sum_error"]) for row in complete), Fraction())
            analytic = Q(whole[moment]["partial_analytic_remainder"]) + sum((Q(row["analytic_remainder"]) for row in all_sources), Fraction())
            ranked = sorted(complete, key=lambda row: Q(row["native_geometry_pairwise_sum_error"]), reverse=True)
            partial.append({
                "moment": moment,
                "scope": "known error terms only; excludes missing forward contribution, global addition and final conversion; no completed target verdict",
                "completed_new_sources": len(complete),
                "retained_whole4_forward_over_original_tau": float(retained_forward / tau),
                "completed_new_forward_over_original_tau": float(new_forward / tau),
                "all871_plus_retained_analytic_remainder_over_original_tau": float(analytic / tau),
                "known_terms_over_original_tau": float((retained_forward + new_forward + analytic) / tau),
                "top_completed_new_forward_sources": [{"source": int(row["source"]), "over_original_tau": float(Q(row["native_geometry_pairwise_sum_error"]) / tau)} for row in ranked[:5]],
            })
        report["incomplete_call_known_error_terms"] = partial
        report["accepted_RK4_evolution"] = False
        report["production_CUDA_transfer_certified"] = False
        report["long_Kepler_validated"] = False
        print(json.dumps(report, indent=2))
        return
    whole = {r["moment"]: r for r in rows(ROOT / "shared_atlas_v4_mpi4/PARTIAL_four_source_all12_errors.csv")}
    original = {}
    for row in rows(ROOT / "curved_whole_target_v3_mpi4/whole_target_ledger.csv"):
        if row["target"] == "68526":
            original[row["moment"]] = row["unchanged_physical_tolerance"]
    aliases = {"kinetic": "K", "physical_U": "U", "physical_Udot": "thermal_rate"}
    budget = rows(path / "all12_ORIGINAL_whole_target_budget.csv")
    assert len(budget) == 12
    caller = (path / "caller_target_packet.bin").read_bytes()
    assert len(caller) == 400
    published = report["whole_call"]["contract"]["target_published"]
    violations = []
    breakdown = []
    positive = True
    for k, row in enumerate(budget):
        moment = row["moment"]
        assert row["all384_complete"] == "1" and bool(int(row["published"])) == published
        tau_bits = quad_bits(original[aliases.get(moment, moment)])
        assert Q(row["ORIGINAL_tau"]) == decode_quad(tau_bits - 1)
        selected = [contributions[s, moment] for s in sorted(sources)]
        sum_value = Q(whole[moment]["partial_value"]) + sum((Q(r["value"]) for r in selected), Fraction())
        sum_error = Q(whole[moment]["partial_native_geometry_sum_error"]) + sum((Q(r["native_geometry_pairwise_sum_error"]) for r in selected), Fraction())
        sum_remainder = Q(whole[moment]["partial_analytic_remainder"]) + sum((Q(r["analytic_remainder"]) for r in selected), Fraction())
        candidate = Q(row["candidate"])
        error = Q(row["completed_native_geometry_sum_error"])
        remainder = Q(row["completed_remainder"])
        conversion = Q(row["final_R80_conversion_error"])
        reported = Q(row["total_reported_error"])
        assert error >= sum_error and abs(candidate - sum_value) <= error - sum_error
        assert remainder >= sum_remainder
        rounded80 = binary_round(candidate, 64, -16382)[0]
        assert conversion == abs(rounded80 - candidate)
        assert reported >= error + remainder + conversion
        tau = Q(row["ORIGINAL_tau"])
        if reported >= tau:
            violations.append({"moment": moment, "reported_error_over_original_tau": float(reported / Q(row["ORIGINAL_tau"])), "forward_error_over_tau": float(error / Q(row["ORIGINAL_tau"])), "analytic_remainder_over_tau": float(remainder / Q(row["ORIGINAL_tau"]))})
        if k in (0, 5):
            positive = positive and rounded80 > reported
        ranked_forward = sorted(selected, key=lambda r: Q(r["native_geometry_pairwise_sum_error"]), reverse=True)
        ranked_remainder = sorted(selected, key=lambda r: Q(r["analytic_remainder"]), reverse=True)
        breakdown.append({
            "moment": moment,
            "reported_error_over_original_tau": float(reported / tau),
            "retained_whole4_forward_over_tau": float(Q(whole[moment]["partial_native_geometry_sum_error"]) / tau),
            "new26_forward_over_tau": float(sum((Q(r["native_geometry_pairwise_sum_error"]) for r in selected), Fraction()) / tau),
            "global_Q128_forward_over_tau": float((error - sum_error) / tau),
            "retained_whole4_analytic_remainder_over_tau": float(Q(whole[moment]["partial_analytic_remainder"]) / tau),
            "new26_analytic_remainder_over_tau": float(sum((Q(r["analytic_remainder"]) for r in selected), Fraction()) / tau),
            "final_R80_conversion_over_tau": float(conversion / tau),
            "top_new_forward_sources": [{"source": int(r["source"]), "over_original_tau": float(Q(r["native_geometry_pairwise_sum_error"]) / tau)} for r in ranked_forward[:5]],
            "top_new_analytic_remainder_sources": [{"source": int(r["source"]), "over_original_tau": float(Q(r["analytic_remainder"]) / tau)} for r in ranked_remainder[:5]],
        })
        if published:
            assert extended(caller[16+16*k:32+16*k]) == rounded80
            assert extended(caller[208+16*k:224+16*k]) == reported
    assert published == (len(violations) == 0 and positive)
    assert report["whole_call"]["native_exit_code"] == (0 if published else 2)
    report["exact_source_group_global_sum_and_serialization_accounting"] = True
    report["original_twelve_tau_failures"] = violations
    report["saved_contribution_budget_breakdown"] = breakdown
    report["serialized_mass_and_physical_heat_strictly_positive"] = positive
    report["whole_target_published"] = published

    # These are independent samples of the retained packet representation;
    # the actual per-query field/all12 oracle bounds were checked natively.
    sampled = 0
    total = 0
    for rank in range(4):
        filename = path / f"NEW_clipped_positive_queries_rank{rank}.bin"
        count, tail = divmod(filename.stat().st_size, 1008)
        assert tail == 0
        total += count
        if not count:
            continue
        with filename.open("rb") as handle:
            for index in sorted({0, count-1, *(i*count//8 for i in range(8))}):
                handle.seek(index*1008)
                raw = handle.read(1008)
                source = struct.unpack("<Q", raw[:8])[0]
                assert source in sources and source % 4 == rank
                values = [decode_quad(int.from_bytes(raw[16+16*j:32+16*j], "little")) for j in range(62)]
                assert values[2] > 0 and values[26] > 0 and values[31] > 0
                assert all(v >= 0 for v in values[3:26])
                assert all(v >= 0 for v in values[38:62])
                for j in range(26, 50):
                    bits = int.from_bytes(raw[16+16*j:32+16*j], "little")
                    assert bits & ((1 << 49)-1) == 0
                sampled += 1
    assert total == 222976
    report["new_positive_query_records"] = total
    report["independently_sampled_packet_records"] = sampled
    report["accepted_RK4_evolution"] = False
    report["production_CUDA_transfer_certified"] = False
    report["long_Kepler_validated"] = False
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
