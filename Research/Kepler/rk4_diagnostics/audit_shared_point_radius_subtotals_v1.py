"""Exact sum of retained native point radii, without field/quadrature work.

A partial prefix is allowed and explicitly reported. Its nonnegative terms
are a lower bound on the eventual *forward error certificate*, not on the
actual physical error. There is no whole-target acceptance verdict here.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import json
from pathlib import Path
import struct
import sys

from audit_shared_regular_q_transfer_v1 import Q, decode_quad, quad_bits, rows, scale_two

ROOT = Path(__file__).resolve().parent
NAMES = ("mass", "px", "py", "angular", "kinetic", "physical_U", "volume",
         "qx", "qy", "polar", "physical_Udot", "potential")
WIDTH = 1008


def add_positive_quad(total, bits):
    assert not bits >> 127
    exponent = (bits >> 112) & 0x7fff
    assert exponent < 0x7fff
    mantissa = bits & ((1 << 112) - 1)
    if exponent:
        mantissa |= 1 << 112
    if not mantissa:
        return total
    power = (exponent - 16383 if exponent else -16382) - 112
    trailing = (mantissa & -mantissa).bit_length() - 1
    mantissa >>= trailing
    power += trailing
    value, unit = total
    if not value:
        return mantissa, power
    if power < unit:
        value <<= unit - power
        unit = power
    return value + (mantissa << (power - unit)), unit


def exact(dyadic):
    return scale_two(Fraction(dyadic[0]), dyadic[1])


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "shared_clipped_target_v8_mpi4"
    run = ROOT / name
    assert run.parent == ROOT and run.is_dir()
    start = datetime.now(timezone.utc).isoformat()
    totals, counts, ranks = {}, Counter(), []
    planned_rows = rows(ROOT / "shared_clipped_target_v7_mpi4/all384_source_completion.csv")
    planned = {int(r["source"]): int(r["planned_panels"]) * 256 for r in planned_rows if int(r["planned_panels"])}
    assert len(planned) == 26 and sum(planned.values()) == 222976
    for rank in range(4):
        packet = run / f"NEW_clipped_positive_queries_rank{rank}.bin"
        count = packet.stat().st_size // WIDTH
        ranks.append(count)
        with packet.open("rb") as handle:
            for first in range(0, count, 4096):
                number = min(4096, count - first)
                raw = handle.read(number * WIDTH)
                assert len(raw) == number * WIDTH
                for i in range(number):
                    offset = i * WIDTH
                    source = struct.unpack_from("<Q", raw, offset)[0]
                    assert source in planned and source % 4 == rank
                    counts[source] += 1
                    total = totals.setdefault(source, [(0, 0)] * 12)
                    for k in range(12):
                        field = offset + 16 + 16 * (38 + k)
                        bits = int.from_bytes(raw[field:field+16], "little")
                        total[k] = add_positive_quad(total[k], bits)
    assert all(n <= planned[s] for s, n in counts.items())
    original = {r["moment"]: r["unchanged_physical_tolerance"] for r in rows(ROOT / "curved_whole_target_v3_mpi4/whole_target_ledger.csv") if r["target"] == "68526"}
    aliases = {"kinetic": "K", "physical_U": "U", "physical_Udot": "thermal_rate"}
    whole = {r["moment"]: r for r in rows(ROOT / "shared_atlas_v4_mpi4/PARTIAL_four_source_all12_errors.csv")}
    native = {}
    native_path = run / "actual_clipped_source_integrals.csv"
    if native_path.exists():
        native = {(int(r["source"]), r["moment"]): r for r in rows(native_path)}
    moment_rows = []
    source_rows = []
    for k, moment in enumerate(NAMES):
        tau = decode_quad(quad_bits(original[aliases.get(moment, moment)]) - 1)
        subtotal = sum((exact(value[k]) for value in totals.values()), Fraction())
        retained = Q(whole[moment]["partial_native_geometry_sum_error"])
        known = subtotal + retained
        ranked = sorted(totals, key=lambda s: exact(totals[s][k]), reverse=True)
        moment_rows.append({
            "moment": moment,
            "new_recorded_point_radii_over_original_tau": float(subtotal / tau),
            "retained_whole4_forward_over_original_tau": float(retained / tau),
            "known_forward_certificate_subtotal_over_original_tau": float(known / tau),
            "known_terms_already_exceed_original_tau": known >= tau,
            "top_recorded_point_radius_sources": [{"source": s, "over_original_tau": float(exact(totals[s][k]) / tau)} for s in ranked[:5]],
        })
        for source, value in sorted(totals.items()):
            item = {"source": source, "moment": moment, "points": counts[source],
                    "exact_radius_integer": str(value[k][0]), "binary_exponent": value[k][1]}
            row = native.get((source, moment))
            if row is not None and row["all12_complete"] == "1":
                assert counts[source] == planned[source]
                tree_error = Q(row["native_geometry_pairwise_sum_error"]) - exact(value[k])
                assert tree_error >= 0
                item["source_tree_rounding_and_bound_propagation_over_original_tau"] = float(tree_error / tau)
            source_rows.append(item)
    print(json.dumps({
        "scope": "saved point-radius accounting only; subtotal of a forward-error certificate, not actual physical error or whole-target publication",
        "run": str(run), "snapshot_start_UTC": start,
        "snapshot_end_UTC": datetime.now(timezone.utc).isoformat(),
        "recorded_prefix_by_rank": ranks, "recorded_prefix_total": sum(ranks),
        "planned_new_records": 222976,
        "all_planned_points_recorded": sum(ranks) == 222976,
        "native_integral_rows_available": bool(native),
        "omitted_nonnegative_terms": ["unrecorded point errors", "new-source tree additions and bound propagation", "analytic remainder", "global addition", "final R80 conversion"],
        "moments": moment_rows, "exact_source_radius_sums": source_rows,
        "whole_target_verdict": "not evaluated by this audit",
    }, indent=2))


if __name__ == "__main__":
    main()
