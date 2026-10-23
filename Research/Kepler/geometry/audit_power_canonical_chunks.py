#!/usr/bin/env python3
"""Independently compare sealed canonical events with every incidence wire."""

import argparse
import csv
from decimal import Decimal
from pathlib import Path
import re

from audit_power_edge_chunks import (
    audit_chunk,
    check_sha_file,
    check_source_snapshot,
    require,
)

KEY = (
    "event_edge",
    "event_cell_low",
    "event_cell_high",
    "event_ordinal",
    "target_segment",
)
SIMPLE = (
    "event_derivative_sign",
    "target_endpoint",
    "source_vertex",
    "edge_owner",
)
INCIDENCE = (
    "source_rank",
    "incident_triangle",
    "local_edge",
    "orientation",
)


def event_key(row):
    return tuple(int(row[field]) for field in KEY)


def audit_one(chunk):
    first, last, total, seal, ranks = audit_chunk(chunk)
    variant = "_root60_v1" if "_root60" in chunk.name else "_v2"
    canonical = chunk.parent / f"power_edge_reconcile_{first}_{last}{variant}"
    require((canonical / "native_exit_code.txt").read_text().strip() == "0",
            f"native reconciliation did not exit zero: {canonical}")
    check_source_snapshot(canonical)
    check_sha_file(canonical / "binary.sha256")
    check_sha_file(canonical / "canonical.sha256")
    require("FAILED" not in (canonical / "final_source_check.txt").read_text(),
            f"reconciliation input changed during run: {canonical}")

    groups = {}
    wires = 0
    for rank in range(4):
        with (chunk / f"source_edge_wires_rank{rank}.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                key = event_key(row)
                require(first <= key[0] < last, f"wire outside edge range: {key}")
                groups.setdefault(key, []).append(row)
                wires += 1

    log = (canonical / "native.log").read_text()
    owner_counts = []
    events = 0
    for owner in range(4):
        found = re.findall(
            rf"NATIVE_CANONICAL_RECONCILE_PASS owner={owner} events=(\d+)", log
        )
        require(len(found) == 1, f"missing unique native owner receipt: {canonical}")
        declared = int(found[0])
        actual = 0
        with (canonical / f"canonical_owner{owner}.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                key = event_key(row)
                require(key in groups, f"canonical event without incidence wires: {key}")
                incidence = groups.pop(key)
                require(int(row["edge_owner"]) == owner and
                        len(incidence) == int(row["incidences"]) and
                        len(incidence) in (1, 2),
                        f"canonical owner/incidence count mismatch: {key}")
                for field in SIMPLE:
                    require(all(int(w[field]) == int(row[field]) for w in incidence),
                            f"canonical {field} mismatch: {key}")
                ordered = sorted(incidence, key=lambda w: int(w["incident_triangle"]))
                require(len({int(w["incident_triangle"]) for w in ordered}) ==
                        len(ordered), f"repeated triangle incidence: {key}")
                for index, wire in enumerate(ordered):
                    for field in INCIDENCE:
                        require(int(row[f"{field}{index}"]) == int(wire[field]),
                                f"canonical {field}{index} mismatch: {key}")
                lower = max(Decimal(w["event_lower"]) for w in incidence)
                upper = min(Decimal(w["event_upper"]) for w in incidence)
                require(lower <= upper and
                        Decimal(row["event_lower"]) == lower and
                        Decimal(row["event_upper"]) == upper,
                        f"canonical binary128 root intersection mismatch: {key}")
                actual += 1
        require(actual == declared, f"native/canonical owner count mismatch: {owner}")
        owner_counts.append(actual)
        events += actual
    require(not groups, f"incidence wires absent from canonical ledger: {len(groups)}")
    receipt = re.findall(
        rf"NATIVE_SOURCE_EDGE_CHUNK_PASS first={first} last={last} "
        rf"incidence_wires=(\d+) canonical_events=(\d+) mass_law_solved=0", log
    )
    require(len(receipt) == 1 and
            tuple(map(int, receipt[0])) == (wires, events) and
            wires == sum(count for _, count in ranks),
            f"native totals or wire count mismatch: {canonical}")
    print(f"PASS canonical [{first},{last}) wires={wires} "
          f"events={events} per_owner={owner_counts}")
    return first, last, total, seal, wires, events


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("chunks", nargs="+", type=Path)
    args = parser.parse_args()
    rows = sorted((audit_one(chunk.resolve()) for chunk in args.chunks),
                  key=lambda row: row[0])
    expected = rows[0][0]
    total = rows[0][2]
    seal = rows[0][3]
    for first, last, current_total, current_seal, _, _ in rows:
        require(first == expected and current_total == total and current_seal == seal,
                f"gap, overlap or source/chart/weight seal mismatch before {first}")
        expected = last
    print(f"PASS canonical contiguous=[{rows[0][0]},{expected}) "
          f"wires={sum(row[4] for row in rows)} "
          f"events={sum(row[5] for row in rows)}")


if __name__ == "__main__":
    main()
