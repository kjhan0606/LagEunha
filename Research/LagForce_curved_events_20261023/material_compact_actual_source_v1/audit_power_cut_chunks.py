#!/usr/bin/env python3
"""Audit saved native source-edge cut-order gates and their provenance."""

import argparse
import re
from pathlib import Path

from audit_power_edge_chunks import (
    check_sha_file,
    check_source_snapshot,
    require,
)


def audit_one(directory):
    directory = directory.resolve()
    match = re.fullmatch(r"power_edge_cut_chunk_(\d+)_(\d+)_(root60_)?v1", directory.name)
    require(match is not None, f"invalid cut-order run name: {directory}")
    first, last = map(int, match.groups()[:2])
    variant = "_root60_v1" if match.group(3) else "_v2"
    require(first < last and last - first <= 50000,
            f"invalid bounded original-edge range: {directory}")
    require((directory / "native_exit_code.txt").read_text().strip() == "0",
            f"native cut ordering did not exit zero: {directory}")
    check_source_snapshot(directory)
    check_sha_file(directory / "binary.sha256")
    require("FAILED" not in (directory / "final_source_check.txt").read_text(),
            f"cut input changed during native run: {directory}")
    log = (directory / "native.log").read_text()
    require("FAIL" not in log, f"native cut failure marker: {directory}")
    require(len(re.findall(r"ACTUAL_STEP5_WEIGHTED_DEGREE_PASS covered=1 "
                           r"minimum_J=([0-9.eE+-]+)", log)) == 4,
            f"missing complete positive chart/degree proof: {directory}")
    require(len(re.findall(
        r"ACTUAL_STEP5_TARGET_BOUNDARY_INDEX_PASS segments=642787 "
        r"physical=212217 internal=430570 nodes=262143 "
        r"trace_ties=10224 mass_law_solved=0", log)) == 4,
        f"missing complete target segment index: {directory}")
    events = cuts = grouped = 0
    for rank in range(4):
        pattern = (
            rf"NATIVE_SOURCE_EDGE_CUT_ORDER_PASS rank={rank} "
            rf"edge_range=\[{first},{last}\) events=(\d+) cuts=(\d+) "
            rf"grouped_target_node_events=(\d+) mass_law_solved=0"
        )
        found = re.findall(pattern, log)
        require(len(found) == 1, f"missing unique native cut rank receipt: {directory}")
        rank_events, rank_cuts, rank_grouped = map(int, found[0])
        require(rank_cuts + rank_grouped == rank_events,
                f"canonical event omitted or duplicated by cut order: {directory}")
        canonical = (directory.parent /
                     f"power_edge_reconcile_{first}_{last}{variant}" /
                     f"canonical_owner{rank}.csv")
        with canonical.open() as stream:
            rows = sum(1 for _ in stream) - 1
        require(rows == rank_events,
                f"saved canonical owner event count mismatch: {directory}, rank {rank}")
        events += rank_events
        cuts += rank_cuts
        grouped += rank_grouped
    print(f"PASS cut order [{first},{last}) events={events} cuts={cuts} "
          f"grouped_target_node_events={grouped}")
    return first, last, events, cuts, grouped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=Path)
    args = parser.parse_args()
    rows = sorted((audit_one(directory) for directory in args.directories),
                  key=lambda row: row[0])
    expected = rows[0][0]
    for first, last, _, _, _ in rows:
        require(first == expected, f"cut-order range gap or overlap before {first}")
        expected = last
    print(f"PASS cut order contiguous=[{rows[0][0]},{expected}) "
          f"events={sum(row[2] for row in rows)} "
          f"cuts={sum(row[3] for row in rows)} "
          f"grouped_target_node_events={sum(row[4] for row in rows)}")


if __name__ == "__main__":
    main()
