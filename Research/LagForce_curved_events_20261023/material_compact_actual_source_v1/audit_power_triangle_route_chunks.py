#!/usr/bin/env python3
"""Audit native original-triangle cut routing receipts for a sealed edge range."""

import argparse
import re
from pathlib import Path

from audit_power_edge_chunks import check_sha_file, check_source_snapshot, require


def audit_one(directory):
    directory = directory.resolve()
    match = re.fullmatch(r"power_edge_triangle_route_chunk_(\d+)_(\d+)_(root60_)?v1", directory.name)
    require(match is not None, f"invalid triangle-route directory: {directory}")
    first, last = map(int, match.group(1, 2))
    canonical_variant = "_root60_v1" if match.group(3) else "_v2"
    require(first < last and last - first <= 50000, f"invalid edge range: {directory}")
    require((directory / "native_exit_code.txt").read_text().strip() == "0",
            f"native route failed: {directory}")
    check_source_snapshot(directory)
    check_sha_file(directory / "binary.sha256")
    require("FAILED" not in (directory / "final_source_check.txt").read_text(),
            f"route input changed during native run: {directory}")
    log = (directory / "native.log").read_text()
    require("FAIL" not in log, f"native route failure marker: {directory}")
    require(len(re.findall(r"ACTUAL_STEP5_TARGET_BOUNDARY_INDEX_PASS segments=642787 "
                           r"physical=212217 internal=430570 nodes=262143 "
                           r"trace_ties=10224 mass_law_solved=0", log)) == 4,
            f"missing target chart boundary gate: {directory}")
    sent = received = events = cuts = grouped = 0
    for rank in range(4):
        cut_match = re.findall(
            rf"NATIVE_SOURCE_EDGE_CUT_ORDER_PASS rank={rank} "
            rf"edge_range=\[{first},{last}\) events=(\d+) cuts=(\d+) "
            rf"grouped_target_node_events=(\d+) mass_law_solved=0", log)
        route_match = re.findall(
            rf"NATIVE_SOURCE_TRIANGLE_CUT_ROUTE_PASS rank={rank} "
            rf"edge_range=\[{first},{last}\) sent=(\d+) received=(\d+) "
            rf"owned_triangles=(\d+) triangle_edge_ranges=(\d+) mass_law_solved=0", log)
        require(len(cut_match) == len(route_match) == 1,
                f"missing unique native cut/route receipt for rank {rank}: {directory}")
        e, c, g = map(int, cut_match[0])
        s, r, t, er = map(int, route_match[0])
        require(e == c + g and er == 3 * t and s >= 0 and r >= 0,
                f"rank {rank} cut/route accounting failed: {directory}")
        canonical = directory.parent / f"power_edge_reconcile_{first}_{last}{canonical_variant}" / f"canonical_owner{rank}.csv"
        with canonical.open() as stream:
            require(sum(1 for _ in stream) - 1 == e,
                    f"rank {rank} canonical input count mismatch: {directory}")
        sent += s
        received += r
        events += e
        cuts += c
        grouped += g
    require(sent == received, f"MPI route lost or duplicated cut wires: {directory}")
    print(f"PASS triangle route [{first},{last}) canonical_events={events} "
          f"ordered_cuts={cuts} grouped_target_nodes={grouped} "
          f"sent={sent} received={received}")
    return first, last, events, cuts, grouped, sent, received


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=Path)
    rows = sorted((audit_one(directory) for directory in parser.parse_args().directories),
                  key=lambda row: row[0])
    expected = rows[0][0]
    for first, last, *_ in rows:
        require(first == expected, f"triangle-route range gap or overlap before {first}")
        expected = last
    print(f"PASS triangle route contiguous=[{rows[0][0]},{expected}) "
          f"canonical_events={sum(row[2] for row in rows)} "
          f"ordered_cuts={sum(row[3] for row in rows)} "
          f"grouped_target_nodes={sum(row[4] for row in rows)} "
          f"sent={sum(row[5] for row in rows)} "
          f"received={sum(row[6] for row in rows)}")


if __name__ == "__main__":
    main()
