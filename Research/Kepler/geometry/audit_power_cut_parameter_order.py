#!/usr/bin/env python3
"""Preflight overlapping canonical source-edge root intervals.

This does not prove target-node identity. The native chart cut-order gate
must still check each overlapping target-node star.
"""

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
from pathlib import Path


def scan(directory):
    edges = defaultdict(list)
    events = 0
    for owner in range(4):
        with (directory / f"canonical_owner{owner}.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                edge = int(row["event_edge"])
                lo, hi = Decimal(row["event_lower"]), Decimal(row["event_upper"])
                if not (0 <= lo <= hi <= 1):
                    raise ValueError(f"invalid root interval in {directory}: {edge}")
                edges[edge].append((lo, hi, int(row["target_endpoint"]),
                                    int(row["source_vertex"]), int(row["target_segment"])))
                events += 1
    overlaps = nonnode = 0
    first_bad = None
    for edge, cuts in edges.items():
        cuts.sort()
        active = []
        for cut in cuts:
            active = [prior for prior in active if prior[1] >= cut[0]]
            for prior in active:
                overlaps += 1
                if prior[2] < 0 or cut[2] < 0:
                    nonnode += 1
                    if first_bad is None:
                        first_bad = (edge, prior, cut)
            active.append(cut)
    print(f"PREFLIGHT {directory.name} edges={len(edges)} events={events} "
          f"overlapping_pairs={overlaps} nonnode_overlapping_pairs={nonnode}")
    if nonnode:
        raise ValueError(f"overlapping intervals at nonnode cuts: {first_bad}")
    return events, overlaps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=Path)
    args = parser.parse_args()
    total_events = total_overlaps = 0
    for directory in args.directories:
        events, overlaps = scan(directory)
        total_events += events
        total_overlaps += overlaps
    print(f"PREFLIGHT_TOTAL events={total_events} "
          f"overlapping_pairs={total_overlaps} nonnode_overlapping_pairs=0")


if __name__ == "__main__":
    main()
