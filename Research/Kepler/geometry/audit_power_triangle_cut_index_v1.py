#!/usr/bin/env python3
"""Independently check sparse per-triangle cut caches against sealed native routes."""

import argparse
import csv
import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path

from audit_power_triangle_route_chunks import audit_one as audit_native_route


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def audit_one(directory):
    directory = directory.resolve()
    match = re.fullmatch(r"power_triangle_cut_index_(\d+)_(\d+)_v1", directory.name)
    require(match is not None, "sparse triangle-cut index directory name")
    first, last = map(int, match.groups())
    receipt = json.loads((directory / "receipt.json").read_text())
    require(receipt["range"] == [first, last] and
            receipt["physical_overlap_mass_solved"] is False,
            "bounded geometry-only sparse index receipt")
    route = Path(receipt["native_route"])
    require(route.parent == directory.parent and
            route.name == f"power_edge_triangle_route_chunk_{first}_{last}" +
            ("_v1" if first < 500000 else "_root60_v1"),
            "matching native route input for sparse cache")
    audit_native_route(route)
    require(digest(route / "native.log") == receipt["native_log_sha256"],
            "native MPI route log unchanged since sparse export")
    canonical = directory.parent / (f"power_edge_reconcile_{first}_{last}" +
        ("_v2" if first < 500000 else "_root60_v1"))
    for name, expected in receipt["canonical_sha256"].items():
        require(re.fullmatch(r"canonical_owner[0-3]\.csv", name) and
                digest(canonical / name) == expected,
                "canonical owner ledger unchanged since sparse export")
    require(set(receipt["canonical_sha256"]) ==
            {f"canonical_owner{rank}.csv" for rank in range(4)},
            "all four canonical owner ledgers sealed")
    log = (route / "native.log").read_text()
    expected_sent = {}
    expected_received = {}
    for rank, sent, received in re.findall(
            rf"NATIVE_SOURCE_TRIANGLE_CUT_ROUTE_PASS rank=(\d+) "
            rf"edge_range=\[{first},{last}\) sent=(\d+) received=(\d+)", log):
        expected_sent[int(rank)] = int(sent)
        expected_received[int(rank)] = int(received)
    require(set(expected_sent) == set(expected_received) == set(range(4)),
            "complete native MPI route count receipt")
    require({int(k): v for k, v in receipt["sent_per_owner"].items()} ==
            expected_sent and
            {int(k): v for k, v in receipt["received_per_triangle_rank"].items()} ==
            expected_received,
            "sparse cache fanout receipt matches native MPI route")
    counts = {}
    for rank in range(4):
        name = f"triangle_cut_rank{rank}.csv"
        path = directory / name
        require(digest(path) == receipt["output_sha256"][name],
                "sparse triangle-cut CSV digest unchanged")
        count = 0
        previous = None
        with path.open(newline="") as source:
            reader = csv.DictReader(source)
            require(set(reader.fieldnames or ()) == {
                "source_triangle", "local_edge", "source_edge", "orientation",
                "lower", "upper", "target_segment", "target_endpoint",
                "source_vertex", "star_count", "source_rank", "edge_owner"},
                "complete sparse triangle-cut fields")
            for row in reader:
                triangle = int(row["source_triangle"])
                side = int(row["local_edge"])
                edge = int(row["source_edge"])
                orientation = int(row["orientation"])
                lower = Decimal(row["lower"])
                upper = Decimal(row["upper"])
                require(triangle >= 0 and side in range(3) and
                        first <= edge < last and orientation in (-1, 1) and
                        lower.is_finite() and upper.is_finite() and
                        Decimal(0) <= lower <= upper <= Decimal(1) and
                        int(row["source_rank"]) == rank and
                        int(row["edge_owner"]) in range(4) and
                        int(row["star_count"]) >= 1 and
                        int(row["target_segment"]) >= 0 and
                        int(row["target_endpoint"]) in (-1, 0, 1),
                        "bounded sealed oriented triangle cut row")
                current = (triangle, side, edge, orientation, lower, upper)
                if previous is not None:
                    require(previous[:2] <= current[:2],
                            "sorted sparse source-triangle side order")
                    if previous[:2] == current[:2]:
                        require(previous[2:4] == current[2:4] and
                                (previous[5] < lower if orientation == 1
                                 else upper < previous[4]),
                                "strict ordered cuts on same oriented triangle side")
                previous = current
                count += 1
        require(count == expected_received[rank],
                "sparse triangle-cut count agrees with native MPI receiver")
        counts[rank] = count
    require(set(receipt["output_sha256"]) ==
            {f"triangle_cut_rank{rank}.csv" for rank in range(4)},
            "all four sparse receiver CSV digests")
    print(f"PASS triangle cut index [{first},{last}) "
          f"wires={sum(counts.values())} received={list(counts.values())}")
    return first, last, sum(counts.values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=Path)
    rows = sorted((audit_one(path) for path in parser.parse_args().directories),
                  key=lambda row: row[0])
    expected = rows[0][0]
    for first, last, _ in rows:
        require(first == expected, "sparse triangle-cut index range gap or overlap")
        expected = last
    print(f"PASS triangle cut index contiguous=[{rows[0][0]},{expected}) "
          f"wires={sum(row[2] for row in rows)}")


if __name__ == "__main__":
    main()
