#!/usr/bin/env python3
"""Materialize the already audited per-triangle boundary cuts without rerunning MPI.

This is a sparse geometry cache for later source/target overlap construction.
It is not an overlap-region or mass certificate. Native routing remains the
authority for target-node identity and original-triangle ownership.
"""

import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from decimal import Decimal, localcontext
from pathlib import Path

from audit_power_triangle_route_chunks import audit_one


HEADER = ("source_triangle", "local_edge", "source_edge", "orientation",
          "lower", "upper", "target_segment", "target_endpoint",
          "source_vertex", "star_count", "source_rank", "edge_owner")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def native_counts(log, first, last):
    sent = {}
    received = {}
    for rank, s, r in re.findall(
            rf"NATIVE_SOURCE_TRIANGLE_CUT_ROUTE_PASS rank=(\d+) "
            rf"edge_range=\[{first},{last}\) sent=(\d+) received=(\d+)", log):
        rank = int(rank)
        require(rank not in sent, "duplicate native route rank")
        sent[rank] = int(s)
        received[rank] = int(r)
    require(set(sent) == set(received) == set(range(4)), "complete native route counts")
    return sent, received


def canonical_items(path):
    with path.open(newline="") as source:
        reader = csv.DictReader(source)
        for row in reader:
            lower = Decimal(row["event_lower"])
            upper = Decimal(row["event_upper"])
            require(lower.is_finite() and upper.is_finite() and
                    Decimal(0) <= lower <= upper <= Decimal(1),
                    "bounded canonical binary128 interval")
            edge = int(row["event_edge"])
            target = int(row["target_segment"])
            endpoint = int(row["target_endpoint"])
            vertex = int(row["source_vertex"])
            count = int(row["incidences"])
            require(count in (1, 2) and endpoint in (-1, 0, 1),
                    "canonical incidence and target endpoint")
            incidence = tuple((int(row[f"source_rank{k}"]),
                               int(row[f"incident_triangle{k}"]),
                               int(row[f"local_edge{k}"]),
                               int(row[f"orientation{k}"])) for k in range(count))
            for rank, triangle, local_edge, orientation in incidence:
                require(rank in range(4) and triangle >= 0 and
                        local_edge in range(3) and orientation in (-1, 1),
                        "bounded original triangle incidence")
            yield (edge, lower, upper, target, endpoint, vertex,
                   int(row["edge_owner"]), incidence,
                   row["event_lower"], row["event_upper"])


def flush_group(group, by_rank, sent):
    if not group:
        return
    first = group[0]
    edge = first[0]
    owner = first[6]
    incidence = first[7]
    lower = max(item[1] for item in group)
    upper = min(item[2] for item in group)
    require(lower <= upper and owner in range(4), "one nonempty canonical cut")
    require(all(item[0] == edge and item[6] == owner and
                item[7] == incidence for item in group),
            "one canonical edge owner and incidence set per grouped cut")
    if len(group) > 1:
        require(all(item[4] >= 0 and item[5] == first[5] for item in group),
                "overlapping canonical intervals require same target-node vertex")
    members = sorted((item[3], item[4]) for item in group)
    require(len(members) == len(set(members)), "unique target-node star members")
    representative = members[0]
    lower_text = next(item[8] for item in group if item[1] == lower)
    upper_text = next(item[9] for item in group if item[2] == upper)
    for rank, triangle, local_edge, orientation in incidence:
        by_rank[rank].append((triangle, local_edge, edge, orientation,
                              lower_text, upper_text, representative[0],
                              representative[1], first[5], len(members),
                              rank, owner))
        sent[owner] += 1


def main(run):
    run = run.resolve()
    match = re.fullmatch(
        r"power_edge_triangle_route_chunk_(\d+)_(\d+)_(root60_)?v1", run.name)
    require(match is not None, "native route run directory name")
    first, last = map(int, match.group(1, 2))
    variant = "_root60_v1" if match.group(3) else "_v2"
    canonical = run.parent / f"power_edge_reconcile_{first}_{last}{variant}"
    output = run.parent / f"power_triangle_cut_index_{first}_{last}_v1"
    temporary = run.parent / f"power_triangle_cut_index_{first}_{last}_v1.tmp"
    require(not output.exists() and not temporary.exists(), "refuse overwrite cut index")
    audit_one(run)
    expected_sent, expected_received = native_counts(
        (run / "native.log").read_text(), first, last)
    by_rank = defaultdict(list)
    sent = defaultdict(int)
    input_sha = {}
    with localcontext() as context:
        context.prec = 100
        for owner in range(4):
            path = canonical / f"canonical_owner{owner}.csv"
            input_sha[path.name] = digest(path)
            items = sorted(canonical_items(path), key=lambda item:
                           (item[0], item[1], item[2], item[3]))
            group = []
            for item in items:
                if group and (item[0] != group[0][0] or
                              item[1] > min(member[2] for member in group)):
                    flush_group(group, by_rank, sent)
                    group = []
                if group:
                    require(item[4] >= 0 and all(member[4] >= 0 for member in group),
                            "non-node intervals may not overlap")
                group.append(item)
            flush_group(group, by_rank, sent)
    received = {rank: len(by_rank[rank]) for rank in range(4)}
    require({rank: sent[rank] for rank in range(4)} == expected_sent and
            received == expected_received,
            "materialized cut fanout matches all four native MPI route receipts")
    temporary.mkdir()
    output_sha = {}
    for rank in range(4):
        rows = by_rank[rank]
        rows.sort(key=lambda row: (row[0], row[1],
                                  Decimal(row[4]) if row[3] == 1 else Decimal(row[5]).copy_negate(),
                                  Decimal(row[5]) if row[3] == 1 else Decimal(row[4]).copy_negate(),
                                  row[6]))
        for previous, current in zip(rows, rows[1:]):
            if previous[:2] != current[:2]:
                continue
            require(previous[2] == current[2] and previous[3] == current[3],
                    "one oriented source edge per original triangle side")
            if current[3] == 1:
                require(Decimal(previous[5]) < Decimal(current[4]),
                        "strict forward cut order on original triangle side")
            else:
                require(Decimal(current[5]) < Decimal(previous[4]),
                        "strict reverse cut order on original triangle side")
        path = temporary / f"triangle_cut_rank{rank}.csv"
        with path.open("w", newline="") as destination:
            writer = csv.writer(destination)
            writer.writerow(HEADER)
            writer.writerows(rows)
        output_sha[path.name] = digest(path)
    receipt = {"range": [first, last], "native_route": str(run),
               "exporter_sha256": digest(Path(__file__)),
               "native_log_sha256": digest(run / "native.log"),
               "canonical_sha256": input_sha, "output_sha256": output_sha,
               "sent_per_owner": expected_sent,
               "received_per_triangle_rank": expected_received,
               "physical_overlap_mass_solved": False}
    (temporary / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    temporary.rename(output)
    print(f"TRIANGLE_CUT_INDEX_PASS [{first},{last}) "
          f"wires={sum(received.values())} output={output}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: export_power_triangle_cuts_v1.py NATIVE_ROUTE_RUN")
    main(Path(sys.argv[1]))
