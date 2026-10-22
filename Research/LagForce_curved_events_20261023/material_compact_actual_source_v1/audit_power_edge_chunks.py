#!/usr/bin/env python3
"""Check sealed, deterministic source-edge event chunks before reconciliation."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def check_sha_file(path, relative_base=None):
    for line in path.read_text().splitlines():
        expected, filename = line.split("  ", 1)
        source = Path(filename)
        if not source.is_absolute():
            source = (relative_base or path.parent) / source
        require(source.is_file(), f"missing sealed file {source}")
        with source.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        require(digest == expected, f"SHA256 mismatch {source}")


def read_one_csv(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 1, f"expected one manifest row: {path}")
    return rows[0]


def check_source_snapshot(directory):
    path = directory / "source_snapshot.json"
    if not path.exists():
        check_sha_file(directory / "source.sha256", directory.parent)
        return
    with path.open() as stream:
        snapshot = json.load(stream)
    with (directory / "source.sha256").open("rb") as stream:
        actual = hashlib.file_digest(stream, "sha256").hexdigest()
    require(actual == snapshot["source_manifest_sha256"],
            f"source manifest changed: {directory}")
    originals = (directory / "source.sha256").read_text().splitlines()
    require(len(originals) == len(snapshot["files"]),
            f"source snapshot length mismatch: {directory}")
    for line, row in zip(originals, snapshot["files"]):
        expected, filename = line.split("  ", 1)
        original = Path(filename)
        if not original.is_absolute():
            original = directory.parent / original
        retained = directory / "source_snapshot" / row["retained"]
        require(str(original) == row["source"] and
                expected == row["sha256"] and retained.is_file(),
                f"source snapshot identity mismatch: {directory}")
        with retained.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        require(actual == expected, f"source snapshot changed: {retained}")


def audit_chunk(directory):
    directory = directory.resolve()
    require((directory / "native_exit_code.txt").read_text().strip() == "0",
            f"native run did not exit 0: {directory}")
    check_sha_file(directory / "source_edge_wires.sha256")
    check_source_snapshot(directory)
    check_sha_file(directory / "binary.sha256")
    if (directory / "h200_cuda_binary.sha256").exists():
        check_sha_file(directory / "h200_cuda_binary.sha256")
    require("FAILED" not in (directory / "final_source_check.txt").read_text(),
            f"source changed during run: {directory}")
    manifests = [read_one_csv(directory / f"source_edge_manifest_rank{rank}.csv")
                 for rank in range(4)]
    log = (directory / "native.log").read_text()
    require("WEIGHTED_FAIL" not in log and "SOURCE_EDGE_EVENT_PRODUCTION_FAIL" not in log,
            f"native failure marker: {directory}")
    common = ("edge_range_first", "edge_range_end", "global_edges",
              "source_registry_hash", "source_epoch", "source_steps",
              "source_time_bits", "mass_law_hash", "target_chart_hash",
              "weight_hash")
    seal = tuple(manifests[0][field] for field in common)
    first = int(manifests[0]["edge_range_first"])
    last = int(manifests[0]["edge_range_end"])
    total = int(manifests[0]["global_edges"])
    require(0 <= first < last <= total, f"invalid edge range: {directory}")
    require(directory.name in (f"power_edge_chunk_{first}_{last}",
                               f"power_edge_chunk_{first}_{last}_h200sm90",
                               f"power_edge_chunk_{first}_{last}_h200sm90_v2",
                               f"power_edge_chunk_{first}_{last}_h200sm90_v2b",
                               f"power_edge_chunk_{first}_{last}_h200sm90_root60",
                               f"power_edge_chunk_{first}_{last}_lageunha_root60"),
            f"directory/range mismatch: {directory}")
    output = []
    for rank, manifest in enumerate(manifests):
        pattern = (rf"ACTUAL_STEP5_PARALLEL_SOURCE_EDGE_PREFIX rank={rank} "
                   rf"edge_range=\[{first},{last}\) unique=(\d+) candidates=(\d+) "
                   rf"root_calls=(\d+) root_boxes=(\d+) wires=(\d+) "
                   rf"complete_global=(\d+)")
        match = re.findall(pattern, log)
        require(len(match) == 1 and
                match[0] == tuple(manifest[field] for field in
                    ("unique_edges", "candidate_segments", "root_calls",
                     "root_boxes", "wire_rows", "complete_global")),
                f"rank {rank} log/manifest count mismatch: {directory}")
        require(int(manifest["rank"]) == rank and
                tuple(manifest[field] for field in common) == seal,
                f"rank {rank} manifest identity mismatch: {directory}")
        require(int(manifest["complete_global"]) == int(first == 0 and last == total),
                f"rank {rank} global completion mismatch: {directory}")
        wire_path = directory / f"source_edge_wires_rank{rank}.csv"
        previous_edge = first
        count = 0
        with wire_path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                edge = int(row["event_edge"])
                require(previous_edge <= edge < last,
                        f"rank {rank} nonmonotone or out-of-range edge: {edge}")
                previous_edge = edge
                require(int(row["source_rank"]) == rank and
                        row["source_hash"] == manifest["source_shard_hash"] and
                        row["source_registry_hash"] == manifest["source_registry_hash"] and
                        row["source_epoch"] == manifest["source_epoch"] and
                        row["target_chart_hash"] == manifest["target_chart_hash"] and
                        row["weight_hash"] == manifest["weight_hash"] and
                        int(row["edge_range_first"]) == first and
                        int(row["edge_range_end"]) == last,
                        f"rank {rank} wire/manifest identity mismatch at edge {edge}")
                lower = Decimal(row["event_lower"])
                upper = Decimal(row["event_upper"])
                require(Decimal(0) <= lower <= upper <= Decimal(1),
                        f"rank {rank} invalid event interval at edge {edge}")
                require(0 <= int(row["edge_owner"]) < 4,
                        f"rank {rank} invalid canonical owner at edge {edge}")
                count += 1
        require(count == int(manifest["wire_rows"]),
                f"rank {rank} wire row count mismatch: {directory}")
        output.append((int(manifest["unique_edges"]), count))
    shards = tuple(row["source_shard_hash"] for row in manifests)
    return first, last, total, (seal[3:], shards), output


def audit_presence(directory, total, shards):
    directory = directory.resolve()
    require((directory / "native_exit_code.txt").read_text().strip() == "0",
            f"source presence run failed: {directory}")
    check_sha_file(directory / "presence.sha256")
    check_sha_file(directory / "binary.sha256")
    check_source_snapshot(directory)
    log = (directory / "native.log").read_text()
    present = []
    for rank in range(4):
        pattern = (rf"SOURCE_EDGE_PRESENCE_PASS rank={rank} "
                   rf"global_edges={total} local_unique=(\d+) "
                   rf"triangle_incidences=(\d+) global_covered={total} "
                   rf"source_epoch=(\d+) source_hash={shards[rank]}")
        require(len(re.findall(pattern, log)) == 1,
                f"source presence rank {rank} log/seal mismatch")
        raw = (directory / f"presence_rank{rank}.bin").read_bytes()
        require(len(raw) == (total + 7) // 8,
                f"source presence bitmap size mismatch rank {rank}")
        if total % 8:
            require(raw[-1] >> (total % 8) == 0,
                    f"nonzero unused source edge bits rank {rank}")
        present.append(raw)
    for byte in range((total + 7) // 8):
        covered = present[0][byte] | present[1][byte] | present[2][byte] | present[3][byte]
        expected = (1 << (total % 8)) - 1 if byte == total // 8 and total % 8 else 255
        require(covered == expected,
                f"missing original global source edge near {8 * byte}")
    exterior_path = directory / "external_edge.bin"
    exterior = exterior_path.read_bytes() if exterior_path.exists() else None
    if exterior is not None:
        require(len(exterior) == total and set(exterior) <= {0, 1},
                f"invalid accepted source exterior edge flags: {directory}")
    return present, exterior


def bit_count_range(raw, first, last):
    return sum(bool(raw[edge // 8] & (1 << (edge % 8)))
               for edge in range(first, last))


def audit_event_pairs(directory, exterior, first, last):
    groups = defaultdict(list)
    for rank in range(4):
        with (directory / f"source_edge_wires_rank{rank}.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                key = tuple(int(row[field]) for field in
                            ("event_edge", "event_cell_low", "event_cell_high",
                             "event_ordinal", "target_segment"))
                groups[key].append(row)
    for key, rows in groups.items():
        edge, low, high, _, _ = key
        require(first <= edge < last and low < high,
                f"invalid canonical source event key {key}")
        needed = 1 if exterior[edge] else 2
        require(len(rows) == needed,
                f"incomplete source triangle incidences for event {key}")
        require(len({row["edge_owner"] for row in rows}) == 1 and
                len({row["incident_triangle"] for row in rows}) == needed and
                len({row["event_derivative_sign"] for row in rows}) == 1 and
                len({row["target_endpoint"] for row in rows}) == 1 and
                len({row["source_vertex"] for row in rows}) == 1,
                f"inconsistent source event incidence identity {key}")
        if needed == 2:
            require(sum(int(row["orientation"]) for row in rows) == 0 and
                    max(Decimal(row["event_lower"]) for row in rows) <=
                    min(Decimal(row["event_upper"]) for row in rows),
                    f"opposing source edge orientation/root intersection failed {key}")
    return len(groups)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=Path)
    parser.add_argument("--require-full", action="store_true")
    parser.add_argument("--presence-directory", type=Path)
    args = parser.parse_args()
    chunks = sorted((audit_chunk(path) for path in args.directories),
                    key=lambda item: item[0])
    reference_total, reference_seal = chunks[0][2:4]
    present = None
    exterior = None
    if args.presence_directory:
        present, exterior = audit_presence(args.presence_directory, reference_total,
                                           reference_seal[1])
    expected = chunks[0][0]
    for first, last, total, seal, ranks in chunks:
        require(first == expected, f"gap/overlap before edge {first}")
        require(total == reference_total and seal == reference_seal,
                f"cross-chunk source/target seal mismatch at edge {first}")
        if present is not None:
            for rank, (unique, _) in enumerate(ranks):
                require(unique == bit_count_range(present[rank], first, last),
                        f"rank {rank} processed-edge count/presence mismatch in [{first},{last})")
        if exterior is not None:
            events = audit_event_pairs(
                next(path.resolve() for path in args.directories
                     if path.name.startswith(f"power_edge_chunk_{first}_{last}")),
                exterior, first, last)
            print(f"PASS canonical event incidence candidates [{first},{last})={events}")
        print(f"PASS [{first},{last}) global={total} "
              f"per_rank={ranks} wires={sum(row[1] for row in ranks)}")
        expected = last
    if args.require_full:
        require(chunks[0][0] == 0 and expected == reference_total,
                f"incomplete global edge tiling: [0,{expected}) of {reference_total}")
    print(f"PASS contiguous chunks={len(chunks)} covered=[{chunks[0][0]},{expected})")


if __name__ == "__main__":
    main()
