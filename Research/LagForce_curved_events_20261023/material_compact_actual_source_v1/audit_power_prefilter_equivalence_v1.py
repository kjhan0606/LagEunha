#!/usr/bin/env python3
"""Compare a certified plane-box skip with the unchanged root-search control."""

import argparse
import csv
import re
from pathlib import Path

from audit_power_edge_chunks import check_sha_file, check_source_snapshot, read_one_csv, require


METRICS = ("unique_edges", "candidate_segments", "root_calls", "root_boxes", "wire_rows")
SEAL = ("global_edges", "source_registry_hash", "source_shard_hash",
        "source_epoch", "source_steps", "source_time_bits", "mass_law_hash",
        "target_chart_hash", "weight_hash")


def audit_run(directory):
    directory = directory.resolve()
    require((directory / "native_exit_code.txt").read_text().strip() == "0",
            "native source-edge run exited zero")
    check_source_snapshot(directory)
    for name in ("source_edge_wires.sha256", "binary.sha256", "h200_cuda_binary.sha256"):
        check_sha_file(directory / name)
    require("FAILED" not in (directory / "final_source_check.txt").read_text(),
            "native run preserved exact inputs")
    manifests = [read_one_csv(directory / f"source_edge_manifest_rank{rank}.csv")
                 for rank in range(4)]
    first = int(manifests[0]["edge_range_first"])
    last = int(manifests[0]["edge_range_end"])
    require(0 <= first < last <= int(manifests[0]["global_edges"]),
            "bounded original-edge control range")
    log = (directory / "native.log").read_text()
    for rank, manifest in enumerate(manifests):
        require(int(manifest["rank"]) == rank and
                int(manifest["edge_range_first"]) == first and
                int(manifest["edge_range_end"]) == last,
                "four complete native source-edge rank manifests")
        pattern = (rf"ACTUAL_STEP5_PARALLEL_SOURCE_EDGE_PREFIX rank={rank} "
                   rf"edge_range=\[{first},{last}\) unique=(\d+) "
                   rf"candidates=(\d+) root_calls=(\d+) root_boxes=(\d+) "
                   rf"wires=(\d+) complete_global=(\d+)")
        found = re.findall(pattern, log)
        require(len(found) == 1 and
                found[0] == tuple(manifest[key] for key in METRICS) +
                (manifest["complete_global"],),
                "native source-edge counters agree with rank manifest")
    return first, last, manifests


def rows(directory, rank):
    with (directory / f"source_edge_wires_rank{rank}.csv").open(newline="") as source:
        yield from csv.DictReader(source)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("prefilter", type=Path)
    args = parser.parse_args()
    base_first, base_last, baseline = audit_run(args.baseline)
    trial_first, trial_last, trial = audit_run(args.prefilter)
    require((base_first, base_last) == (trial_first, trial_last),
            "same exact original-edge range in control and prefilter")
    require(args.baseline.name ==
            f"power_edge_chunk_{base_first}_{base_last}_h200sm90_root60" and
            args.prefilter.name ==
            f"power_edge_chunk_{base_first}_{base_last}_h200sm90_root60_prefilter_v1",
            "distinct baseline and certified-filter output directories")
    totals = {key: [0, 0] for key in METRICS}
    for rank in range(4):
        left, right = baseline[rank], trial[rank]
        require(tuple(left[key] for key in SEAL) ==
                tuple(right[key] for key in SEAL),
                "same accepted source, target chart, weight and shard seal")
        for key in METRICS:
            totals[key][0] += int(left[key])
            totals[key][1] += int(right[key])
        require(left["unique_edges"] == right["unique_edges"] and
                left["candidate_segments"] == right["candidate_segments"] and
                left["wire_rows"] == right["wire_rows"] and
                int(right["root_calls"]) <= int(left["root_calls"]) and
                int(right["root_boxes"]) <= int(left["root_boxes"]),
                "filter retains candidates and wires while only skipping certified roots")
        original = list(rows(args.baseline, rank))
        filtered = list(rows(args.prefilter, rank))
        require(original == filtered,
                f"rank {rank} original-edge wire CSV rows match exactly")
        require(len(original) == int(left["wire_rows"]),
                "every same-edge event wire compared")
    print(f"PASS plane-box prefilter exact wires range=[{base_first},{base_last}) "
          f"rows={totals['wire_rows'][0]} "
          f"candidates={totals['candidate_segments'][0]} "
          f"root_calls={totals['root_calls'][0]}->{totals['root_calls'][1]} "
          f"root_boxes={totals['root_boxes'][0]}->{totals['root_boxes'][1]}")


if __name__ == "__main__":
    main()
