#!/usr/bin/env python3
"""Compare duplicated owner/halo particle inputs captured on GPU MPI ranks."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path


NAME = re.compile(r"face_input_r(\d+)_call(\d+)\.tsv$")
CAPTURE_VERSION = 3
FIELDS = (
    "x", "y", "vx", "vy", "mass", "V", "rho", "P", "cs", "U", "W", "Wold",
    "Bxx", "Bxy", "Byx", "Byy", "dPdx", "dPdy", "dRhodx", "dRhody",
    "tauxx", "tauxy", "tauyy", "alpha_cd", "divv", "lap_vx", "lap_vy", "U64", "V64",
)
PARAM_FIELDS = (
    "av_mode", "use_muscl", "OoA", "Courant", "Gamma", "alphavis", "betavis",
    "etavis", "epsvis", "nu_phys", "prandtl", "cd_amax", "blend_theta", "dtold",
    "hyperv_alpha", "hyperv_force_cap", "phase1",
    "laguerre_rot", "face_charge", "geom_fv", "w2lin", "thermal64",
    "var_geom_cfl", "var_geom_audit", "GFS_VARPRESS_W0",
)


def read_files(paths):
    copies = defaultdict(lambda: defaultdict(list))
    calls = set()
    ranks = set()
    params_by_rank = {}
    total_rows = 0
    for path in paths:
        match = NAME.search(path.name)
        if not match:
            raise ValueError(f"unexpected input filename: {path}")
        rank, call = map(int, match.groups())
        ranks.add(rank)
        calls.add(call)
        meta_path = path.with_name(f"face_input_meta_r{rank}_call{call}.json")
        if not meta_path.is_file():
            raise ValueError(f"missing metadata: {meta_path}")
        meta = json.loads(meta_path.read_text())
        if meta.get("rank") != rank or meta.get("call") != call or meta.get("all_inputs") != 1:
            raise ValueError(f"{meta_path}: capture is not a matching all-input dump")
        if meta.get("capture_version") != CAPTURE_VERSION:
            raise ValueError(f"{meta_path}: expected capture version {CAPTURE_VERSION}")
        if rank in params_by_rank:
            raise ValueError(f"multiple input files supplied for rank {rank}")
        missing_params = set(PARAM_FIELDS) - set(meta)
        if missing_params:
            raise ValueError(f"{meta_path}: missing force parameters {sorted(missing_params)}")
        params = {field: meta[field] for field in PARAM_FIELDS}
        if any(not isinstance(value, (int, float)) or not math.isfinite(float(value))
               for value in params.values()):
            raise ValueError(f"{meta_path}: non-finite or non-numeric force parameter")
        params_by_rank[rank] = params
        file_rows = 0
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            missing = {"index", "id", *FIELDS} - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path}: missing columns {sorted(missing)}")
            for line, raw in enumerate(reader, start=2):
                try:
                    particle_id = int(raw["id"])
                    values = {field: float(raw[field]) for field in FIELDS}
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"{path}:{line}: invalid row") from exc
                if not all(math.isfinite(value) for value in values.values()):
                    raise ValueError(f"{path}:{line}: non-finite state")
                copies[particle_id][rank].append(values)
                total_rows += 1
                file_rows += 1
        if meta.get("rows") != file_rows or file_rows == 0:
            raise ValueError(f"{path}: row count disagrees with metadata or is empty")
    if not copies:
        raise ValueError("no particle input rows found")
    if len(ranks) < 2:
        raise ValueError("at least two MPI-rank input files are required")
    if len(calls) != 1:
        raise ValueError(f"mixed force-call numbers: {sorted(calls)}")
    return copies, sorted(ranks), calls.pop(), total_rows, params_by_rank


def audit(copies, ranks, call, total_rows, params_by_rank, top):
    duplicate_ids = {
        particle_id: {rank: rows for rank, rows in by_rank.items() if len(rows) != 1}
        for particle_id, by_rank in copies.items()
        if any(len(rows) != 1 for rows in by_rank.values())
    }
    shared = {
        particle_id: by_rank
        for particle_id, by_rank in copies.items()
        if len(by_rank) > 1 and particle_id not in duplicate_ids
    }
    differences = {}
    for field in FIELDS:
        mismatches = []
        max_abs = 0.0
        max_rel = 0.0
        for particle_id, by_rank in shared.items():
            values = {rank: rows[0][field] for rank, rows in by_rank.items()}
            span = max(values.values()) - min(values.values())
            if span == 0:
                continue
            scale = max(abs(value) for value in values.values())
            rel = span / scale if scale else math.inf
            max_abs = max(max_abs, span)
            max_rel = max(max_rel, rel)
            mismatches.append({"id": particle_id, "ranks": sorted(values),
                               "values": values, "abs_span": span,
                               "relative_span": rel})
        mismatches.sort(key=lambda item: item["abs_span"], reverse=True)
        differences[field] = {
            "mismatching_ids": len(mismatches),
            "max_abs_span": max_abs,
            "max_relative_span": max_rel,
            "largest_mismatches": mismatches[:top],
        }
    parameter_mismatches = {
        field: {str(rank): params_by_rank[rank][field] for rank in ranks}
        for field in PARAM_FIELDS
        if len({params_by_rank[rank][field] for rank in ranks}) > 1
    }
    return {
        "capture_version": CAPTURE_VERSION,
        "force_call": call,
        "ranks": ranks,
        "input_rows": total_rows,
        "unique_global_ids": len(copies),
        "ids_present_on_multiple_ranks": len(shared),
        "ids_duplicated_within_a_rank": len(duplicate_ids),
        "rank_parameter_mismatches": parameter_mismatches,
        "field_copy_differences": differences,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path,
                        help="face_input_r*_call*.tsv files and matching metadata")
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--expected-ranks", type=int,
                        help="require exactly ranks 0 through N-1")
    args = parser.parse_args()
    if args.top < 0:
        parser.error("--top must be nonnegative")
    if args.expected_ranks is not None and args.expected_ranks < 1:
        parser.error("--expected-ranks must be positive")
    data = read_files(args.files)
    if args.expected_ranks is not None:
        expected = list(range(args.expected_ranks))
        if data[1] != expected:
            parser.error(f"expected ranks {expected}, found {data[1]}")
    print(json.dumps(audit(*data, args.top), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
