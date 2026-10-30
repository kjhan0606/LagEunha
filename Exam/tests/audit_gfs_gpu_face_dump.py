#!/usr/bin/env python3
"""Pair per-rank GPU face dumps by global particle IDs and audit cancellation."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median


FIELDS = (
    "rank", "call", "id_i", "id_j", "wall", "pstar", "udotds", "work",
    "dSx", "dSy", "x1", "y1", "x2", "y2",
)
STENCIL_FIELDS = ("kp_idx", "km_idx", "id_kp", "id_km")


def read_rows(paths):
    rows = []
    stencil_schemas = set()
    for path in paths:
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            missing = set(FIELDS) - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path}: missing columns {sorted(missing)}")
            present_stencil = set(STENCIL_FIELDS) & set(reader.fieldnames or ())
            if present_stencil and present_stencil != set(STENCIL_FIELDS):
                raise ValueError(f"{path}: incomplete face-stencil metadata")
            stencil_schemas.add(bool(present_stencil))
            for line, raw in enumerate(reader, start=2):
                try:
                    row = {key: int(raw[key]) for key in ("rank", "call", "id_i", "id_j", "wall")}
                    row.update({key: float(raw[key]) for key in FIELDS[5:]})
                    if present_stencil:
                        row.update({key: int(raw[key]) for key in STENCIL_FIELDS})
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"{path}:{line}: invalid row") from exc
                if not all(math.isfinite(row[key]) for key in FIELDS[5:]):
                    raise ValueError(f"{path}:{line}: non-finite measurement")
                rows.append(row)
    if not rows:
        raise ValueError("no face rows found")
    calls = {row["call"] for row in rows}
    if len(calls) != 1:
        raise ValueError(f"mixed force-call numbers: {sorted(calls)}")
    if len(stencil_schemas) > 1:
        raise ValueError("mixed face-dump schemas: some ranks lack stencil metadata")
    return rows


def summarize_group(pairs):
    work_sum = sum(pair["work_sum"] for pair in pairs)
    work_abs = sum(pair["work_abs"] for pair in pairs)
    return {
        "pairs": len(pairs),
        "work_residual_signed_sum": work_sum,
        "work_residual_abs_sum": sum(abs(pair["work_sum"]) for pair in pairs),
        "work_abs_contribution_sum": work_abs,
        "aggregate_cancellation_fraction": abs(work_sum) / max(work_abs, 1e-300),
        "median_pstar_rel_difference": median(pair["pstar_rel"] for pair in pairs),
        "max_pstar_rel_difference": max(pair["pstar_rel"] for pair in pairs),
        "median_normal_rel_mismatch": median(pair["normal_rel"] for pair in pairs),
        "max_normal_rel_mismatch": max(pair["normal_rel"] for pair in pairs),
        "median_udotds_antisymmetry_error": median(pair["udotds_rel"] for pair in pairs),
        "max_udotds_antisymmetry_error": max(pair["udotds_rel"] for pair in pairs),
    }


def audit(rows, top):
    groups = defaultdict(list)
    unpaired_boundaries = 0
    for row in rows:
        if row["wall"] != 0 or row["id_j"] < 0:
            unpaired_boundaries += 1
            continue
        groups[tuple(sorted((row["id_i"], row["id_j"])))].append(row)

    paired = {"same_rank": [], "rank_boundary": []}
    malformed = []
    for (id_a, id_b), face_rows in groups.items():
        if len(face_rows) != 2 or {row["id_i"] for row in face_rows} != {id_a, id_b}:
            malformed.append({"ids": [id_a, id_b], "rows": len(face_rows),
                              "ranks": sorted(row["rank"] for row in face_rows)})
            continue
        first, second = face_rows
        norm_a = math.hypot(first["dSx"], first["dSy"])
        norm_b = math.hypot(second["dSx"], second["dSy"])
        pair = {
            "ids": [id_a, id_b],
            "ranks": sorted((first["rank"], second["rank"])),
            "work_sum": first["work"] + second["work"],
            "work_abs": abs(first["work"]) + abs(second["work"]),
            "pstar_rel": abs(first["pstar"] - second["pstar"]) /
                         max(abs(first["pstar"]), abs(second["pstar"]), 1e-300),
            "normal_rel": math.hypot(first["dSx"] + second["dSx"],
                                     first["dSy"] + second["dSy"]) /
                          max(0.5 * (norm_a + norm_b), 1e-300),
            "udotds_rel": abs(first["udotds"] + second["udotds"]) /
                          max(abs(first["udotds"]), abs(second["udotds"]), 1e-300),
        }
        relation = "same_rank" if first["rank"] == second["rank"] else "rank_boundary"
        paired[relation].append(pair)

    stencil_available = all(field in rows[0] for field in STENCIL_FIELDS)
    return {
        "force_call": rows[0]["call"],
        "ranks": sorted({row["rank"] for row in rows}),
        "face_rows": len(rows),
        "stencil_metadata_available": stencil_available,
        "faces_with_missing_stencil_neighbor": (
            sum(row["kp_idx"] < 0 or row["km_idx"] < 0 for row in rows)
            if stencil_available else None
        ),
        "unique_interior_pairs": len(groups),
        "unpaired_wall_rows": unpaired_boundaries,
        "malformed_interior_pairs": malformed,
        "by_pair_location": {
            name: summarize_group(items) if items else {"pairs": 0}
            for name, items in paired.items()
        },
        "largest_pair_work_residuals": {
            name: sorted(items, key=lambda pair: abs(pair["work_sum"]), reverse=True)[:top]
            for name, items in paired.items()
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path, help="face_dump_r*_call*.tsv files")
    parser.add_argument("--top", type=int, default=10, help="largest pair residuals to include")
    args = parser.parse_args()
    if args.top < 0:
        parser.error("--top must be nonnegative")
    print(json.dumps(audit(read_rows(args.files), args.top), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
