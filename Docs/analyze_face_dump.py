#!/usr/bin/env python3
"""Pair directed GPU face rows by global particle IDs for one force call.

Usage: python3 Docs/analyze_face_dump.py face_dump_r{0,1,2,3}_callN.tsv
Pass exactly one call at a time; no directory scan is performed.
"""

import csv
import math
import sys
from collections import defaultdict


def main(paths):
    if not paths:
        raise SystemExit("provide the four rank TSV paths for one force call")
    pairs = defaultdict(list)
    calls = set()
    skipped_walls = 0
    for path in paths:
        with open(path, newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                calls.add(int(row["call"]))
                if int(row["wall"]):
                    skipped_walls += 1
                    continue
                a, b = int(row["id_i"]), int(row["id_j"])
                pairs[tuple(sorted((a, b)))].append(row)
    if len(calls) != 1:
        raise SystemExit(f"expected one force call, found {sorted(calls)}")

    stats = {"same": [0, 0.0, 0.0], "cross": [0, 0.0, 0.0]}
    bad = []
    worst = []
    for ids, rows in pairs.items():
        if len(rows) != 2 or rows[0]["id_i"] != rows[1]["id_j"] or rows[0]["id_j"] != rows[1]["id_i"]:
            bad.append((ids, len(rows)))
            continue
        x, y = rows
        work = float(x["work"]) + float(y["work"])
        px, py = float(x["pstar"]), float(y["pstar"])
        p_rel = abs(px - py) / max(abs(px), abs(py), 1e-30)
        ds = math.hypot(float(x["dSx"]) + float(y["dSx"]),
                        float(x["dSy"]) + float(y["dSy"]))
        ds_rel = ds / max(math.hypot(float(x["dSx"]), float(x["dSy"])), 1e-30)
        u_sum = float(x["udotds"]) + float(y["udotds"])
        kind = "same" if x["rank"] == y["rank"] else "cross"
        s = stats[kind]
        s[0] += 1
        s[1] += work
        s[2] += abs(work)
        worst.append((abs(work), ids, kind, work, p_rel, ds_rel, u_sum))

    print(f"call={next(iter(calls))} walls_skipped={skipped_walls} unpaired_or_duplicate={len(bad)}")
    for kind, (count, signed, absolute) in stats.items():
        print(f"{kind}: faces={count} signed_work={signed:.9e} abs_work={absolute:.9e}")
    if bad:
        print("first_unpaired_or_duplicate:", bad[:10])
    print("top |pair work|: ids kind work pstar_rel dS_rel udotds_sum")
    for _, ids, kind, work, p_rel, ds_rel, u_sum in sorted(worst, reverse=True)[:10]:
        print(ids, kind, f"{work:.9e}", f"{p_rel:.6e}",
              f"{ds_rel:.6e}", f"{u_sum:.6e}")


if __name__ == "__main__":
    main(sys.argv[1:])
