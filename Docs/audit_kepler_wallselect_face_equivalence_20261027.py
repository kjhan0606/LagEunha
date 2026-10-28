"""Compare directed physical faces after changing only fixed-wall preclips."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered


FIELDS = ("x", "y", "vx", "vy", "mass", "volume", "ie", "den", "pressure", "w2")
FACE_NAME = "face_dump_r0_call1.tsv"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def face_groups(path):
    groups = {}
    row_count = 0
    total_length = 0.0
    tiny_rows = 0
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"id_i", "id_j", "wall", "dSx", "dSy"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("Missing native face columns")
        for row in reader:
            key = (int(row["id_i"]), int(row["id_j"]), int(row["wall"]))
            sx, sy = float(row["dSx"]), float(row["dSy"])
            length = math.hypot(sx, sy)
            if not math.isfinite(length) or length <= 0:
                raise ValueError(f"Nonpositive or nonfinite physical face: {key}")
            current = groups.get(key)
            if current is None:
                groups[key] = [1, length, sx, sy]
            else:
                current[0] += 1
                current[1] += length
                current[2] += sx
                current[3] += sy
            row_count += 1
            total_length += length
            tiny_rows += length <= 1e-14
    return groups, {"rows": row_count, "unique_directed_supports": len(groups),
                    "total_directed_face_length": total_length,
                    "rows_length_at_most_1e-14": tiny_rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    paths = {name: root / name / FACE_NAME for name in ("control", "candidate")}
    audits = {name: json.loads((root / name / "audit_t0826.json").read_text())
              for name in paths}
    if audits["control"]["restart_sha256"] != audits["candidate"]["restart_sha256"]:
        raise ValueError("Different accepted restart snapshots")
    groups, summaries = {}, {}
    for name, path in paths.items():
        groups[name], summaries[name] = face_groups(path)
        summaries[name]["path"] = str(path)
        summaries[name]["sha256"] = digest(path)
    a, b = groups["control"], groups["candidate"]
    keys_a, keys_b = set(a), set(b)
    common = keys_a & keys_b
    max_length_difference = max((abs(a[k][1]-b[k][1]) for k in common), default=0.0)
    max_vector_difference = max((math.hypot(a[k][2]-b[k][2], a[k][3]-b[k][3])
                                 for k in common), default=0.0)
    states = {name: ordered(Path(audits[name]["endpoint"]), 55351)[0]
              for name in paths}
    field_equal = {name: bool(np.array_equal(states["control"][name],
                                             states["candidate"][name]))
                   for name in FIELDS}
    checks = {
        "both_native_physics_gates": all(all(audit["checks"].values()) for audit in audits.values()),
        "same_directed_support_keys": keys_a == keys_b,
        "same_integrated_face_lengths": max_length_difference < 1e-12,
        "same_integrated_face_area_vectors": max_vector_difference < 1e-12,
        "same_physical_endpoint_fields_bitwise": all(field_equal.values()),
    }
    result = {
        "run": str(root), "restart_sha256": audits["control"]["restart_sha256"],
        "faces": summaries,
        "support_keys_only_control": len(keys_a-keys_b),
        "support_keys_only_candidate": len(keys_b-keys_a),
        "same_support_different_segment_count": sum(a[k][0] != b[k][0] for k in common),
        "max_aggregated_length_difference": max_length_difference,
        "max_aggregated_area_vector_difference": max_vector_difference,
        "field_bitwise_equal": field_equal,
        "checks": checks,
        "scope": ("One accepted-state native force call. Segments are grouped by "
                  "directed owner/neighbor/wall support; this proves equivalence "
                  "here, not arbitrary late-time geometry."),
    }
    rendered = json.dumps(result, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as stream:
        stream.write(rendered)
    print(rendered)
    if not all(checks.values()):
        raise SystemExit("Selective wall preclip face-equivalence gate failed")


if __name__ == "__main__":
    main()
