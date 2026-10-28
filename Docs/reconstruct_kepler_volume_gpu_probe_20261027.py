"""Reconstruct the accepted bounded-Kepler source in an isolated worktree."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS = Path(__file__).resolve().parent
BASE = DOCS / "KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27"
BOUNDED = DOCS / "KEPLER_FIXEDPOLY_NEAR32_SOURCE_2026-10-27"
COMMIT = "472e7fcc3419b060d8fa0c74f0961f6c1f9f188b"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def copy_checked(source, dest, expected):
    if digest(source) != expected:
        raise ValueError(f"Archived source hash changed: {source}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, dest)
    if digest(dest) != expected:
        raise ValueError(f"Copied source hash changed: {dest}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("worktree", type=Path)
    args = ap.parse_args()
    dest = args.worktree.resolve()
    if subprocess.check_output(["git", "-C", str(dest), "rev-parse", "HEAD"],
                               text=True).strip() != COMMIT:
        raise ValueError("Wrong base Git commit")
    if (dest / ".kepler_bounded_overlay.json").exists():
        raise FileExistsError("This worktree was already overlaid")
    base = json.loads((BASE / "MANIFEST.json").read_text())
    if base["base_git_commit"] != COMMIT:
        raise ValueError("Base archive manifest differs from worktree")
    records = []
    for entry in base["source_files"]:
        relative = Path(entry["path"])
        if not relative.parts or relative.parts[0] != "Exam" or ".." in relative.parts:
            raise ValueError(f"Unexpected overlay path {relative}")
        copy_checked(BASE / relative, dest / relative, entry["sha256"])
        records.append({"path": str(relative), "sha256": entry["sha256"],
                        "overlay": "base"})
    bounded = json.loads((BOUNDED / "MANIFEST.json").read_text())
    for relative_string, expected in bounded["source_sha256"].items():
        relative = Path(relative_string)
        if not relative.parts or relative.parts[0] != "Exam" or ".." in relative.parts:
            raise ValueError(f"Unexpected bounded overlay path {relative}")
        copy_checked(BOUNDED / relative, dest / relative, expected)
        records.append({"path": str(relative), "sha256": expected,
                        "overlay": "bounded-near32"})
    result = {"worktree": str(dest), "base_commit": COMMIT,
              "accepted_binary_sha256": bounded["binary_sha256"],
              "archived_source_records": records}
    with (dest / ".kepler_bounded_overlay.json").open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"worktree": str(dest), "base_records": len(base["source_files"]),
                      "bounded_records": len(bounded["source_sha256"]),
                      "binary_sha256": bounded["binary_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
