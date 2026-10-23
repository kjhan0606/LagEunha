#!/usr/bin/env python3
"""Retain exact build inputs for a finished source-edge chunk."""

import hashlib
import json
import shutil
import sys
from pathlib import Path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main(directory):
    directory = directory.resolve()
    if (directory / "native_exit_code.txt").read_text().strip() != "0":
        raise ValueError("freeze only a completed native success")
    manifest = directory / "source.sha256"
    temporary = directory / "source_snapshot.tmp"
    destination = directory / "source_snapshot"
    if temporary.exists() or destination.exists():
        raise FileExistsError("source snapshot already started or completed")
    temporary.mkdir()
    rows = []
    for index, line in enumerate(manifest.read_text().splitlines()):
        expected, filename = line.split("  ", 1)
        original = Path(filename)
        if not original.is_absolute():
            original = directory.parent / original
        if digest(original) != expected:
            raise ValueError(f"build input changed before freeze: {original}")
        retained = temporary / f"{index:03d}_{original.name}"
        shutil.copyfile(original, retained)
        if digest(retained) != expected:
            raise ValueError(f"snapshot copy mismatch: {original}")
        rows.append({"source": str(original), "retained": retained.name,
                     "sha256": expected})
    temporary.rename(destination)
    receipt = {"source_manifest_sha256": digest(manifest), "files": rows}
    (directory / "source_snapshot.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(f"FROZEN {directory} inputs={len(rows)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: freeze_power_chunk_source.py RUN_DIRECTORY")
    main(Path(sys.argv[1]))
