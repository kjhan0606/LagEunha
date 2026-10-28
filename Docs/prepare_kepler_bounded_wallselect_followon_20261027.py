"""Stage the fast CPU and selective-wall CUDA binary from an accepted state."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys


PARENT_BINARY_SHA = "d92b6cd3b8491678cdfafa7610608378db89ca9bd16585a52a51aaf26f3b3367"
FAST_BINARY = Path(
    "/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/"
    "wallselect_pair_t58_a10_v1/candidate/eunha2.exe"
)
FAST_BINARY_SHA = "87d60c1cb059e1710dd2749c84af5c98ba17fe077e6fe080b9c8b42d5fd40968"
COPY_FILES = (
    "params.dat", "glass.bin", "areas.bin",
    "audit_kepler_bounded_native_evolution_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("parent_audit")
    parser.add_argument("parent_time", type=float)
    parser.add_argument("runner")
    args = parser.parse_args()
    parent = args.parent.resolve()
    dest = Path.cwd().resolve()
    # The parser is staged from the accepted parent only after this gate runs.
    # Import its parent copy for the native completeness check before copying.
    sys.path.insert(0, str(parent))
    from prepare_kepler_sparse_area_restart_20261027 import snapshot_time
    if (dest / "khout.000001.dat").exists() or (dest / "inputs.sha256").exists():
        raise FileExistsError("Follow-on stage already contains restart input")
    audit_path = parent / args.parent_audit
    audit = json.loads(audit_path.read_text())
    if abs(float(audit["time"]) - args.parent_time) > 1e-10 or not all(audit["checks"].values()):
        raise ValueError("Parent physics gate did not accept the requested state")
    endpoint = Path(audit["endpoint"]).resolve()
    initial = Path(audit["initial"]).resolve()
    if endpoint.parent != parent or not endpoint.is_file():
        raise ValueError("Accepted endpoint is not a complete file in parent run")
    if digest(endpoint) != audit["endpoint_sha256"] or digest(initial) != audit["initial_sha256"]:
        raise ValueError("Accepted snapshot SHA differs from parent physics audit")
    if abs(snapshot_time(endpoint, 55351) - args.parent_time) > 1e-10:
        raise ValueError("Accepted endpoint has wrong or incomplete native time")
    if digest(parent / "eunha2.exe") != PARENT_BINARY_SHA:
        raise ValueError("Parent executable differs from accepted bounded binary")
    if digest(FAST_BINARY) != FAST_BINARY_SHA:
        raise ValueError("Faster executable differs from bitwise-tested binary")
    if not (dest / args.runner).is_file():
        raise FileNotFoundError("Submitted runner is missing")
    for name in COPY_FILES:
        target = dest / name
        if target.exists():
            if digest(target) != digest(parent / name):
                raise FileExistsError(f"Refusing different existing file: {target}")
        else:
            shutil.copy2(parent / name, target)
    shutil.copy2(FAST_BINARY, dest / "eunha2.exe")
    shutil.copy2(endpoint, dest / "khout.000001.dat")
    manifest = {
        "parent": str(parent),
        "parent_audit": str(audit_path),
        "parent_audit_sha256": digest(audit_path),
        "accepted_parent_time": args.parent_time,
        "parent_endpoint": str(endpoint),
        "parent_endpoint_sha256": audit["endpoint_sha256"],
        "initial": str(initial),
        "initial_sha256": audit["initial_sha256"],
        "parent_binary_sha256": PARENT_BINARY_SHA,
        "binary_source": str(FAST_BINARY),
        "binary_sha256": FAST_BINARY_SHA,
        "runner": args.runner,
    }
    with (dest / "restart_origin.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    inputs = (*COPY_FILES, "eunha2.exe", "khout.000001.dat", args.runner,
              "prepare_kepler_bounded_wallselect_followon_20261027.py", "restart_origin.json")
    with (dest / "inputs.sha256").open("x") as stream:
        for name in inputs:
            stream.write(f"{digest(dest / name)}  {name}\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
