"""Stage one immutable H200 fixed-polygon wall/preflight gate on syntax."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil


PROJECT = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce")
DOCS = PROJECT / "Docs"
PARAMS = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/params.dat")
RUN_ROOT = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026")
EXPECTED_BINARY = "3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
EXPECTED_GLASS = "37fcc5f1600bca9083b90aa2a24f98f973acff517266fe623d07e43fb92558c2"
EXPECTED_AREAS = "37583489abdfd5fa8617cb4cc74e5dae7e42fd55b40619852e5a6e426726ec20"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path,
                        default=Path("/tmp/Kepler_sparse_fixedpoly_20261027/eunha2.exe"))
    parser.add_argument("--name", default="q2_bounded_r275_wallpreflight_h200_20261027")
    parser.add_argument("--gpu", choices=("h200", "a10"), default="h200")
    args = parser.parse_args()
    if "/" in args.name or args.name in ("", ".", ".."):
        raise ValueError("Run name must be a new single directory name")
    inputs = {
        "eunha2.exe": args.binary,
        "params.dat": PARAMS,
        "glass.bin": DOCS / "KEPLER_Q2_BOUNDED_R275_N1024_2026-10-27.glass.bin",
        "areas.bin": DOCS / "KEPLER_Q2_BOUNDED_R275_N1024_2026-10-27.areas.bin",
        "job.slurm": DOCS / f"kepler_q2_bounded_wallpreflight_{args.gpu}_20261027.slurm",
        "audit_kepler_bounded_native_preflight_20261027.py":
            DOCS / "audit_kepler_bounded_native_preflight_20261027.py",
        "audit_kepler_fixedpoly_wall_faces_20261027.py":
            DOCS / "audit_kepler_fixedpoly_wall_faces_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py":
            DOCS / "audit_kepler_sparse_area_t0825_20261027.py",
        "audit_kepler_sparse_state_consistency_20261027.py":
            DOCS / "audit_kepler_sparse_state_consistency_20261027.py",
        "prepare_kepler_sparse_area_restart_20261027.py":
            DOCS / "prepare_kepler_sparse_area_restart_20261027.py",
    }
    expected = {"eunha2.exe": EXPECTED_BINARY,
                "glass.bin": EXPECTED_GLASS, "areas.bin": EXPECTED_AREAS}
    original_hashes = {name: sha(path) for name, path in inputs.items()}
    for name, digest in expected.items():
        if original_hashes[name] != digest:
            raise ValueError(f"Accepted {name} SHA mismatch: {original_hashes[name]}")
    destination = RUN_ROOT / args.name
    destination.mkdir(exist_ok=False)
    staged_hashes = {}
    for name, source in inputs.items():
        target = destination / name
        shutil.copyfile(source, target)
        target.chmod(source.stat().st_mode & 0o777)
        staged_hashes[name] = sha(target)
        if staged_hashes[name] != original_hashes[name]:
            raise IOError(f"Copy changed {name}")
    (destination / "inputs.sha256").write_text("".join(
        f"{digest}  {name}\n" for name, digest in sorted(staged_hashes.items())))
    receipt = {"run": str(destination), "source_binary": str(args.binary),
               "source_params": str(PARAMS), "input_sha256": staged_hashes,
               "purpose": "one-rank native fixed-polygon t=0.001 and first-call wall-face gate"}
    (destination / "staging.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
