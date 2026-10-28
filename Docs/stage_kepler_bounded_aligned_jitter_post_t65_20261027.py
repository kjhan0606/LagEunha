"""Stage a CPU-only comparison gated by both native t=6.5 jobs."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess


DOCS = Path(__file__).resolve().parent
DEST = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/aligned_jitter_post_t65_cpu_v2")
PARENTS = (437317, 437327)
RUNNER = "kepler_bounded_aligned_jitter_post_t65_cpu_20261027.slurm"
FILES = (
    RUNNER,
    "audit_kepler_bounded_radial_tails_20261027.py",
    "audit_kepler_bounded_aligned_jitter_pair_20261027.py",
    "plot_kepler_bounded_aligned_jitter_density_20261027.py",
    "plot_kepler_bounded_periodic_density_20261027.py",
    "plot_kepler_native_density_long_20261026.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_native_ballistic_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
PARSER = Path("/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/audit_snapshot.py")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    if not PARSER.is_file():
        raise FileNotFoundError(PARSER)
    DEST.mkdir()
    for name in FILES:
        shutil.copy2(DOCS / name, DEST / name)
    hashes = {name: digest(DEST / name) for name in FILES}
    with (DEST / "prestaged.sha256").open("x") as stream:
        for name, sha in sorted(hashes.items()):
            stream.write(f"{sha}  {name}\n")
        stream.write(f"{digest(PARSER)}  {PARSER}\n")
    manifest = {
        "destination": str(DEST),
        "parents": PARENTS,
        "dependency": "afterok:" + ":".join(str(x) for x in PARENTS),
        "runner": RUNNER,
        "prestaged_sha256": hashes,
        "external_snapshot_parser": {"path": str(PARSER), "sha256": digest(PARSER)},
        "physics_state": "Read only after both parent native t=6.5 auditors pass",
    }
    (DEST / "prestage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
