"""Stage an immutable v5-based build with a diagnostic projection switch."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
V5 = ROOT / "gpu_energy_neutral_capacity64_build_v5"
DEST = ROOT / "gpu_energy_neutral_raw_rate_capacity64_build_v11b"
PARENT_NAMES = ("source.tar", "patch_guard.py", "patch_iterations.py",
                "patch_zero.py", "patch_capacity.py",
                "gfs_pressure_neutral_projection_20261027.h", "test_projection.cpp")


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    parent = json.loads((V5 / "stage_manifest.json").read_text())
    sources = {name: V5 / name for name in PARENT_NAMES}
    sources.update({
        "patch_raw_rate.py": DOCS / "apply_kepler_weight_projection_toggle_20261027.py",
        "test_raw_rate.cpp": DOCS / "test_kepler_weight_projection_toggle_20261027.cpp",
        "job.slurm": DOCS / "build_kepler_weighted_raw_rate_v11_20261027.slurm",
    })
    for name in PARENT_NAMES:
        if digest(sources[name]) != parent["inputs"][name]["sha256"]:
            raise ValueError(f"Changed v5 input: {name}")
    subprocess.run(["bash", "-n", str(sources["job.slurm"])], check=True)
    DEST.mkdir()
    hashes = {}
    for name, source in sources.items():
        if name == "source.tar":
            (DEST / name).symlink_to(source.resolve())
        else:
            shutil.copy2(source, DEST / name)
        hashes[name] = digest(DEST / name)
        if hashes[name] != digest(source):
            raise ValueError(f"Staged input changed: {name}")
    with (DEST / "inputs.sha256").open("x") as out:
        for name, sha in sorted(hashes.items()):
            out.write(f"{sha}  {name}\n")
    with (DEST / "stage_manifest.json").open("x") as out:
        json.dump({"destination": str(DEST), "base": str(V5),
                   "inputs": hashes,
                   "scope": "Isolated optional raw Jacobi weight rate; "
                            "default global energy projection remains v5."},
                  out, indent=2)
        out.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()
