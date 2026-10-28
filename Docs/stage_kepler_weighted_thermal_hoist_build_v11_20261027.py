"""Build an isolated weighted v10 plus RK4 thermal64-check hoist candidate."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
DOCS = Path(__file__).resolve().parent
PARENT = ROOT / "gpu_energy_neutral_wallselect_cache_capacity64_build_v10"
RUN = ROOT / "gpu_energy_neutral_wallselect_thermalhoist_capacity64_build_v11"
THERMAL_PATCH = "KEPLER_RK4_DIAG_THERMAL_HOIST_2026-10-27.patch"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parent = json.loads((PARENT / "stage_manifest.json").read_text())
    template = (PARENT / "job.slurm").read_text()
    if digest(PARENT / "job.slurm") != parent["inputs"]["job.slurm"]:
        raise ValueError("v10 build runner changed")
    replacements = {
        "LagQ2W64WallBuild": "LagQ2W64ThermBuild",
        "lag_q2_wall_capacity64_": "lag_q2_thermal_capacity64_",
        "--partition=a10": "--partition=h200",
        "--gres=gpu:0": "--gres=gpu:H200:1",
        "for patch_name in accel.patch wallselect.patch; do":
            "for patch_name in accel.patch wallselect.patch thermal.patch; do",
        "eunha2_wallselect_cache_capacity64_v10.exe":
            "eunha2_wallselect_cache_thermal_capacity64_v11.exe",
    }
    for old, new in replacements.items():
        expected_count = 2 if old.startswith("eunha2_") else 1
        if template.count(old) != expected_count:
            raise ValueError(f"Unexpected v10 runner anchor: {old!r}")
        template = template.replace(old, new)

    names = (
        "source.tar", "patch_guard.py", "patch_iterations.py", "patch_zero.py",
        "patch_capacity.py", "gfs_pressure_neutral_projection_20261027.h",
        "test_projection.cpp", "accel.patch", "wallselect.patch",
    )
    sources = {name: PARENT / name for name in names}
    sources["thermal.patch"] = DOCS / THERMAL_PATCH
    for name in names:
        if digest(sources[name]) != parent["inputs"][name]:
            raise ValueError(f"Changed v10 input: {name}")
    if RUN.exists():
        raise FileExistsError(RUN)

    RUN.mkdir()
    (RUN / "job.slurm").write_text(template)
    subprocess.run(["bash", "-n", str(RUN / "job.slurm")], check=True)
    hashes = {"job.slurm": digest(RUN / "job.slurm")}
    for name, source in sources.items():
        destination = RUN / name
        if name == "source.tar":
            destination.symlink_to(source.resolve())
        else:
            shutil.copy2(source, destination)
        hashes[name] = digest(destination)
        if hashes[name] != digest(source):
            raise ValueError(f"Changed staged input: {name}")
    with (RUN / "inputs.sha256").open("x") as output:
        for name, value in sorted(hashes.items()):
            output.write(f"{value}  {name}\n")
    with (RUN / "stage_manifest.json").open("x") as output:
        json.dump({"parent": str(PARENT), "inputs": hashes,
                   "scope": "Weighted v10 plus the diagnostic-only thermal64 check hoist; "
                            "the new binary requires a same-state native pair gate."},
                  output, indent=2)
        output.write("\n")
    print(RUN)


if __name__ == "__main__":
    main()
