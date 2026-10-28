"""Stage a same-binary raw/projected q2 pair from the accepted t=0 IC."""

import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
IC = ROOT / "jitter005_t001_a10_v1"
PILOT = ROOT / "gpu_energy_neutral_weighted_fresh_jitter64_t001_a10_v1"
BUILD = ROOT / "gpu_energy_neutral_raw_rate_capacity64_build_v11b"
INITIAL_SHA = "8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b"
BINARY_SHA = "5d10fd4c17ac67bfafcc461fda03937bc645159758d1fff1ab73bd27ecf54f72"
MASK_SHA = "b560aca44c81e4560f8371c6d979e1afbe853f8cd18b592d0fc19a58becc8f4c"
OLD_BINARY_SHA = "7358c54996b820a7b45bd53bcf6f788f125754b3158845d2b4265afead570f38"


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Unexpected job template anchor: {old}")
    return text.replace(old, new)


def main():
    old_receipt = json.loads((PILOT / "audit_t001.json").read_text())
    if (old_receipt["time"] != .001 or len(old_receipt["checks"]) != 18
            or not all(old_receipt["checks"].values())
            or old_receipt["reference_sha256"] != INITIAL_SHA
            or digest(IC / "khout.000000.dat") != INITIAL_SHA
            or digest(BUILD / "eunha2_raw_rate_v11.exe") != BINARY_SHA
            or digest(PILOT / "disk_active_mask.bin") != MASK_SHA):
        raise ValueError("Accepted t=0 IC, mask, pilot, or diagnostic binary changed")
    for mode in ("projected", "raw"):
        if (ROOT / f"gpu_energy_neutral_v11b_fresh_{mode}_t001_a10_v1").exists():
            raise FileExistsError(mode)

    sources = {
        "eunha2.exe": BUILD / "eunha2_raw_rate_v11.exe",
        "glass.bin": IC / "glass.bin",
        "areas.bin": IC / "areas.bin",
        "params.dat": IC / "params.dat",
        "accepted_initial.dat": IC / "khout.000000.dat",
        "ic_manifest.json": IC / "ic_manifest.json",
        "disk_active_mask.bin": PILOT / "disk_active_mask.bin",
    }
    for name in (
        "audit_kepler_q2_weighted_fresh_pilot_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py",
        "audit_kepler_sparse_state_consistency_20261027.py",
        "laguerre_nonlinear_volume_probe_t6_20261026.py",
        "prepare_kepler_sparse_area_restart_20261027.py",
    ):
        sources[name] = DOCS / name
    template = (DOCS / "kepler_q2_weighted_fresh_jitter64_t001_20261027.slurm").read_text()

    for mode, projection in (("projected", 1), ("raw", 0)):
        dest = ROOT / f"gpu_energy_neutral_v11b_fresh_{mode}_t001_a10_v1"
        dest.mkdir()
        hashes = {}
        for name, source in sources.items():
            shutil.copy2(source, dest / name)
            hashes[name] = digest(dest / name)
            if hashes[name] != digest(source):
                raise ValueError(f"Changed copy: {name}")
        runner = replace_once(template, "LagQ2W64Fresh001", f"LagQ2Fresh{mode[:3]}001")
        runner = replace_once(runner, OLD_BINARY_SHA, BINARY_SHA)
        runner = replace_once(runner, "export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64",
                              "export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64 "
                              f"GFS_GPU_VOLW2_PROJECTION={projection}")
        (dest / "job.slurm").write_text(runner)
        subprocess.run(["bash", "-n", str(dest / "job.slurm")], check=True)
        hashes["job.slurm"] = digest(dest / "job.slurm")
        with (dest / "inputs.sha256").open("x") as output:
            for name, sha in sorted(hashes.items()):
                output.write(f"{sha}  {name}\n")
        with (dest / "stage_manifest.json").open("x") as output:
            json.dump({"run": str(dest), "mode": mode, "projection": projection,
                       "accepted_initial_sha256": INITIAL_SHA,
                       "source_pilot_receipt": str(PILOT / "audit_t001.json"),
                       "binary_sha256": BINARY_SHA, "mask_sha256": MASK_SHA,
                       "target_time": .001, "inputs": hashes,
                       "scope": "Isolated same-binary native t=0→.001 causal comparison; "
                                "does not replace the ongoing later raw/projected runs."},
                      output, indent=2)
            output.write("\n")
        print(dest)


if __name__ == "__main__":
    main()
