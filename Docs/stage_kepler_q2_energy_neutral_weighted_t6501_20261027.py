"""Stage a native one-step weighted replay from both accepted t=6.5 states.

All inputs are checked before creating a new destination. Existing runs are
never overwritten. The weighted mode remains an experimental diagnostic.
"""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
BUILD = ROOT / "gpu_energy_neutral_guard_build_v2"
BINARY = BUILD / "eunha2_energy_neutral_guard.exe"
RUNNER = "kepler_q2_energy_neutral_weighted_t6501_a10_20261027.slurm"
EXPECTED_PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
CASES = {
    "aligned": ("to65_near32_a10_v1", "KEPLER_Q2_BOUNDED_NATIVE_T65_2026-10-27.json"),
    "jitter": ("jitter005_to65_a10_v1", "KEPLER_Q2_BOUNDED_JITTER005_NATIVE_T65_2026-10-27.json"),
}
HELPERS = (
    "audit_kepler_q2_weighted_short_replay_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "laguerre_nonlinear_volume_probe_t6_20261026.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)


def copy_inputs(destination, sources):
    destination.mkdir(exist_ok=True)
    result = {}
    for name, source in sources.items():
        shutil.copy2(source, destination / name)
        result[name] = digest(destination / name)
    with (destination / "inputs.sha256").open("x") as output:
        for name, value in sorted(result.items()):
            output.write(f"{value}  {name}\n")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    destination = args.destination.resolve()
    if destination.exists():
        raise FileExistsError(destination)
    subprocess.run(["bash", "-n", str(DOCS / RUNNER)], check=True)
    binary_record = (BUILD / "binary.sha256").read_text().split()[0]
    if digest(BINARY) != binary_record:
        raise ValueError("Built executable no longer matches its receipt")
    build_record = json.loads((BUILD / "guard_receipt.json").read_text())
    if (build_record["original_sha256"] !=
            "5ebc8d575bbb8bb361069923384b450682f950b36625c82c7396e50cb1c33f67"
            or build_record["patched_sha256"] !=
            "75e2a75b3069052407d6c7ff20514850fe18872d5227430334b0673ab4b6164f"):
        raise ValueError("Wrong guarded native source build")
    mask_record = json.loads(
        (DOCS / "KEPLER_Q2_BOUNDED_NATIVE_WDOT_PROBE_INPUT_T65_2026-10-27.json").read_text()
    )
    mask = Path(mask_record["mask"])
    if mask.stat().st_size != 55351 or digest(mask) != mask_record["mask_sha256"]:
        raise ValueError("Changed 49,038-cell disk activity mask")
    accepted = {}
    for label, (subdir, receipt_name) in CASES.items():
        receipt = json.loads((DOCS / receipt_name).read_text())
        source = ROOT / subdir
        endpoint = Path(receipt["endpoint"])
        if (receipt["time"] != 6.5 or len(receipt["checks"]) != 17
                or not all(receipt["checks"].values())
                or endpoint.parent != source
                or digest(endpoint) != receipt["endpoint_sha256"]
                or abs(snapshot_time(endpoint, 55351) - 6.5) > 1e-10
                or digest(source / "params.dat") != EXPECTED_PARAMS_SHA):
            raise ValueError(f"{label}: accepted parent changed")
        accepted[label] = (source, endpoint, receipt)
    destination.mkdir(parents=True)
    root_inputs = copy_inputs(destination, {
        "disk_active_mask.bin": mask,
        "job.slurm": DOCS / RUNNER,
    })
    staged = {}
    for label, (source, endpoint, _) in accepted.items():
        staged[label] = copy_inputs(destination / label, {
            "eunha2.exe": BINARY,
            "params.dat": source / "params.dat",
            "khout.000001.dat": endpoint,
            **{name: DOCS / name for name in HELPERS},
        })
    manifest = {
        "run": str(destination),
        "mode": "experimental fixed-mass pressure-work-neutral weighted RK4",
        "start_time": 6.5,
        "target_time": 6.501,
        "binary_sha256": binary_record,
        "guard_source_sha256": build_record["patched_sha256"],
        "mask_sha256": mask_record["mask_sha256"],
        "parent_sha256": {label: accepted[label][2]["endpoint_sha256"] for label in CASES},
        "root_inputs": root_inputs,
        "case_inputs": staged,
        "scope": "One native weighted step only; not a release or long-term Kepler result.",
    }
    with (destination / "stage_manifest.json").open("x") as output:
        json.dump(manifest, output, indent=2)
        output.write("\n")
    print(json.dumps({"run": str(destination), "binary_sha256": binary_record,
                      "parent_sha256": manifest["parent_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
