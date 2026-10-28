"""Stage an isolated accepted-weighted one-step first-ring face dump."""

import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
SOURCE = ROOT/"gpu_energy_neutral_capacity64_v5_fresh_localwork_t6501_ada_v1"
DEST = ROOT/"gpu_energy_neutral_capacity64_v5_first_ring_faces_t6501_ada_v1"
RECEIPT = DOCS/"KEPLER_FRESH_WEIGHT_LOCAL_T6501_2026-10-27/audit_t6501.json"
TEMPLATE = DOCS/"kepler_q2_weighted_fresh_localwork_t0826_ada_20261027.sh"
NAMES = ("eunha2.exe", "params.dat", "disk_active_mask.bin", "initial.dat",
         "khout.000001.dat", "audit_kepler_q2_weighted_equilibrium_20261027.py",
         "audit_kepler_bounded_native_evolution_20261027.py",
         "audit_kepler_sparse_area_t0825_20261027.py",
         "audit_kepler_sparse_native_ballistic_20261027.py",
         "audit_kepler_sparse_state_consistency_20261027.py",
         "laguerre_nonlinear_volume_probe_t6_20261026.py",
         "prepare_kepler_sparse_area_restart_20261027.py")


def main():
    if DEST.exists():
        raise FileExistsError(DEST)
    audit = json.loads(RECEIPT.read_text())
    source_manifest = json.loads((SOURCE/"stage_manifest.json").read_text())
    if (audit["time"] != 6.501 or audit["parent_time"] != 6.5
            or len(audit["checks"]) != 21 or not all(audit["checks"].values())
            or audit["binary_sha256"] != source_manifest["binary_sha256"]
            or audit["parent_sha256"] != source_manifest["parent_endpoint_sha256"]):
        raise ValueError("Original accepted one-step audit is not intact")
    for key, sha_key in (("initial", "initial_sha256"),
                         ("parent", "parent_sha256"),
                         ("endpoint", "endpoint_sha256")):
        if digest(Path(audit[key])) != audit[sha_key]:
            raise ValueError(f"Changed {key} snapshot")
    for name in NAMES:
        if digest(SOURCE/name) != source_manifest["inputs"][name]:
            raise ValueError(f"Changed source input: {name}")

    p0, _ = ordered(Path(audit["initial"]), audit["n"])
    start, _ = ordered(Path(audit["parent"]), audit["n"])
    r0 = np.hypot(p0["x"].astype(float)-3., p0["y"].astype(float)-3.)
    ring = np.flatnonzero((r0 >= .5) & (r0 < .52))
    if len(ring) != 294:
        raise ValueError("Unexpected first-ring site count")
    relative_pressure = np.abs(start["pressure"].astype(float)
                               / p0["pressure"].astype(float)-1.)
    high = ring[np.argsort(relative_pressure[ring], kind="stable")[-32:]]
    rest = np.setdiff1d(ring, high)
    angle = np.arctan2(p0["y"].astype(float)-3.,
                       p0["x"].astype(float)-3.)
    ordered_rest = rest[np.argsort(angle[rest])]
    uniform = ordered_rest[np.linspace(0, len(rest)-1, 32, dtype=int)]
    selected = np.sort(np.concatenate((high, uniform)))
    if len(selected) != 64 or len(np.unique(selected)) != 64:
        raise ValueError("Failed to select 64 distinct first-ring IDs")

    runner = TEMPLATE.read_text()
    replacements = {
        "HYDRO_TSTOP=0.826": "HYDRO_TSTOP=6.501",
        "--parent-time 0.825 --time 0.826": "--parent-time 6.5 --time 6.501",
        "--output audit_t0826.json": "--output audit_t6501.json",
        "unset GFS_GPU_WDOT_PROBE_MASK GFS_GPU_FACE_DUMP_CALL GFS_GPU_FACE_DUMP_IDS GFS_GPU_TESS_DUMP":
            ("unset GFS_GPU_WDOT_PROBE_MASK GFS_GPU_TESS_DUMP "
             "GFS_GPU_FACE_DUMP_FIRST GFS_GPU_FACE_DUMP_LAST GFS_GPU_FACE_DUMP_WALL_ONLY\n"
             "export GFS_GPU_FACE_DUMP_CALL=1\n"
             "read -r GFS_GPU_FACE_DUMP_IDS < selected_ids.txt\n"
             "export GFS_GPU_FACE_DUMP_IDS"),
        "if [ \"$run_rc\" -ne 0 ]; then exit \"$run_rc\"; fi":
            ("if [ \"$run_rc\" -ne 0 ]; then exit \"$run_rc\"; fi\n"
             "test -s face_dump_r0_call1.tsv\n"
             "test -s face_input_r0_call1.tsv\n"
             "test -s face_input_meta_r0_call1.json"),
    }
    for old, new in replacements.items():
        if runner.count(old) != 1:
            raise ValueError(f"Unexpected runner template anchor: {old}")
        runner = runner.replace(old, new)
    DEST.mkdir()
    hashes = {}
    for name in NAMES:
        shutil.copy2(SOURCE/name, DEST/name)
        hashes[name] = digest(DEST/name)
        if hashes[name] != source_manifest["inputs"][name]:
            raise ValueError(f"Copy mismatch: {name}")
    (DEST/"selected_ids.txt").write_text(
        ",".join(str(int(i)) for i in selected)+"\n")
    (DEST/"job_ada.sh").write_text(runner)
    subprocess.run(["bash", "-n", str(DEST/"job_ada.sh")], check=True)
    for name in ("selected_ids.txt", "job_ada.sh"):
        hashes[name] = digest(DEST/name)
    with (DEST/"inputs.sha256").open("x") as stream:
        for name, sha in sorted(hashes.items()):
            stream.write(f"{sha}  {name}\n")
    manifest = {"destination": str(DEST), "source_receipt": str(RECEIPT),
                "source_receipt_sha256": digest(RECEIPT),
                "source_endpoint_sha256": audit["endpoint_sha256"],
                "parent_endpoint_sha256": audit["parent_sha256"],
                "runner_template_sha256": digest(TEMPLATE),
                "selected_ids": [int(i) for i in selected],
                "selected_high_pressure_ids": [int(i) for i in high],
                "selected_uniform_angle_ids": [int(i) for i in uniform],
                "inputs": hashes,
                "scope": "One accepted weighted step on Ada with read-only first-RK4-call face dump for 64 first-ring sites."}
    with (DEST/"stage_manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(DEST)


if __name__ == "__main__":
    main()
