"""Stage two independent one-step A10 restarts for GPU kernel-A work profiling.

Run on syntax. All paths are explicit, and existing output directories are
refused. The accepted clipcompact binary and the instrumented binary get the
same t=2.55 checkpoint, parameters, physical settings and Slurm script.
"""

import hashlib
import json
import shutil
from pathlib import Path


DOCS = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs")
PARENT = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/") / (
    "q2_profile_wide_area_t001/to0p825"
)
ACCEPTED = PARENT / "clipcompact_from2517696_a10"
PROBE_BIN = Path("/tmp/Kepler_sparse_workprofile_20261027/eunha2.exe")
RUNNER = "kepler_sparse_tesswork_one_step_a10_20261027.slurm"
SCRIPTS = (
    "audit_kepler_sparse_area_t10_20261027.py",
    "audit_kepler_sparse_state_consistency_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)
CHECKED = {
    ACCEPTED / "khout.001349.dat":
        "503dad7e489fdbf86b8869d3297ddb343675d6d200e73423a59d95e8420fea46",
    ACCEPTED / "params.dat":
        "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43",
    ACCEPTED / "eunha2.exe":
        "886638100a4ae057ccd6c3624c1b712141dea0656ccaa791934d0f00d44834fc",
    PROBE_BIN:
        "6bf4af2b5fb09d524aca6e2a5b90625e1c7ac278d58a3a92b99e0b889e65aebc",
}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    outputs = {
        "baseline": PARENT / "tessbase_t2551_a10",
        "instrumented": PARENT / "tesswork_t2551_a10",
    }
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(path)
    for path, expected in CHECKED.items():
        actual = digest(path)
        if actual != expected:
            raise ValueError(f"SHA mismatch: {path}: {actual} != {expected}")
    for path in (DOCS / RUNNER, *(DOCS / name for name in SCRIPTS)):
        if not path.is_file():
            raise FileNotFoundError(path)

    for name, out in outputs.items():
        out.mkdir()
        for source in (
            ACCEPTED / "khout.001349.dat",
            ACCEPTED / "params.dat",
            DOCS / RUNNER,
            *(DOCS / script for script in SCRIPTS),
        ):
            shutil.copy2(source, out / source.name)
        binary = (ACCEPTED / "eunha2.exe") if name == "baseline" else PROBE_BIN
        shutil.copy2(binary, out / "eunha2.exe")
        staged = (
            "khout.001349.dat", "params.dat", "eunha2.exe", RUNNER, *SCRIPTS
        )
        with (out / "inputs.sha256").open("x") as stream:
            for filename in staged:
                stream.write(f"{digest(out / filename)}  {filename}\n")
        receipt = {
            "case": name,
            "run_dir": str(out),
            "checkpoint_time": 2.55,
            "stop_time": 2.551,
            "binary_source": str(binary),
            "binary_sha256": digest(out / "eunha2.exe"),
            "checkpoint_sha256": digest(out / "khout.001349.dat"),
            "runner_sha256": digest(out / RUNNER),
            "scope": "One-step GPU diagnostic; no long-time Kepler evidence.",
        }
        with (out / "stage_receipt.json").open("x") as stream:
            json.dump(receipt, stream, indent=2)
            stream.write("\n")
        print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
