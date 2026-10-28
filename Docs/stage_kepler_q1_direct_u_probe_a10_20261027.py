"""Stage one isolated native q1 direct-U arithmetic run from accepted t=.825."""

import hashlib
import json
from pathlib import Path
import shutil
import struct


N = 384 * 384
STRIDE = 576
ACCEPTED_BINARY_SHA = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"
PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
CHECKPOINT_SHA = "23fa751f5a52f27f8e8d17b51d8c85a8d334ce0e1aa57f56f5fbf2678053f7c7"
BASE = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1")
PARENT = BASE / "halo12_from061876_to0825"
CHECKPOINT = PARENT / "khout.010856.dat"
DESTINATION = Path("/gpfs/kjhan/LagForce/Kepler_DirectU_A10_20261027")
DOCS = Path(__file__).resolve().parent
SCRIPTS = (
    "kepler_q1_direct_u_to0845_a10_20261027.slurm",
    "audit_kepler_polar_long_gate_20261027.py",
    "audit_kepler_native_radial_20261026.py",
    "audit_kepler_q1_direct_u_anchor_20261027.py",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check_snapshot(path, expected_time):
    if path.stat().st_size != 16 + N * STRIDE + 52:
        raise ValueError(f"Incomplete snapshot: {path}")
    with path.open("rb") as stream:
        head = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        actual_time = struct.unpack("<d", stream.read(8))[0]
    if head != (-1, N, STRIDE, 1) or abs(actual_time - expected_time) > 1e-10:
        raise ValueError(f"Wrong snapshot format or time: {path}")


def main():
    if DESTINATION.exists():
        raise FileExistsError(f"Refusing to reuse {DESTINATION}")
    receipt = json.loads((PARENT / "audit_long_t0825_strict.json").read_text())
    if (not all(receipt["checks"].values()) or
            Path(receipt["endpoint"]).resolve() != CHECKPOINT.resolve() or
            abs(receipt["time"] - .825) > 1e-10):
        raise ValueError("Accepted parent strict gate is invalid")
    check_snapshot(BASE / "khout.000000.dat", 0.)
    check_snapshot(CHECKPOINT, .825)
    if (digest(CHECKPOINT) != CHECKPOINT_SHA or
            digest(PARENT / "eunha2.exe") != ACCEPTED_BINARY_SHA or
            digest(PARENT / "params.dat") != PARAMS_SHA):
        raise ValueError("Accepted source checkpoint or native inputs changed")
    params = (PARENT / "params.dat").read_text()
    if params.count("define GAS entropy_mode   = 0") != 1:
        raise ValueError("Direct-U requires energy-mode q1 params")
    DESTINATION.mkdir(parents=True)
    for name in SCRIPTS:
        shutil.copy2(DOCS / name, DESTINATION / name)
    for source in (PARENT / "eunha2.exe", PARENT / "params.dat", CHECKPOINT):
        shutil.copy2(source, DESTINATION / source.name)
    (DESTINATION / SCRIPTS[0]).chmod(0o755)
    (DESTINATION / "scripts.sha256").write_text("".join(
        f"{digest(DESTINATION / name)}  {name}\n" for name in SCRIPTS))
    (DESTINATION / "inputs.sha256").write_text("".join(
        f"{digest(DESTINATION / name)}  {name}\n" for name in
        ("eunha2.exe", "params.dat", CHECKPOINT.name)))
    report = {
        "destination": str(DESTINATION),
        "accepted_start": str(CHECKPOINT), "start_sha256": CHECKPOINT_SHA,
        "start_time": .825, "target_time": .845,
        "binary_sha256": ACCEPTED_BINARY_SHA, "params_sha256": PARAMS_SHA,
        "mode": "GFS_RK4_DIRECT_U=1, entropy_mode=0, same native RK4 equations",
        "anchor_offset": 472,
        "scripts_sha256": digest(DESTINATION / "scripts.sha256"),
        "inputs_sha256": digest(DESTINATION / "inputs.sha256"),
        "scope": "Short arithmetic attribution probe only; no t>6 or performance claim",
    }
    with (DESTINATION / "manifest.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
