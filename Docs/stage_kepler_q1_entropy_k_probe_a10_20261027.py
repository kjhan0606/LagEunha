"""Stage two isolated q1 t=.825→.845 native A10 entropy-mode controls.

The baseline executable is the accepted q1 binary. The entropy executable
differs only by the opt-in K initialization patch. No run is submitted here.
"""

import hashlib
import json
from pathlib import Path
import shutil
import struct

import numpy as np


N = 384 * 384
STRIDE = 576
K_OFFSET = 552  # C offsetof(treevorostressrk4particletype, stress.K)
OLD_BIN_SHA = "7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284"
NEW_BIN_SHA = "f836416a58b4743c7b81d348c54a365f03c5bafc0e76c1cc492d8d95182f40c0"
PARAMS_SHA = "2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43"
CHECKPOINT_SHA = "23fa751f5a52f27f8e8d17b51d8c85a8d334ce0e1aa57f56f5fbf2678053f7c7"
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_EntropyK_A10_20261027")
Q1 = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1")
PARENT = Q1 / "halo12_from061876_to0825"
CHECKPOINT = PARENT / "khout.010856.dat"
DOCS = Path(__file__).resolve().parent
NEW_BINARY = Path("/tmp/Kepler_entropy_k_20261027/eunha2.exe")
SCRIPTS = (
    "kepler_q1_entropy_k_to0845_a10_20261027.slurm",
    "audit_kepler_polar_long_gate_20261027.py",
    "audit_kepler_native_radial_20261026.py",
)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require_snapshot(path, time):
    if path.stat().st_size != 16 + N * STRIDE + 52:
        raise ValueError(f"Incomplete snapshot: {path}")
    with path.open("rb") as stream:
        magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
        stream.seek(-52, 2)
        actual = struct.unpack("<d", stream.read(8))[0]
    if (magic, n, stride, stress) != (-1, N, STRIDE, 1) or abs(actual - time) > 1e-10:
        raise ValueError(f"Incompatible snapshot time or format: {path}")


def sha_file(path):
    return f"{digest(path)}  {path.name}\n"


def main():
    if ROOT.exists():
        raise FileExistsError(f"Refusing to reuse staged root: {ROOT}")
    receipt = json.loads((PARENT / "audit_long_t0825_strict.json").read_text())
    if (not all(receipt["checks"].values()) or receipt["n"] != N
            or abs(receipt["time"] - .825) > 1e-10
            or Path(receipt["endpoint"]).resolve() != CHECKPOINT.resolve()):
        raise ValueError("q1 accepted checkpoint strict audit is not valid")
    require_snapshot(Q1 / "khout.000000.dat", 0.)
    require_snapshot(CHECKPOINT, .825)
    if digest(CHECKPOINT) != CHECKPOINT_SHA:
        raise ValueError("Accepted checkpoint hash changed")
    if digest(PARENT / "eunha2.exe") != OLD_BIN_SHA:
        raise ValueError("Accepted q1 binary hash changed")
    if digest(NEW_BINARY) != NEW_BIN_SHA:
        raise ValueError("Experimental binary hash changed")
    if digest(PARENT / "params.dat") != PARAMS_SHA:
        raise ValueError("Accepted params hash changed")
    raw = np.memmap(CHECKPOINT, dtype="u1", mode="r", offset=16,
                    shape=(N * STRIDE,))
    K = np.ndarray((N,), dtype="<f8", buffer=raw, offset=K_OFFSET,
                   strides=(STRIDE,))
    if not np.isfinite(K).all() or np.count_nonzero(K):
        raise ValueError("Accepted energy-mode checkpoint has nonzero K")

    params = (PARENT / "params.dat").read_text()
    old_line = "define GAS entropy_mode   = 0"
    if params.count(old_line) != 1:
        raise ValueError("Cannot uniquely select native entropy_mode parameter")
    ROOT.mkdir(parents=True)
    manifest = {
        "parent_strict_receipt": str(PARENT / "audit_long_t0825_strict.json"),
        "parent_checkpoint": str(CHECKPOINT),
        "parent_checkpoint_sha256": CHECKPOINT_SHA,
        "parent_K_nonzero": 0,
        "checkpoint_K_offset": K_OFFSET,
        "initial": str(Q1 / "khout.000000.dat"),
        "start_time": .825, "target_time": .845,
        "old_binary_sha256": OLD_BIN_SHA, "new_binary_sha256": NEW_BIN_SHA,
        "source_patch": str(DOCS / "KEPLER_ENTROPY_K_INIT_2026-10-27.patch"),
        "source_patch_sha256": digest(DOCS / "KEPLER_ENTROPY_K_INIT_2026-10-27.patch"),
        "scope": "Matched native q1 short attribution probe; no long-time acceptance claim",
        "runs": {},
    }
    for label, mode, binary in (
        ("energy", 0, PARENT / "eunha2.exe"),
        ("entropy", 1, NEW_BINARY),
    ):
        run = ROOT / label
        run.mkdir()
        for name in SCRIPTS:
            shutil.copy2(DOCS / name, run / name)
        shutil.copy2(binary, run / "eunha2.exe")
        shutil.copy2(CHECKPOINT, run / CHECKPOINT.name)
        (run / "params.dat").write_text(
            params.replace(old_line, f"define GAS entropy_mode   = {mode}"))
        (run / "expected_mode.txt").write_text(f"{mode}\n")
        (run / "scripts.sha256").write_text(
            "".join(sha_file(run / name) for name in SCRIPTS))
        (run / "inputs.sha256").write_text(
            "".join(sha_file(run / name) for name in
                    ("eunha2.exe", CHECKPOINT.name, "params.dat", "expected_mode.txt")))
        (run / SCRIPTS[0]).chmod(0o755)
        manifest["runs"][label] = {
            "path": str(run), "entropy_mode": mode,
            "binary_sha256": digest(run / "eunha2.exe"),
            "params_sha256": digest(run / "params.dat"),
            "scripts_sha256": digest(run / "scripts.sha256"),
            "inputs_sha256": digest(run / "inputs.sha256"),
        }
    with (ROOT / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
