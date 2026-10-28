"""Stage a guarded native A10 bounded Kepler continuation to exact t=3."""

import hashlib
import json
from pathlib import Path
import shutil


DOCS = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs")
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
PREVIOUS = ROOT / "to0825_near32_v1"
REFERENCE = Path("/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_bounded_r275_wallpreflight_near32_a10_20261027/khout.000000.dat")
DEST = ROOT / "to3_near32_v1"
BINARY_SHA = "3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e"
REFERENCE_SHA = "8871ce22054655461aadb6809771db4abacc4ef6e7a1cc82e4fb6547931c1355"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    receipt_path = PREVIOUS / "audit_t0825.json"
    receipt = json.loads(receipt_path.read_text())
    if abs(receipt["time"]-.825) > 1e-10 or not all(receipt["checks"].values()):
        raise ValueError("Predecessor t=.825 physics gate has not passed")
    endpoint = Path(receipt["endpoint"])
    if endpoint.parent.resolve() != PREVIOUS.resolve():
        raise ValueError("Predecessor endpoint outside accepted run")
    if digest(endpoint) != receipt["endpoint_sha256"]:
        raise ValueError("Predecessor endpoint changed after audit")
    if digest(PREVIOUS/"eunha2.exe") != BINARY_SHA or digest(REFERENCE) != REFERENCE_SHA:
        raise ValueError("Native binary or reference initial state changed")
    files = {
        "eunha2.exe": PREVIOUS/"eunha2.exe",
        "params.dat": PREVIOUS/"params.dat",
        "khout.000001.dat": endpoint,
        "reference_initial.dat": REFERENCE,
        "job.slurm": DOCS/"kepler_q2_bounded_poly1024_to3_a10_20261027.slurm",
        "audit_kepler_bounded_native_evolution_20261027.py":
            DOCS/"audit_kepler_bounded_native_evolution_20261027.py",
        "audit_kepler_sparse_area_t0825_20261027.py":
            DOCS/"audit_kepler_sparse_area_t0825_20261027.py",
        "audit_kepler_sparse_native_ballistic_20261027.py":
            DOCS/"audit_kepler_sparse_native_ballistic_20261027.py",
        "audit_kepler_sparse_state_consistency_20261027.py":
            DOCS/"audit_kepler_sparse_state_consistency_20261027.py",
        "prepare_kepler_sparse_area_restart_20261027.py":
            DOCS/"prepare_kepler_sparse_area_restart_20261027.py",
    }
    hashes = {name: digest(path) for name, path in files.items()}
    DEST.mkdir(exist_ok=False)
    for name, source in files.items():
        shutil.copy2(source, DEST/name)
        if digest(DEST/name) != hashes[name]:
            raise IOError(f"Staged {name} changed")
    (DEST/"inputs.sha256").write_text("".join(
        f"{value}  {name}\n" for name, value in sorted(hashes.items())))
    stage = {"run": str(DEST), "predecessor": str(PREVIOUS),
             "predecessor_audit_sha256": digest(receipt_path),
             "predecessor_endpoint_sha256": receipt["endpoint_sha256"],
             "target_time": 3.0, "input_sha256": hashes}
    (DEST/"staging.json").write_text(json.dumps(stage, indent=2)+"\n")
    print(json.dumps(stage, indent=2))


if __name__ == "__main__":
    main()
