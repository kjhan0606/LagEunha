"""Stage reproducible density rendering for the same-IC t≈8.252 comparison."""
from pathlib import Path
import shutil
import subprocess

from prepare_kepler_sparse_area_restart_20261027 import digest

DOCS = Path(__file__).resolve().parent
RUN = Path(
    "/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/"
    "fresh_weighted_density_step904_t8252_cpu_v1"
)
FILES = (
    "kepler_fresh_weighted_density_t8252_20261027.slurm",
    "plot_kepler_q2_weighted_vs_w0_density_20261027.py",
    "audit_kepler_sparse_area_t0825_20261027.py",
    "prepare_kepler_sparse_area_restart_20261027.py",
)


def main():
    if RUN.exists():
        raise FileExistsError(RUN)
    subprocess.run(["bash", "-n", str(DOCS / FILES[0])], check=True)
    RUN.mkdir()
    hashes = {}
    for name in FILES:
        source = DOCS / name
        destination = RUN / name
        shutil.copy2(source, destination)
        hashes[name] = digest(destination)
        if hashes[name] != digest(source):
            raise ValueError(f"Changed staged file: {name}")
    with (RUN / "inputs.sha256").open("x") as output:
        for name, sha256 in sorted(hashes.items()):
            output.write(f"{sha256}  {name}\n")
    print(RUN)


if __name__ == "__main__":
    main()
