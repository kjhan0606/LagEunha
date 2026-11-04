"""Saved native planar FV cells versus an independently integrated exact Sod solution."""
from pathlib import Path
import argparse
import io
import json
import subprocess

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("reference", type=Path)
    args = parser.parse_args()
    records = []
    for case in sorted(args.root.glob("D*_sod_N*")):
        result = json.loads((case/"result.json").read_text())
        if result["failed"]:
            raise RuntimeError(f"Nonterminal scientific endpoint: {case}")
        final, birth = np.loadtxt(case/"final.dat"), np.loadtxt(case/"initial.dat")
        d, side, time = result["dim"], result["side"], result["time"]
        inferred = 1+birth[:, 12]*birth[:, 10]/(birth[:, 9]-.5*np.sum(birth[:, 6:9]**2, axis=1)/birth[:, 5])
        assert np.max(abs(inferred-1.4)) < 1e-12
        # Planar extrusion is checked rather than assumed from the filename.
        assert np.max(abs(final[:, 2:4]-birth[:, 2:4])) < 1e-10
        assert np.max(abs(final[:, 17:19])) < 1e-10
        transverse = float(np.max(abs(final[:, 14:16])))
        assert transverse < 1e-9
        width = final[:, 10]*side**(d-1)
        center = final[:, 1]+final[:, 16]
        left, right = center-.5*width, center+.5*width
        assert left.min() > -1e-11 and right.max() < 1+1e-11
        intervals = "".join(f"{a:.17g} {b:.17g}\n" for a, b in zip(left, right))
        process = subprocess.run([str(args.reference), str(time), "1.4"], input=intervals,
                                 text=True, capture_output=True, check=True)
        reference = np.loadtxt(io.StringIO(process.stdout), delimiter=",")
        norm = np.sum(final[:, 10, None]*abs(final[:, [11, 13, 12]]-reference[:, :3]), axis=0)
        record = {"case": case.name, "reference": "Native exact Riemann solution integrated separately across all waves on actual planar cell intervals",
                  "gamma": 1.4, "time": time, "density_L1": float(norm[0]),
                  "mass_mean_velocity_L1": float(norm[1]), "physical_mean_pressure_L1": float(norm[2]),
                  "transverse_velocity_max": transverse, "planar_extrusion_checked": True,
                  "pressure_interpretation": "Actual physical pressure average; continuum subcell kinetic variance is not relabelled as thermal energy."}
        (case/"exact_cell_reference.csv").write_text(process.stdout)
        (case/"sod_accuracy.json").write_text(json.dumps(record, indent=2)+"\n")
        records.append(record)
    (args.root/"sod_accuracy.json").write_text(json.dumps(records, indent=2)+"\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
