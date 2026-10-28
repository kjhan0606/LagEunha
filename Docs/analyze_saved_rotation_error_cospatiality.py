"""Compare saved intrinsic-spin errors with density/pressure errors by cell.

Read-only join of existing actual-cell steady-reference reports and saved
rotation diagnostics at t=.025/.05. Co-location/correlation is descriptive,
not a causal attribution or a physics update.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from compare_full_moment_checkpoint_errors import load


CASES = (
    ("t0025_coarse", "coarse", "coarse_0005"),
    ("t0025_fine", "fine", "fine_0010"),
    ("t005_coarse", "coarse", "coarse_0010"),
    ("t005_fine", "fine", "fine_0020"),
)
REGIONS = ("cavity", "inner_contact", "disk_core")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def weighted_corr(a, b, weight):
    a = np.asarray(a, dtype=np.longdouble)
    b = np.asarray(b, dtype=np.longdouble)
    w = np.asarray(weight, dtype=np.longdouble)
    w /= np.sum(w)
    a -= np.sum(w*a)
    b -= np.sum(w*b)
    aa = np.sum(w*a*a)
    bb = np.sum(w*b*b)
    return None if aa == 0 or bb == 0 else float(np.sum(w*a*b)/np.sqrt(aa*bb))


def summarize_region(ids, mass, mask, spin_error, density_error, pressure_error):
    indices = np.flatnonzero(mask)
    if not len(indices):
        return {"cells": 0}
    w = np.asarray(mass[indices], dtype=np.longdouble)
    spin = np.asarray(spin_error[indices], dtype=np.longdouble)
    density = np.asarray(density_error[indices], dtype=np.longdouble)
    pressure = np.asarray(pressure_error[indices], dtype=np.longdouble)
    scores = {
        "spin": w*spin*spin,
        "density": w*density*density,
        "pressure": w*pressure*pressure,
    }
    top = {}
    sets = {}
    for name, score in scores.items():
        k = max(1, int(math.ceil(.1*len(indices))))
        local = np.argsort(score)[-k:][::-1]
        sets[name] = set(int(x) for x in ids[indices[local]])
        top[name] = dict(
            cells=k,
            ids=ids[indices[local]].astype(int).tolist(),
            fraction_of_region_squared_error=float(np.sum(score[local])/np.sum(score))
            if np.sum(score)>0 else None,
        )
    overlaps = {}
    for a, b in (("spin", "density"), ("spin", "pressure"), ("density", "pressure")):
        intersection=sets[a]&sets[b];union=sets[a]|sets[b]
        overlaps[f"top10pct_{a}_{b}"]=dict(
            intersection_cells=len(intersection),
            union_cells=len(union),
            jaccard=float(len(intersection)/len(union)) if union else None,
            shared_ids=sorted(intersection),
        )
    return dict(
        cells=len(indices),
        mass_weighted_correlation=dict(
            signed_spin_error_vs_signed_density_error=weighted_corr(spin,density,w),
            signed_spin_error_vs_signed_pressure_error=weighted_corr(spin,pressure,w),
            absolute_spin_error_vs_absolute_density_error=weighted_corr(abs(spin),abs(density),w),
            absolute_spin_error_vs_absolute_pressure_error=weighted_corr(abs(spin),abs(pressure),w),
        ),
        top10pct_squared_error_contributors=top,
        top10pct_overlap=overlaps,
    )


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root",type=Path,default=Path(
        "/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019"))
    parser.add_argument("--rotation-root",type=Path,default=Path("Docs"))
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"refusing to overwrite {args.output}")
    cases={}
    for tag,branch,checkpoint_tag in CASES:
        folder=args.run_root/branch/"checkpoints"/checkpoint_tag
        fields,report,source=load(folder/"physics")
        rot_dir=args.rotation_root/f"canonical_saved_cell_rotation_{tag}_2026-10-19"
        with np.load(rot_dir/"rotation.npz") as z:
            rotation={key:z[key].copy() for key in z.files}
        ids=fields["ids"]
        if not np.array_equal(ids,rotation["ids"]):
            raise SystemExit(f"ID mismatch in {tag}")
        spin_error=rotation["actual_ell"]-rotation["reference_ell"]
        physical_spin_error=fields["actual_intrinsic_ell"]-fields["reference_intrinsic_ell"]
        spin_crosscheck_max_abs=float(abs(spin_error-physical_spin_error).max())
        if spin_crosscheck_max_abs>1e-18:
            raise SystemExit(f"spin diagnostic disagrees with physical comparison in {tag}")
        pressure_error=fields["relative_pressure_error"]
        density_error=fields["relative_density_reference_error"]
        cases[tag]=dict(
            time=float(report["time"]),
            rotation_error_crosscheck_max_absolute=spin_crosscheck_max_abs,
            source=dict(physical_comparison=source,
                        rotation_report=str((rot_dir/"report.json").resolve()),
                        rotation_report_sha256=sha256(rot_dir/"report.json"),
                        rotation_payload_sha256=sha256(rot_dir/"rotation.npz")),
            regions={name:summarize_region(ids,fields["actual_mass"],fields[name],
                spin_error,density_error,pressure_error) for name in REGIONS},
        )
    result=dict(
        status="completed_saved_rotation_error_cospatiality",
        scope="Four saved endpoints; cellwise steady-reference spin errors joined to density/pressure errors.",
        spin_error_definition="actual intrinsic ell minus the steady Kepler reference ell on that endpoint's actual moving cell",
        physical_error_definitions=dict(
            density="relative_density_reference_error from the accepted physical comparison",
            pressure="relative_pressure_error from the accepted physical comparison"),
        cases=cases,
        limits=[
            "Mass-weighted correlation and top-cell overlap do not prove causal coupling.",
            "Each endpoint's Kepler reference is evaluated on its own actual moving-cell geometry.",
            "No RHS, provider, accepted step, or physical model change was made.",
        ],
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    brief={tag:{region:row for region,row in case["regions"].items()}
           for tag,case in cases.items()}
    print(json.dumps({"status":result["status"],"cases":brief},sort_keys=True))


if __name__=="__main__":
    main()
