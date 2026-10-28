"""Reconstruct the native W=0 pressure-force onset from existing q2 snapshots."""

import json
import os
import subprocess
import sys
from pathlib import Path

from audit_kepler_sparse_area_t0825_20261027 import ordered
from prepare_kepler_sparse_area_restart_20261027 import digest


ROOT = Path(__file__).resolve().parent
CASES = ROOT / "KEPLER_Q2_THREE_IC_NATIVE_T65_2026-10-27.json"
AUDITOR = ROOT / "audit_kepler_bounded_pressure_force_t10_20261027.py"
OUTPUT = ROOT / "KEPLER_Q2_INITIAL_PRESSURE_SERIES_2026-10-27"
TARGETS = (0.825, 2.475, 4.126, 5.775, 6.5)


def main():
    OUTPUT.mkdir(exist_ok=False)
    manifest = json.loads(CASES.read_text())
    rows = []
    for case, item in manifest["cases"].items():
        receipt = Path(item["source_audit"])
        if digest(receipt) != item["source_audit_sha256"]:
            raise ValueError(f"Changed accepted audit for {case}")
        audit = json.loads(receipt.read_text())
        directory = Path(audit["endpoint"]).parent
        snapshots = []
        for entry in os.scandir(directory):
            if entry.name.startswith("khout.") and entry.name.endswith(".dat"):
                _, meta = ordered(Path(entry.path), 55351)
                snapshots.append((float(meta["t"]), Path(entry.path)))
        for target in TARGETS:
            when, checkpoint = min(snapshots, key=lambda row: abs(row[0]-target))
            if abs(when-target) > 0.06:
                raise ValueError(f"No near-{target} checkpoint for {case}")
            result_file = OUTPUT / f"{case}_t{target:.3f}.json"
            command = [sys.executable, str(AUDITOR), str(receipt),
                       "--checkpoint", str(checkpoint), "--output", str(result_file)]
            subprocess.run(command, check=True)
            result = json.loads(result_file.read_text())
            if result["source_snapshot_sha256"] != digest(checkpoint):
                raise ValueError(f"Changed source snapshot for {case} {target}")
            disk = result["whole_disk"]
            rows.append({"case":case,"target_time":target,"actual_time":when,
                         "source_snapshot":str(checkpoint),
                         "source_sha256":result["source_snapshot_sha256"],
                         "pressure_relative_rms":disk["pressure_relative_rms"],
                         "acceleration_mass_weighted_rms":disk["acceleration_mass_weighted_rms"],
                         "tangential_acceleration_mass_weighted_rms":disk["tangential_acceleration_mass_weighted_rms"],
                         "area_max_relative_difference":result["reference_area_max_relative_difference_from_native"]})
            print(f"{case} t={when:.6f} acceleration={disk['acceleration_mass_weighted_rms']:.6e}",flush=True)
    summary={"source_three_ic":str(CASES),"source_three_ic_sha256":digest(CASES),
             "auditor":str(AUDITOR),"auditor_sha256":digest(AUDITOR),"rows":rows,
             "scope":"Frozen-state pressure acceleration at five saved times in each accepted t6.5 run. The intermediate snapshots have not each passed an independent full native audit; this is not an RK4 force-history integral or a causal force-off run."}
    with (OUTPUT/"summary.json").open("x") as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
