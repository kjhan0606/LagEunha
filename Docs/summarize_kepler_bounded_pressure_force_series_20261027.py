"""Pair aligned/jitter frozen-state pressure-force measurements by time."""

import json
from pathlib import Path


INDICES=("000053","000478","000904","001329","001754","001806")


def main():
    rows=[]
    for idx in INDICES:
        a=json.loads(Path(f"aligned_{idx}.json").read_text())
        j=json.loads(Path(f"jitter_{idx}.json").read_text())
        if abs(a["time"]-j["time"])>1e-8:
            raise ValueError(f"Different physical times at {idx}")
        if a["native_failed_checks"] or j["native_failed_checks"] != ["disk_ballistic_residual"]:
            raise ValueError("Final native failure labels changed")
        aw,jw=a["whole_disk"],j["whole_disk"]
        rows.append({"snapshot_index":idx,"time":a["time"],
                     "aligned_sha256":a["source_snapshot_sha256"],
                     "jitter_sha256":j["source_snapshot_sha256"],
                     "aligned_force_rms":aw["acceleration_mass_weighted_rms"],
                     "jitter_force_rms":jw["acceleration_mass_weighted_rms"],
                     "jitter_over_aligned_force_rms":jw["acceleration_mass_weighted_rms"]/aw["acceleration_mass_weighted_rms"],
                     "aligned_tangential_force_rms":aw["tangential_acceleration_mass_weighted_rms"],
                     "jitter_tangential_force_rms":jw["tangential_acceleration_mass_weighted_rms"],
                     "aligned_pressure_relative_rms":aw["pressure_relative_rms"],
                     "jitter_pressure_relative_rms":jw["pressure_relative_rms"]})
    if any(rows[k]["time"] >= rows[k+1]["time"] for k in range(len(rows)-1)):
        raise ValueError("Snapshot times are not increasing")
    result={"rows":rows,
            "source":"Frozen-state force reconstruction with t6.5 64-cell GPU cross-check; t10 jitter parent fails disk_ballistic_residual",
            "scope":"Intermediate checkpoints have no separate 17-check native gate. Force snapshots show growth but do not isolate the causal mechanism of integrated density error."}
    with Path("pressure_force_series_t10.json").open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    for row in rows:
        print(row["time"],row["aligned_force_rms"],row["jitter_force_rms"],row["jitter_over_aligned_force_rms"])


if __name__=="__main__":
    main()
