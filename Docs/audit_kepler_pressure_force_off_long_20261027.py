"""Evaluate the t=10 pressure-force-off causal control without accepting it as hydro."""

import argparse
import json
from pathlib import Path

from prepare_kepler_sparse_area_restart_20261027 import digest


SAFE_CHECKS=("finite_positive","fixed_mass_exact","W0_exact",
             "initial_area_tiles_polygon","endpoint_area_tiles_polygon",
             "density_and_EOS","zero_injection","inside_polygon_incircle",
             "no_wall_crossing","fixed_polygon_enabled")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("run",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    run=args.run.resolve()
    stage=json.loads((run/"staging.json").read_text())
    short=json.loads(Path(stage["one_step_pair"]).read_text())
    original=json.loads(Path(stage["original_t10_receipt"]).read_text())
    diagnostic=json.loads((run/"audit_t10_native.json").read_text())
    if (not all(short["checks"].values())
            or short["accepted_parent_sha256"]!=stage["accepted_parent_sha256"]
            or diagnostic["time"]!=10 or original["time"]!=10
            or digest(run/"eunha2.exe")!=stage["diagnostic_binary_sha256"]
            or digest(run/"khout.000001.dat")!=stage["accepted_parent_sha256"]
            or digest(Path(diagnostic["endpoint"]))!=diagnostic["endpoint_sha256"]
            or digest(Path(original["endpoint"]))!=original["endpoint_sha256"]
            or diagnostic["initial_sha256"]!=original["initial_sha256"]):
        raise ValueError("Causal control or original reference changed")
    log=(run/"log").read_text(errors="replace")
    if ("[DIAG_PRESSURE_FORCE_OFF] pressure force suppressed" not in log
            or "EXIT:0" not in log):
        raise ValueError("Pressure-force-off native evolution did not finish cleanly")
    safe={key:bool(diagnostic["checks"].get(key,False)) for key in SAFE_CHECKS}
    keys=("disk_density_abs_mass_weighted",
          "disk_ballistic_density_abs_mass_weighted",
          "disk_native_minus_ballistic_abs_mass_weighted",
          "disk_entropy_log_abs_mass_weighted",
          "disk_position_error_abs_mass_weighted",
          "disk_radial_speed_rms_over_circular",
          "disk_angular_momentum_relative",
          "total_energy_relative")
    comparison={key:{"original":float(original[key]),
                     "pressure_force_off":float(diagnostic[key]),
                     "diagnostic_over_original_magnitude":(
                         abs(float(diagnostic[key]))/abs(float(original[key]))
                         if original[key] else None)} for key in keys}
    result={"time":10.0,"stage":str(run/"staging.json"),
            "one_step_pair":stage["one_step_pair"],
            "one_step_pressure_acceleration_ratio":short["difference_acceleration_over_frozen_pressure_force_rms"],
            "original_t10_receipt":stage["original_t10_receipt"],
            "original_t10_failed_checks":[k for k,v in original["checks"].items() if not v],
            "diagnostic_native_receipt":str(run/"audit_t10_native.json"),
            "diagnostic_native_failed_checks":[k for k,v in diagnostic["checks"].items() if not v],
            "accepted_parent_sha256":stage["accepted_parent_sha256"],
            "diagnostic_endpoint_sha256":diagnostic["endpoint_sha256"],
            "safety_checks":safe,"metric_comparison":comparison,
            "hypothesis_density_residual_reduced":bool(
                diagnostic["disk_native_minus_ballistic_abs_mass_weighted"]
                <original["disk_native_minus_ballistic_abs_mass_weighted"]),
            "scope":"Causal diagnostic with native pressure force suppressed from accepted t6.5 jitter state. Geometric PdV/EOS retained, but physical momentum-energy adjoint broken; this is not an accepted hydro solution or a method to use at long time."}
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2))
    if not all(safe.values()):
        raise SystemExit("Pressure-force-off safety gate failed")


if __name__=="__main__":
    main()
