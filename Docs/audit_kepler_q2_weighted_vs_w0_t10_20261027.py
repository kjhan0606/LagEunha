"""Compare accepted-time W0 and experimental weighted native t=10 endpoints."""

import json
from pathlib import Path

from prepare_kepler_sparse_area_restart_20261027 import digest


DOCS = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027")
OUTPUT = DOCS / "KEPLER_Q2_WEIGHTED_NATIVE_V2_VS_W0_T10_2026-10-27.json"
W0 = {"aligned": ROOT / "to10_near32_a10_v1/audit_t10.json",
      "jitter": ROOT / "jitter005_to10_a10_v1/audit_t10.json"}
WEIGHTED = ROOT / "gpu_energy_neutral_weighted_t10_a10_v1"
METRICS = (
    "disk_density_abs_mass_weighted",
    "disk_entropy_log_abs_mass_weighted",
    "disk_position_error_abs_mass_weighted",
    "disk_radial_speed_rms_over_circular",
    "total_energy_relative",
)
REQUIRED_WEIGHTED = (
    "finite_positive", "fixed_mass_exact", "finite_nonnegative_variable_weights",
    "independent_power_area", "density_and_EOS", "zero_injection",
    "total_energy", "weight_stage_sequence", "zero_stage_weight_power",
    "inside_polygon_incircle", "no_wall_crossing", "fixed_polygon_enabled",
)


def read(path):
    value = json.loads(path.read_text())
    if digest(Path(value["endpoint"])) != value["endpoint_sha256"]:
        raise ValueError(f"Changed native endpoint: {path}")
    return value


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    rows = {}
    for label in ("aligned", "jitter"):
        base = read(W0[label])
        candidate = read(WEIGHTED / label / "audit_t10.json")
        if (abs(base["time"]-10.)>1e-12 or abs(candidate["time"]-10.)>1e-12
                or base["initial_sha256"] != candidate["initial_sha256"]):
            raise ValueError(f"{label}: unequal IC or physical time")
        if {name for name,good in base["checks"].items() if not good} != (
                set() if label=="aligned" else {"disk_ballistic_residual"}):
            raise ValueError(f"{label}: unexpected W0 gate status")
        if any(not candidate["checks"][name] for name in REQUIRED_WEIGHTED):
            raise ValueError(f"{label}: candidate violates a conservation or geometry prerequisite")
        comparison = {}
        for name in METRICS:
            a,b = base[name],candidate[name]
            comparison[name] = {
                "W0": a, "weighted": b,
                "W0_over_weighted": (a/b if b else None),
            }
        rows[label] = {
            "initial_sha256": base["initial_sha256"],
            "W0_endpoint_sha256": base["endpoint_sha256"],
            "weighted_endpoint_sha256": candidate["endpoint_sha256"],
            "W0_failed_checks": [name for name,good in base["checks"].items() if not good],
            "weighted_failed_checks": [name for name,good in candidate["checks"].items() if not good],
            "metrics": comparison,
            "weighted_stage_mean_seconds": candidate["mean_stage_rate_seconds"],
            "weighted_stage_max_abs_power": candidate["max_abs_stage_weight_power"],
        }
    result = {
        "time": 10., "cases": rows,
        "scope": "Same IC and exact physical time, bounded q2 controls. Weighted density is measured against its own initial field; an exact-orbit W0 ballistic-area residual is not a weighted-cell gate. This does not establish GIZMO t=600 stability.",
    }
    with OUTPUT.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
