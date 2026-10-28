"""Audit optional RK4 timing rows against full KH evolution wall timing."""

import json
from pathlib import Path
import re
import statistics
import sys


def main():
    log = Path(sys.argv[1]).read_text(errors="replace")
    output = Path(sys.argv[2])
    pattern = re.compile(
        r"\[RK_FINE\] step=(\d+) wall=([0-9.]+) named=([0-9.]+) "
        r"stage_diag=([0-9.]+) tx_valid=([0-9.]+) residual=([0-9.]+)")
    rows = [tuple(float(value) for value in match.groups()[1:])
            for match in pattern.finditer(log)]
    steps = [int(match.group(1)) for match in pattern.finditer(log)]
    kh = [float(value) for value in re.findall(
        r"\[KH_WALL\] step=\d+ t=[0-9.eE+-]+ evol=([0-9.]+)s", log)]
    if len(rows) != 14 or len(kh) != 14 or len(set(steps)) != 14:
        raise ValueError("Expected 14 unique RK_FINE and KH_WALL rows")
    names = ("wall", "named", "stage_diag", "tx_valid", "residual")
    median = {name: statistics.median(row[index] for row in rows)
              for index, name in enumerate(names)}
    checks = {
        "all_components_nonnegative": all(all(value >= -1e-7 for value in row)
                                          for row in rows),
        "components_close_each_step": all(
            abs(row[0] - sum(row[1:])) < 1e-6 for row in rows),
        "rk_wall_matches_kh_evolution": all(
            abs(row[0] - value) < .1 for row, value in zip(rows, kh)),
        "native_run_completed": "EXIT:0" in log,
    }
    result = {"steps": steps, "median_seconds": median,
              "median_kh_evolution_seconds": statistics.median(kh),
              "stage_diag_wall_fraction": median["stage_diag"] / median["wall"],
              "tx_valid_wall_fraction": median["tx_valid"] / median["wall"],
              "residual_wall_fraction": median["residual"] / median["wall"],
              "checks": checks,
              "scope": "One-rank q1 t=.001 profiling only; residual includes RK stage arithmetic, other checks, and setup"}
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not all(checks.values()):
        raise SystemExit("RK fine profile gate failed")


if __name__ == "__main__":
    main()
