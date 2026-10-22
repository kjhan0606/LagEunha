#!/usr/bin/env python3
"""Lightweight saved-output audit only; no physical recomputation."""
import csv
import hashlib
import json
import math
import struct
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNS = {n: ROOT / f"geometric_pressure_v2_mpi{n}" for n in (1, 4)}
SAVED = ROOT / "shared_trace_rhs_v1_mpi4"

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

def frozen(run, state, ranks):
    hashes = {}
    counts = {}
    header = struct.Struct("=6Q")
    row = struct.Struct("=Q5d")
    for rank in range(ranks):
        with (run / state / f"frozen_physical_rank{rank}.bin").open("rb") as stream:
            h = header.unpack(stream.read(header.size))
            require(h[:3] == (0x47465346524F5A32, 2, 1) and h[4:] == (row.size, rank), "frozen thermal header")
            count = 0
            while True:
                block = stream.read(row.size * 8192)
                if not block:
                    break
                require(len(block) % row.size == 0, "complete frozen physical rows")
                for offset in range(0, len(block), row.size):
                    raw = block[offset:offset + row.size]
                    tri, X, Y, weight, rho, U = row.unpack(raw)
                    require(tri % ranks == rank and all(map(math.isfinite, (X, Y, weight, rho, U))) and min(weight, rho, U) > 0, "one-owner positive frozen physical point")
                    hashes.setdefault(tri, hashlib.sha256()).update(raw)
                    counts[tri] = counts.get(tri, 0) + 1
                    count += 1
            require(count == h[3], "frozen physical point count")
    require(len(hashes) == 384 and sum(counts.values()) == 817152, "same original positive rule")
    return {tri: (counts[tri], digest.hexdigest()) for tri, digest in hashes.items()}

ledger = {"scope": "two immutable SAVED-state native controls, no physical evolution", "runs": {}, "comparisons": {}}
for ranks, run in RUNS.items():
    require((run / "native_exit_code.txt").read_text().strip() == "0", "native exit success")
    for name in ("final_input_integrity.txt", "final_source_integrity.txt"):
        lines = (run / name).read_text().splitlines()
        require(lines and all(line.endswith(": OK") for line in lines), "input/source SHA preservation")
    result = json.loads((run / "contract.json").read_text())
    require(result["same_state_bundle_pass"] and result["actual_mass_equation_residual_max"] < 1e-16, "complete native bundle")
    states = {}
    for state, original in (("initial", "initial_state.bin"), ("finite", "shared_state.bin")):
        data = json.loads((run / state / "result.json").read_text())
        require(data["same_state_bundle_pass"], "physical state gate")
        require((run / state / "shared_state.bin").read_bytes() == (SAVED / original).read_bytes(), "x/v/basis/owner checkpoint remains bitwise original")
        constants = table(run / state / "constant_pressure_native.csv")
        require(len(constants) == 3 and all(max(float(v) for k, v in row.items() if k != "constant") < 1e-11 for row in constants), "actual EOS/native stress/assembly constant-pressure gate")
        for gravity in (0, 2):
            g = json.loads((run / state / f"gravity_{gravity}.json").read_text())
            require(g["energy_rate_error_over_U"] < 1e-8 and g["angular_rate_error_over_2K"] < 1e-11, "same-state gravity ledger")
        cells = table(run / state / "cell_physics_rotation.csv")
        require(len(cells) == 24 and all(float(c["inertia"]) > 0 and float(c["rotation_energy"]) >= 0 for c in cells), "actual resolved rotation")
        # Finite result's zero chi is NOT a measurement: its IC curvature
        # formula is deliberately applied only at the initial state.
        data["initial_curvature_evaluated"] = state == "initial"
        if state == "finite":
            data["max_chi"] = None
            data["every_cell_mean_chi_max"] = None
        states[state] = {"native_result": data, "constant_pressure_tests": constants, "frozen_digest": frozen(run, state, ranks)}
    initial = states["initial"]["native_result"]
    require(initial["max_chi"] <= .002 and initial["every_cell_mean_chi_max"] <= .001, "unchanged initial curvature criteria")
    ledger["runs"][str(ranks)] = {"summary": result, "states": states}

for state in ("initial", "finite"):
    one = ledger["runs"]["1"]["states"][state]
    four = ledger["runs"]["4"]["states"][state]
    require(one["frozen_digest"] == four["frozen_digest"], "MPI1/4 original physical m/U point table bitmatch")
    require((RUNS[1] / state / "shared_coefficients.bin").read_bytes() == (RUNS[4] / state / "shared_coefficients.bin").read_bytes(), "same SAVED x/v including signs of zero")
    a, b = table(RUNS[1] / state / "acceleration.csv"), table(RUNS[4] / state / "acceleration.csv")
    maximum, scale = Decimal(0), Decimal(0)
    with localcontext() as context:
        context.prec = 50
        for x, y in zip(a, b):
            require(x["global_id"] == y["global_id"] and x["owner"] == y["owner"], "MPI canonical acceleration IDs")
            for field in ("ax", "ay", "Fx", "Fy"):
                xx, yy = Decimal(x[field]), Decimal(y[field])
                maximum = max(maximum, abs(xx - yy))
                scale = max(scale, abs(xx), abs(yy))
        relative = maximum / scale
    require(len(a) == len(b) == 3797 and relative < Decimal("1e-11"), "MPI1/4 mass/load result agreement")
    ledger["comparisons"][state] = {"physical_m_U_reference_rows_bitwise_equal": True, "saved_x_v_coefficients_bitwise_equal": True, "acceleration_load_combined_max_difference": str(maximum), "acceleration_load_combined_relative_difference": str(relative)}
    # Per-source point hashes are an audit detail; keep an aggregate digest.
    for entry in (one, four):
        encoded = json.dumps(entry.pop("frozen_digest"), sort_keys=True).encode()
        entry["frozen_physical_digest"] = hashlib.sha256(encoded).hexdigest()

ledger["all_same_state_native_checks_pass"] = True
ledger["whole_target_transfer_certified"] = False
ledger["production_time_step_or_long_Kepler_certified"] = False
output = ROOT / "GEOMETRIC_PRESSURE_LEDGER_20261021.json"
with output.open("x") as stream:
    json.dump(ledger, stream, indent=2)
    stream.write("\n")
print(json.dumps({"audit": "PASS", "ledger": str(output), "MPI_comparisons": ledger["comparisons"]}))
