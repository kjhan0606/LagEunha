"""Conditional native weighted-geometry gate after the finite algebra gates.

No fluid integration or production cost claim. Positions and weights are
perturbed together, and the original native tessellator measures volumes.
"""
from pathlib import Path
import argparse
import json
import struct
import subprocess

import numpy as np

from material_volume_candidate_audit import load_geometry, assemble, algebra

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("executable", type=Path)
parser.add_argument("prerequisite", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
status = json.loads((args.prerequisite / "status.json").read_text())
assert status["status"].startswith("necessary algebra/nonaffine gates pass")
args.output.mkdir()
(args.output / "registration.json").write_text(json.dumps({
    "scope": __doc__, "cases": [[d, s, w] for d, s in [(1, 16), (2, 6), (3, 6)] for w in [0, 0.02]],
    "prerequisite": str(args.prerequisite),
    "relative_displacements": [0.01, 0.003, 0.001, 0.0003, 0.0001, 0.00003],
    "gate": "Best relative volume derivative error <1e-5; improvement required unless the first error is already <1e-10 (exact linear maps reach cancellation immediately). All cells remain valid; common-weight shift changes neither geometry nor C."}, indent=2) + "\n")


def direction_file(path, dim, velocity, weight_rate):
    with path.open("xb") as stream:
        stream.write(struct.pack("<ii", dim, len(weight_rate)))
        stream.write(np.asarray(velocity, dtype="<f8").tobytes())
        stream.write(np.asarray(weight_rate, dtype="<f8").tobytes())


def perturb(directory, dim, side, fraction, direction, epsilon, name):
    command = [str(args.executable), str(dim), str(side), str(fraction), str(direction), str(epsilon)]
    (directory / (name + ".command.json")).write_text(json.dumps(command) + "\n")
    with (directory / (name + ".json")).open("x") as out, (directory / (name + ".log")).open("x") as log:
        code = subprocess.run(command, stdout=out, stderr=log).returncode
    assert code == 0, (directory, name, code)
    return json.loads((directory / (name + ".json")).read_text())


rows = []
try:
    for dim, side in [(1, 16), (2, 6), (3, 6)]:
        for fraction in [0, 0.02]:
            directory = args.output / f"D{dim}_N{side}_W{fraction:g}"
            g = load_geometry(args.executable, directory, dim, side, fraction)
            a = assemble(g)
            algebra_result = algebra(g, a)
            velocity = np.empty_like(a["x"])
            for k in range(dim):
                velocity[:, k] = 0.13 * np.sin(2 * np.pi * a["x"][:, k]) + 0.07 * np.cos(2 * np.pi * a["x"][:, (k + 1) % dim])
            rate = a["R"] @ velocity.ravel()
            rate -= rate.mean()  # Constant W-rate changes no cell volume.
            direction = directory / "direction.bin"
            direction_file(direction, dim, velocity, rate)
            expected = a["C"] @ velocity.ravel()
            scale = max(float(np.max(abs(expected))), 1e-30)
            hmin = float(np.min(a["V"])**(1 / dim))
            duration = hmin / np.linalg.norm(velocity, axis=1).max()
            measurements = []
            for k, relative in enumerate([0.01, 0.003, 0.001, 0.0003, 0.0001, 0.00003]):
                epsilon = relative * duration
                plus = perturb(directory, dim, side, fraction, direction, epsilon, f"plus_{k}")
                minus = perturb(directory, dim, side, fraction, direction, -epsilon, f"minus_{k}")
                derivative = (np.array(plus["volume"]) - np.array(minus["volume"])) / (2 * epsilon)
                error = float(np.max(abs(derivative - expected)) / scale)
                measurements.append({"relative_displacement": relative, "epsilon": epsilon, "relative_volume_derivative_error": error})
            shift = directory / "gauge_shift.bin"
            direction_file(shift, dim, np.zeros_like(velocity), np.ones(len(rate)))
            shifted = perturb(directory, dim, side, fraction, shift, 0.37, "common_weight_shift")
            shifted_operator = assemble(shifted)
            geometry_gauge_error = float(np.max(abs(np.array(shifted["volume"]) / a["V"] - 1)))
            operator_gauge_error = float(np.max(abs(shifted_operator["C"] - a["C"])))
            assert geometry_gauge_error < 1e-9 and operator_gauge_error < 1e-9
            errors = [x["relative_volume_derivative_error"] for x in measurements]
            (directory / "measurements.json").write_text(json.dumps(measurements, indent=2) + "\n")
            assert min(errors) < 1e-5
            assert errors[0] < 1e-10 or min(errors[1:]) < errors[0]
            row = {"dimension": dim, "side": side, "weight_fraction": fraction,
                   "necessary_algebra": algebra_result, "measurements": measurements,
                   "common_weight_shift_geometry_error": geometry_gauge_error,
                   "common_weight_shift_operator_error": operator_gauge_error,
                   "best_relative_volume_derivative_error": min(errors),
                   "accepted_steps": 0, "weighted_geometry_gate_pass": True}
            rows.append(row)
            print(json.dumps({key: row[key] for key in ["dimension", "side", "weight_fraction", "best_relative_volume_derivative_error"]}), flush=True)
            (directory / "result.json").write_text(json.dumps(row, indent=2) + "\n")
    result = {"status": "weighted geometry gate passes; native geometry-inclusive cost and PDE evolution unqualified", "cases": rows, "accepted_steps": 0}
except Exception as error:
    result = {"status": "candidate stopped at weighted geometry gate", "error": repr(error), "completed_cases": rows, "accepted_steps": 0}
    (args.output / "status.json").write_text(json.dumps(result, indent=2) + "\n")
    raise
(args.output / "status.json").write_text(json.dumps(result, indent=2) + "\n")
