"""Run a separately built native Athena++ reference and compare saved FV data.

Python schedules native C++ executables and integrates piecewise-constant
output overlaps. It does not evolve a fluid or change an accepted checkpoint.
The reference source preparation and its exact upstream commit are registered
in the caller-owned run directory before this program is invoked.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

import numpy as np


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_case(root, count):
    case = root / f"N{count}"
    case.mkdir(exist_ok=False)
    text = f"""<job>
problem_id = ShuReference
<output1>
file_type = tab
variable = prim
data_format = %24.16e
dt = 1.8
<output2>
file_type = hst
dt = 1.8
<time>
cfl_number = 0.3
nlim = -1
tlim = 1.8
integrator = rk3
xorder = 3
ncycle_out = 10000
<mesh>
nx1 = {count}
x1min = -5
x1max = 5
ix1_bc = outflow
ox1_bc = outflow
nx2 = 1
x2min = 0
x2max = 1
ix2_bc = periodic
ox2_bc = periodic
nx3 = 1
x3min = 0
x3max = 1
ix3_bc = periodic
ox3_bc = periodic
<hydro>
gamma = 1.4
"""
    (case / "athinput").write_text(text)
    binary = root / "athena_source/bin/athena"
    command = [str(binary), "-i", "athinput"]
    write_json(case / "command.json", command)
    start = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS="1", TMPDIR=str(root))
    with (case / "stdout.log").open("x") as out, (case / "stderr.log").open("x") as err:
        process = subprocess.Popen(command, cwd=case, env=env, stdout=out, stderr=err)
        write_json(case / "running.json", {"pid": process.pid, "command": command})
        code = process.wait()
    record = {"N": count, "exit_code": code, "wall_seconds": time.monotonic() - start,
              "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest()}
    write_json(case / "process.json", record)
    if code:
        raise RuntimeError(f"Reference N{count} failed: exit {code}")
    return record


def run_native_case(root, specification):
    kind, problem, count, directory = specification
    directory.mkdir(parents=True, exist_ok=False)
    if problem == "blast":
        (directory / "blast").mkdir()
        shutil.copy2(root / "native_references/blast/ref.dat", directory / "blast/ref.dat")
    if problem == "shuosher":
        (directory / "shuosher").mkdir()
        table = np.loadtxt(root / "N6400/ShuReference.block0.out1.00001.tab")
        values = table[:, [1, 2, 4, 3]]
        energy = values[:, 3] / (.4*values[:, 1])
        np.savetxt(directory / "shuosher/ref.dat", np.column_stack((values, energy)), fmt="%.16e")
    command = [str(root / f"native_{kind}_solver"), "unified_implicit", str(count), problem]
    write_json(directory / "command.json", command)
    start = time.monotonic()
    with (directory / "stdout.log").open("x") as out, (directory / "stderr.log").open("x") as err:
        process = subprocess.Popen(command, cwd=directory, stdout=out, stderr=err,
                                   env=dict(os.environ, OMP_NUM_THREADS="1", TMPDIR=str(root)))
        write_json(directory / "running.json", {"pid": process.pid, "command": command})
        code = process.wait()
    record = {"kind": kind, "problem": problem, "initial_N": count, "exit_code": code,
              "wall_seconds": time.monotonic()-start, "directory": str(directory)}
    write_json(directory / "process.json", record)
    text = (directory / "stdout.log").read_text()
    if code or "failed=0" not in text or "efloor=0 retry=0 failed=0" not in text:
        raise RuntimeError(f"Native test failed: {record}, {text}")
    return record


def run_native_suite(root):
    problems = ("contact", "moving_contact", "sod", "blast", "noh", "lax", "dblfan", "collision")
    specifications = [(kind, problem, 200, root / f"regression_{kind}" / problem)
                      for kind in ("parent", "corrected") for problem in problems]
    with ThreadPoolExecutor(max_workers=4) as executor:
        regression_processes = list(executor.map(lambda s: run_native_case(root, s), specifications))
    regression = []
    for problem in problems:
        parent = root / "regression_parent" / problem / f"unified_{problem}.dat"
        corrected = root / "regression_corrected" / problem / f"unified_{problem}.dat"
        identical = parent.read_bytes() == corrected.read_bytes()
        regression.append({"problem": problem, "byte_identical": identical,
                           "sha256": hashlib.sha256(corrected.read_bytes()).hexdigest()})
        assert identical, problem
    write_json(root / "unchanged_problem_regression.json",
               {"comparisons": regression, "processes": regression_processes})
    specifications = [("corrected", "shuosher", n, root / "native_corrected" / f"N{n}")
                      for n in (200, 400, 800, 1600, 3200)]
    with ThreadPoolExecutor(max_workers=4) as executor:
        processes = list(executor.map(lambda s: run_native_case(root, s), specifications))
    write_json(root / "native_processes.json", processes)


def reference(root, count):
    case = root / f"N{count}"
    record = json.loads((case / "process.json").read_text())
    assert record["exit_code"] == 0
    files = sorted(case.glob("ShuReference.block0.out1.*.tab"))
    assert len(files) == 2, files
    first, final = np.loadtxt(files[0]), np.loadtxt(files[-1])
    header = files[-1].read_text().splitlines()[0]
    found = re.search(r"time=\s*([^\s]+)", header)
    assert found and abs(float(found.group(1)) - 1.8) < 1e-12, header
    # Athena formatted-table 1D columns: i, x1v, rho, pressure, vx, vy, vz.
    assert first.shape == final.shape == (count, 7), final.shape
    x = final[:, 1]
    h = 10. / count
    assert np.max(abs(x - (-5 + (np.arange(count) + .5)*h))) < 2e-14
    values = final[:, [2, 4, 3]]  # rho, vx, P, matching the native comparison.
    assert np.isfinite(values).all() and values[:, [0, 2]].min() > 0
    initial_density = np.where(x < -4, 3.857143,
                               1 + .2*(np.cos(5*np.pi*(x-h/2)) -
                                        np.cos(5*np.pi*(x+h/2)))/(5*np.pi*h))
    assert np.max(abs(first[:, 2] - initial_density)) < 2e-12
    stdout = (case / "stdout.log").read_text()
    assert "Terminating on time limit" in stdout
    final_time = re.findall(r"^time=([^\s]+) cycle=", stdout, re.MULTILINE)
    assert len(final_time) == 1 and abs(float(final_time[0])-1.8) < 1e-14
    stderr = (case / "stderr.log").read_text()
    assert "REFERENCE_FLOOR_WOULD_ACTIVATE" not in stderr
    record.update({"final_time": float(final_time[0]), "final_density_min": float(values[:, 0].min()),
                   "final_pressure_min": float(values[:, 2].min()),
                   "upstream_pressure_max_error": float(np.max(abs(values[x < -4, 2] - 10.33333))),
                   "exact_density_projection_max_error": float(np.max(abs(first[:, 2] - initial_density))),
                   "floor_or_reset_count": 0, "final_table": str(files[-1])})
    edges = np.linspace(-5, 5, count+1)
    return edges[:-1], edges[1:], values, record


def overlap(left, right, values, other_left, other_right, other_values, lo=-4.9, hi=4.9):
    """Unnormalised spatial integral, using the actual cell intersections."""
    assert np.all(right > left) and np.all(other_right > other_left)
    assert left[0] <= lo and right[-1] >= hi
    assert other_left[0] <= lo and other_right[-1] >= hi
    i = j = 0
    norm = np.zeros(3)
    covered = 0.
    while i < len(left) and j < len(other_left):
        a, b = max(lo, left[i], other_left[j]), min(hi, right[i], other_right[j])
        if b > a:
            norm += (b-a)*abs(values[i]-other_values[j])
            covered += b-a
        if right[i] < other_right[j]:
            i += 1
        else:
            j += 1
    assert abs(covered - (hi-lo)) < 1e-10, covered
    return {"density_L1": float(norm[0]), "velocity_L1": float(norm[1]),
            "pressure_L1": float(norm[2]), "covered_length": covered}


def analyze(root, package):
    counts = tuple(json.loads((root / "registration.json").read_text())["N"])
    assert counts in ((1600, 3200, 6400), (1600, 3200, 6400, 12800))
    refs = {n: reference(root, n) for n in counts}
    convergence = []
    for coarse, fine in zip(counts[:-1], counts[1:]):
        convergence.append({"coarse": coarse, "fine": fine,
                            **overlap(*refs[coarse][:3], *refs[fine][:3])})
    comparisons = []
    # The earlier strain_resolution series uses sin(5*x). It must never be
    # compared to this sin(5*pi*x) reference or silently treated as the same
    # test. Only use new native runs with verified source/IC provenance.
    series = [("implicit_KDK_corrected_IC", n, root / "native_corrected" / f"N{n}")
              for n in (200, 400, 800, 1600, 3200)]
    manifest = json.loads((root / "native_manifest.json").read_text())
    assert manifest["density_wavenumber"] == "5*pi"
    for method, count, directory in series:
        native = directory / "unified_shuosher.dat"
        saved = np.loadtxt(native)
        # W0 Voronoi faces are neighboring generator midpoints. A generator
        # is generally not its cell centroid; x +/- stored V/2 is incorrect.
        assert np.all(np.diff(saved[:, 0]) > 0) and np.isfinite(saved).all()
        x = saved[:, 0]
        left, right = .5*(x[:-2]+x[1:-1]), .5*(x[1:-1]+x[2:])
        assert np.max(abs((right-left)-saved[1:-1, 6])) < 2e-9
        values = saved[1:-1, [1, 2, 3]]
        assert values[:, [0, 2]].min() > 0
        text = (directory / "stdout.log").read_text()
        terminal = re.search(r"shuosher N=\d+ t=([^\s]+) target=([^\s]+).* failed=(\d+)", text)
        process = json.loads((directory / "process.json").read_text())
        assert process["exit_code"] == 0 and terminal, text
        assert abs(float(terminal.group(1))-1.8) < 1e-14 and terminal.group(3) == "0", text
        assert "efloor=0 retry=0 failed=0" in text, text
        comparisons.append({"method": method, "initial_N": count, "saved_cells": len(saved),
                            "native_file": str(native), "native_sha256": hashlib.sha256(native.read_bytes()).hexdigest(),
                            **overlap(left, right, values, *refs[counts[-1]][:3])})
    result = {"scope": "Independent Euler reference comparison; no native checkpoints rewritten",
              "reference": "Athena++ HLLC/PPM/RK3 at t=1.8, same sin(5 pi x) initial condition",
              "reference_final_N": counts[-1],
              "norm": "Integral of absolute piecewise-constant cell-state differences over [-4.9,4.9]",
              "fields": "rho, mean velocity, mean-state pressure; no subcell kinetic recovery",
              "reference_processes": [refs[n][3] for n in counts],
              "reference_self_convergence": convergence, "saved_native_vs_reference": comparisons,
              "native_provenance": manifest,
              "limits": ["Reference is a refined numerical approximation with reported uncertainty, not an analytic solution.",
                         "Previous sin(5*x) common runs are different ICs, excluded from this comparison.",
                         "This does not qualify the new multidimensional method or matched GIZMO cost."]}
    write_json(root / "result.json", result)
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("package", type=Path)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--run-native", action="store_true")
    args = parser.parse_args()
    if args.run:
        with ThreadPoolExecutor(max_workers=3) as executor:
            results = list(executor.map(lambda n: run_case(args.root, n), (1600, 3200, 6400)))
        write_json(args.root / "processes.json", results)
    if args.run_native:
        run_native_suite(args.root)
    analyze(args.root, args.package)


if __name__ == "__main__":
    main()
