"""Compile and run bounded causal controls using unchanged archived C code.

Existing accepted results are reused. Each new calculation gets a new directory.
"""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import time

ROOT = Path(__file__).resolve().parent
COLUMNS = ["density_L1_site", "density_L1_cell_average",
           "density_L1_exact_piecewise", "max_pressure_error", "max_speed",
           "max_acceleration", "energy", "mass", "domain_left", "domain_right"]
started = time.time()
compile_run = subprocess.run(
    ["gcc", "-O2", "-Wall", "-Wextra", "-o", str(ROOT / "contact_probe"),
     str(ROOT / "contact_probe.c"), "-lm"], capture_output=True, text=True)
(ROOT / "build.stdout").write_text(compile_run.stdout)
(ROOT / "build.stderr").write_text(compile_run.stderr)
compile_run.check_returncode()
cases = []
for n in (198, 200, 396, 400, 792, 800):
    cases.append(("gfs", "original", n, .25, 1))
for n in (198, 200, 396, 400):
    for scenario in ("equalvolume", "pressure_only", "fv_projected"):
        cases.append(("gfs", scenario, n, .25, 1))
for n in (198, 396):
    cases.append(("gfs", "laguerre_aligned", n, .25, 1))
for method in ("hllc", "ns"):
    for scenario in ("original", "equalvolume", "pressure_only", "fv_projected"):
        cases.append((method, scenario, 200, .25, 1))
for dt in (.5, .25):
    cases.append(("gfs", "original", 200, .25, dt))
for end in (.0001, .005, .05):
    cases.append(("gfs", "original", 200, end, 1))

env = dict(os.environ)
# Ensure controls are independent of switches in an interactive shell.
for key in list(env):
    if key.startswith("GFS_"):
        del env[key]
rows = []
for method, scenario, n, end, dt in cases:
    name = f"{method}_{scenario}_N{n}_t{end}_dt{dt}"
    run = ROOT / "runs" / name
    if (run / "stdout.json").is_file() and (run / "stdout.json").stat().st_size:
        data = json.loads((run / "stdout.json").read_text())
    else:
        # The pre-run wrapper exited without output at the first HLLC failure.
        # Preserve the empty attempt and retry in a distinct directory.
        if run.exists():
            run = ROOT / "runs" / f"{name}_failure_record_v2"
        run.mkdir(parents=True, exist_ok=False)
        argv = [str(ROOT / "contact_probe"), method, scenario, str(n), str(end), str(dt)]
        result = subprocess.run(argv, cwd=run, env=env, capture_output=True, text=True)
        (run / "stdout.json").write_text(result.stdout)
        (run / "stderr.log").write_text(result.stderr)
        result.check_returncode()
        data = json.loads(result.stdout)
    data.setdefault("failed", 0)
    data["run"] = str(run.relative_to(ROOT))
    rows.append(data)
    print(f"{name}: t={data['t']:.7g} fail={data['failed']} "
          f"L1={data['final'][0]:.8g} vmax={data['final'][4]:.5g}", flush=True)
    (ROOT / "results.json").write_text(json.dumps(
        {"columns": COLUMNS, "cases": rows,
         "wall_seconds_this_invocation": time.time() - started}, indent=2) + "\n")

inputs = ["contact_probe.c", "run_audit.py", "source/Hydro1DExam/laguerre_sod.c",
          "source/Exam/gfs_pair.h", "contact_probe", "results.json"]
(ROOT / "audit_manifest.json").write_text(json.dumps({
    name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    for name in inputs}, indent=2) + "\n")
