from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
data = json.loads((ROOT / "results.json").read_text())
original_path = ROOT / "runs/gfs_original_N200_t0.25_dt1/initial.dat"
x = np.loadtxt(original_path)
artifacts = []
archive = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/src/Hydro1DExam")
for name in ["gfs_contact.dat", "gfs_pair_contact.dat", "vr_contact.dat",
             "contact/voronoi.dat", "contact/hllc_cd10.dat", "contact/ns_av.dat",
             "contact/sph.dat"]:
    p = archive / name
    z = np.loadtxt(p)
    artifacts.append({
        "path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        "density_L1_site": float(np.sum(abs(z[:, 1] - np.where(z[:, 0] < .5, 1, .125))*z[:, 6])),
        "max_pressure_error": float(np.max(abs(z[:, 3] - 1))),
        "max_speed": float(np.max(abs(z[:, 2]))),
        "ordered_sites": bool(np.all(np.diff(z[:, 0]) > 0)),
        "sum_volumes": float(z[:, 6].sum())})

summary = {
    "contact_initial_cells": x[177:179].tolist(),
    "contact_initial_face": float(x[177, 9]),
    "initial_domain_right": float(x[-1, 9]),
    "resolution_orders": [], "saved_artifacts": artifacts, "causal_checks": {}}
for n in (198, 396):
    a = next(c for c in data["cases"] if c["method"] == "gfs" and
             c["scenario"] == "original" and c["N"] == n and c["t"] == .25)
    b = next(c for c in data["cases"] if c["method"] == "gfs" and
             c["scenario"] == "original" and c["N"] == 2*n and c["t"] == .25)
    summary["resolution_orders"].append(
        {"N_pair": [n, 2*n], "order": float(np.log2(a["final"][0]/b["final"][0]))})
c = summary["causal_checks"]
c["33_bounded_controls"] = len(data["cases"]) == 33
c["initial_pressure_defect_before_motion"] = bool(
    abs(x[177, 3]-.3789473684210515) < 1e-13 and
    abs(x[178, 3]-1.3031674208144834) < 1e-13)
controls = [a for a in data["cases"] if a["scenario"] != "original"]
c["all_consistent_controls_reach_t025"] = all(
    a["t"] == .25 and not a["failed"] for a in controls)
c["all_consistent_controls_stationary"] = all(a["final"][4] < 1e-10 for a in controls)
c["all_control_masses_fixed_during_evolution"] = all(
    a["max_mass_change"] == 0 for a in controls)
c["all_runs_no_thermal_floor"] = all(a["floor_hits"] == 0 for a in data["cases"])
c["saved_gfs_pair_result_reproduced"] = bool(
    abs(artifacts[1]["density_L1_site"]-next(
        a for a in data["cases"] if a["run"] ==
        "runs/gfs_original_N200_t0.25_dt1")["final"][0]) < 2e-12)
c["one_unrepaired_HLLC_failure_retained"] = sum(a["failed"] for a in data["cases"]) == 1
(ROOT / "causal_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
assert all(c.values()), c

plt.rcParams.update({"font.size": 11, "axes.labelsize": 11, "legend.fontsize": 11})
fig, axes = plt.subplots(2, 3, figsize=(11.6, 5.5), sharex=True, sharey="row",
                         gridspec_kw={"wspace": 0, "hspace": 0})
panels = [
    ("Original IC, t = 0", "gfs_original_N198_t0.25_dt1", "initial.dat", "#d95f02"),
    ("Original IC, t = 0.25", "gfs_original_N198_t0.25_dt1", "final.dat", "#d95f02"),
    ("Aligned Laguerre, t = 0.25", "gfs_laguerre_aligned_N198_t0.25_dt1", "final.dat", "#0072b2")]
for j, (label, run, file, color) in enumerate(panels):
    z = np.loadtxt(ROOT / "runs" / run / file)
    for row, col in enumerate((1, 3)):
        ax = axes[row, j]
        ax.plot([.46, .5, .5, .56],
                [1, 1, .125, .125] if row == 0 else [1, 1, 1, 1],
                color=".55", lw=1.5)
        # Draw cell averages as the actual piecewise-constant field.
        edges = np.r_[z[:, 8], z[-1, 9]]
        ax.stairs(z[:, col], edges, color=color, lw=1.2)
        ax.plot(z[:, 0], z[:, col], ".", color=color, ms=4)
        ax.set_xlim(.46, .56)
        ax.set_ylim(.05, 1.38)
        ax.tick_params(top=True, right=True, direction="in")
        ax.set_xticks([.47, .50, .53])
        ax.text(.05, .93 if row == 0 else .09,
                label if row == 0 else "N = 198, equal mass",
                transform=ax.transAxes, ha="left",
                va="top" if row == 0 else "bottom", fontsize=11)
        if row == 1:
            ax.set_xlabel("x")
axes[0, 0].set_ylabel(r"$\rho$")
axes[1, 0].set_ylabel("P")
fig.subplots_adjust(left=.065, right=.99, bottom=.1, top=.99)
fig.savefig(ROOT / "contact_causal_comparison.png", dpi=180)
fig.savefig(ROOT / "contact_causal_comparison.pdf")
plt.close(fig)
print(json.dumps(c, indent=2))
