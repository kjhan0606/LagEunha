"""Measure the published GIZMO Kepler-disk L1 diagnostics on a native snapshot.

Hopkins 2015, section 4.2.4, defines mean absolute velocity and density
errors in 0.8<R<1.2. The paper's density reference is 1, although its
plateau IC is 1.01; report both that literal metric and a 1.01-corrected
metric so the small initial offset is visible.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from audit_kepler_polar_checkpoint_20261026 import (
    CX, CY, EPS, load_ordered, weighted_mean,
)


def metrics(p, mask):
    x = p["x"][mask].astype(float) - CX
    y = p["y"][mask].astype(float) - CY
    r = np.hypot(x, y)
    vx = p["vx"][mask].astype(float)
    vy = p["vy"][mask].astype(float)
    vphi = (x*vy - y*vx) / r
    vk = r * (r*r + EPS*EPS)**(-.75)
    den = p["den"][mask].astype(float)
    weights = p["mass"][mask].astype(float)
    ev = np.abs(vphi-vk)/vk
    ed_paper = np.abs(den-1.0)
    ed_plateau = np.abs(den-1.01)/1.01
    return {
        "cells": int(np.sum(mask)),
        "mass": float(np.sum(weights)),
        "L1_velocity_paper": float(np.mean(ev)),
        "L1_density_paper": float(np.mean(ed_paper)),
        "L1_density_plateau_corrected": float(np.mean(ed_plateau)),
        "L1_velocity_mass_weighted": weighted_mean(ev, weights),
        "L1_density_paper_mass_weighted": weighted_mean(ed_paper, weights),
        "velocity_p95": float(np.quantile(ev, .95)),
        "density_plateau_corrected_p95": float(np.quantile(ed_plateau, .95)),
        "passes_paper_L1_velocity_0p01": bool(np.mean(ev) < .01),
        "passes_paper_L1_density_0p25": bool(np.mean(ed_paper) < .25),
    }


def audit(initial, endpoint):
    p0, meta0 = load_ordered(initial)
    pt, metat = load_ordered(endpoint)
    if float(meta0["t"]) != 0 or float(metat["t"]) < 0:
        raise ValueError("Expected t=0 initial and nonnegative endpoint")
    r0 = np.hypot(p0["x"].astype(float)-CX, p0["y"].astype(float)-CY)
    rt = np.hypot(pt["x"].astype(float)-CX, pt["y"].astype(float)-CY)
    cohort = (r0 > .8) & (r0 < 1.2)
    current = (rt > .8) & (rt < 1.2)
    if not np.any(current) or not np.any(cohort):
        raise ValueError("No cells in GIZMO reference annulus")
    if not (np.isfinite(pt["den"]).all() and np.isfinite(pt["vx"]).all() and
            np.isfinite(pt["vy"]).all()):
        raise ValueError("Nonfinite native snapshot")
    return {
        "initial": str(initial),
        "endpoint": str(endpoint),
        "time": float(metat["t"]),
        "paper_reference": "https://arxiv.org/html/1409.7395#S4.SS2.SSS4",
        "annulus": "0.8<R<1.2, current radius unless marked initial cohort",
        "paper_density_reference": 1.0,
        "actual_plateau_initial_density": 1.01,
        "initial_current_annulus": metrics(p0, cohort),
        "endpoint_current_annulus": metrics(pt, current),
        "endpoint_initial_cohort": metrics(pt, cohort),
        "scope": "Native fixed-mass Kepler W=0 snapshot; thresholds report paper Fig. 8 gates, not full t=600 stability",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.initial, args.endpoint)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
