"""One-body RK4 control for the native softened Kepler orbital dt limiter.

The force and initial circular speed match the native EUNHA_IC=kepler path.
This excludes hydrodynamics and mesh motion; it only tests whether increasing
the *global* RK4 step to the observed hydro CFL resolves inner orbits.
"""

import argparse
import json
import math
import struct
from pathlib import Path

import numpy as np


def first_nonzero_radius(glass: Path) -> float:
    with glass.open("rb") as stream:
        n, lx, ly = struct.unpack("<idd", stream.read(20))
        if (n, lx, ly) != (384*384, 6.0, 6.0):
            raise ValueError(f"Unexpected glass header: {glass}")
        x = np.fromfile(stream, dtype="<f8", count=n)
        y = np.fromfile(stream, dtype="<f8", count=n)
    r = np.hypot(x-3.0, y-3.0)
    return float(np.min(r[r > 1e-8]))


def rhs(state: tuple[float, float, float, float], eps: float):
    x, y, vx, vy = state
    inv_r3 = (x*x+y*y+eps*eps)**-1.5
    return (vx, vy, -x*inv_r3, -y*inv_r3)


def rk4(state: tuple[float, float, float, float], dt: float, eps: float):
    k1 = rhs(state, eps)
    k2 = rhs(tuple(s+.5*dt*k for s, k in zip(state, k1)), eps)
    k3 = rhs(tuple(s+.5*dt*k for s, k in zip(state, k2)), eps)
    k4 = rhs(tuple(s+dt*k for s, k in zip(state, k3)), eps)
    return tuple(s+dt*(a+2*b+2*c+d)/6
                 for s, a, b, c, d in zip(state, k1, k2, k3, k4))


def integrate(radius: float, eps: float, dt: float, end: float,
              subcycle_limit: float | None = None) -> dict:
    omega = (radius*radius+eps*eps)**-.75
    initial = (radius, 0.0, 0.0, radius*omega)
    state = initial
    t = 0.0
    macro_count = 0
    count = 0
    while t < end-1e-14:
        h = min(dt, end-t)
        nsub = 1 if subcycle_limit is None else int(math.ceil(h/subcycle_limit))
        for _ in range(nsub):
            state = rk4(state, h/nsub, eps)
            count += 1
        t += h
        macro_count += 1
        if not all(math.isfinite(v) for v in state):
            raise ValueError(f"Nonfinite state after {count} steps")
    x, y, vx, vy = state
    exact_x = radius*math.cos(omega*end)
    exact_y = radius*math.sin(omega*end)
    energy0 = .5*(radius*omega)**2-1/math.sqrt(radius*radius+eps*eps)
    energy = .5*(vx*vx+vy*vy)-1/math.sqrt(x*x+y*y+eps*eps)
    L0 = radius*radius*omega
    L = x*vy-y*vx
    return {"radius": radius, "omega_dt": omega*dt, "requested_dt": dt,
            "subcycle_limit": subcycle_limit, "endpoint": end,
            "macro_steps": macro_count, "gravity_substeps": count,
            "position_error_over_initial_radius":
            math.hypot(x-exact_x, y-exact_y)/radius,
            "radius_relative_error": math.hypot(x, y)/radius-1,
            "angular_momentum_relative_error": L/L0-1,
            "energy_error_over_initial_binding_energy": (energy-energy0)/abs(energy0)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("q1_glass", type=Path)
    ap.add_argument("q2_glass", type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    eps = .01175
    orbital_global = .06*eps**1.5
    cases = []
    for name, glass, hydro_dt in (("q1", args.q1_glass, .000679791730362922),
                                   ("q2", args.q2_glass, .002267122268676758)):
        inner = first_nonzero_radius(glass)
        for radius_name, radius in (("innermost_nonzero", inner), ("disk_inner_edge", .5)):
            for dt_name, dt, subcycle_limit in (
                    ("native_orbital_bound", orbital_global, None),
                    ("observed_global_hydro_cfl", hydro_dt, None),
                    ("hydro_macro_with_local_orbital_gravity_subcycle", hydro_dt,
                     .06*(radius*radius+eps*eps)**.75),
                    ("hydro_macro_with_global_orbital_gravity_subcycle", hydro_dt,
                     orbital_global)):
                for end in (.825, 6.0):
                    result = integrate(radius, eps, dt, end, subcycle_limit)
                    result.update(layout=name, radius_case=radius_name,
                                  dt_case=dt_name, glass=str(glass))
                    cases.append(result)
    output = {"scope": "Pure softened central-force RK4 control, exact circular "
              "orbit reference. No mesh or hydro forces. Native double-precision "
              "position/velocity arithmetic modeled; not a full solver result.",
              "eps": eps, "global_orbital_bound": orbital_global, "cases": cases}
    args.output.write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
