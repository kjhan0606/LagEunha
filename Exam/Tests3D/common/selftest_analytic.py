#!/usr/bin/env python3
"""Checks of the analytic/reference solvers used by the 3D suite."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sedov_exact  # noqa: E402
import noh_exact  # noqa: E402
import kh_linear  # noqa: E402

ok = True


def check(name, val, ref, tol):
    global ok
    good = abs(val - ref) <= tol
    ok &= good
    print("%-58s %12.6f  ref %10.6f  tol %.1e  %s" % (name, val, ref, tol, "ok" if good else "FAIL"))


check("Sedov xi0, spherical, gamma=5/3", sedov_exact.sedov_xi0(5 / 3, 3), 1.15167, 5e-4)
check("Sedov xi0, spherical, gamma=1.4", sedov_exact.sedov_xi0(1.4, 3), 1.0328, 5e-4)
check("Sedov xi0, cylindrical, gamma=1.4 (alpha=0.984)", sedov_exact.sedov_xi0(1.4, 2), 0.984 ** -0.25, 1e-3)
rho, v, p, R = sedov_exact.sedov_profile(np.array([0.999999 * sedov_exact.shock_radius(0.05)]), 0.05)
check("Sedov post-shock density (gamma=5/3)", rho[0], 4.0, 1e-3)
rho, v, p, rs = noh_exact.noh_profile(np.array([0.2, 1.0]), 2.0)
check("Noh 3D post-shock density", rho[0], 64.0, 1e-12)
check("Noh 3D post-shock pressure", p[0], 64.0 / 3.0, 1e-12)
check("Noh 3D shock radius at t=2", rs, 2.0 / 3.0, 1e-12)
check("Noh 3D pre-shock density at r=1, t=2", rho[1], 9.0, 1e-12)
L = 0.01
s, w, _ = kh_linear.growth_rate(k=0.4446 / L, P0=1000.0, n=768,
                                profile=kh_linear.tanh_pair_profile(0.5, L))
check("KH solver vs Michalke 1964: sigma L/U0 at kL=0.4446", s * L / 0.5, 0.1897, 1e-3)
s, w, s0 = kh_linear.growth_rate(n=512)
s2, _, _ = kh_linear.growth_rate(n=256)
check("KH McNally setup sigma (n=512 vs n=256 grid convergence)", s, s2, 0.02)
print("McNally 2012 setup, compressible linear sigma = %.4f (sharp interface %.4f)" % (s, s0))
print("ANALYTIC", "OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
