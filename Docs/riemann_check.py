"""Compare the code's face pressure P* with (i) the exact ideal-gas Riemann
solution and (ii) the same Davis-HLLC with face-state sound speeds and no
early vacuum branch, on the faces recorded by GFS_EDGE_DUMP (av_mode 5,
reconstructed states). Torques are re-summed with the face lever arm
recovered from the record. Usage: python riemann_check.py RUNDIR [gamma]
"""
import os
import sys

import numpy as np
from scipy.optimize import brentq

RUN = sys.argv[1]
G = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0 / 3.0
cols = ("t call id_i r_i id_j r_j A d toff du_n_cell du_n_rec cs_i P_i Pstar extreme tq_n tq_t "
        "rhoL rhoR pL pR uL uR ua_n flux_E force_x vx_i").split()
a = np.loadtxt(os.path.join(RUN, 'edge_faces.txt'))
c = {k: i for i, k in enumerate(cols)}


def f_side(p, r, pk, ck):
    if p > pk:
        A, B = 2 / ((G + 1) * r), (G - 1) / (G + 1) * pk
        return (p - pk) * np.sqrt(A / (p + B))
    return 2 * ck / (G - 1) * ((p / pk) ** ((G - 1) / (2 * G)) - 1)


def exact_pstar(rl, pl, ul, rr, pr, ur):
    cl, cr = np.sqrt(G * pl / rl), np.sqrt(G * pr / rr)
    if ur - ul >= 2 * (cl + cr) / (G - 1):
        return 0.0
    fun = lambda p: f_side(p, rl, pl, cl) + f_side(p, rr, pr, cr) + ur - ul
    hi = max(pl, pr, 1e-30)
    while fun(hi) < 0: hi *= 10
    return brentq(fun, 0.0, hi, xtol=1e-14 * hi, rtol=1e-12)


def hllc_davis(rl, pl, ul, rr, pr, ur, cl, cr):
    SL, SR = min(ul - cl, ur - cr), max(ul + cl, ur + cr)
    den = rl * (SL - ul) - rr * (SR - ur)
    SM = (pr - pl + rl * ul * (SL - ul) - rr * ur * (SR - ur)) / den
    return pl + rl * (SL - ul) * (SM - ul)


for w0 in sorted(set(np.round(a[:, 0]))):
    sel = np.abs(a[:, 0] - w0) < 0.5
    b = a[sel]
    calls = len(set(b[:, c['call']]))
    rl, rr, pl, pr, ul, ur = (b[:, c[k]] for k in ('rhoL', 'rhoR', 'pL', 'pR', 'uL', 'uR'))
    Pc = b[:, c['Pstar']]
    arm = np.where(np.abs(Pc) > 0, -b[:, c['tq_n']] / np.where(Pc != 0, Pc, 1), 0.0)
    Pex = np.array([exact_pstar(*z) for z in zip(rl, pl, ul, rr, pr, ur)])
    clr, crr = np.sqrt(G * pl / rl), np.sqrt(G * pr / rr)
    Phl = np.array([max(hllc_davis(*z), 0.0) for z in zip(rl, pl, ul, rr, pr, ur, clr, crr)])
    du = ur - ul
    vac_code = Pc < 1e-20
    vac_exact = du >= 2 * (clr + crr) / (G - 1)
    print(f"== window t~{w0:.0f}: faces {len(b)} in {calls} force calls (r_i < GFS_EDGE_DUMP_R)")
    print(f"   code vacuum branch (P*~0): {vac_code.mean():.3f} of faces; exact vacuum: {vac_exact.mean():.3f}; "
          f"faces with c < du < 6c-type exact threshold: {(vac_code & ~vac_exact).mean():.3f}")
    for lo, hi in ((0, 0.49), (0.49, 0.53), (0.53, 0.6)):
        s = (b[:, c['r_i']] >= lo) & (b[:, c['r_i']] < hi)
        if not s.any(): continue
        tq = lambda P: -(P[s] * arm[s]).sum() / calls
        print(f"   r in [{lo:.2f},{hi:.2f}) n={s.sum():5d}: torque/call  code={tq(Pc):+.3e}  exact={tq(Pex):+.3e}  "
              f"HLLC(face c, no early vacuum)={tq(Phl):+.3e}   median P*/P_i code={np.median(Pc[s] / b[s, c['P_i']]):.3f} "
              f"exact={np.median(Pex[s] / b[s, c['P_i']]):.3f}")
    comp = du < 0
    print(f"   compression faces: <P*_code/P*_exact> = {np.median(Pc[comp] / np.maximum(Pex[comp], 1e-300)):.3f}; "
          f"expansion faces: <P*_code - P*_exact>/P_i = {np.median((Pc[~comp] - Pex[~comp]) / b[~comp, c['P_i']]):+.3e}")
