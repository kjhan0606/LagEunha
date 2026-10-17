"""Summarise edge_faces.txt (GFS_EDGE_DUMP): per time window and radial band of
the owner particle, the reconstructed normal jump over sound speed, P*/P,
and the torque split between the normal Riemann pressure and the (H) VPH
tangential term. Also groups faces by the tangential offset of the face
centre (long cap-clipped faces have large offsets).

Usage: python analyze_edge_faces.py RUNDIR [RUNDIR ...]
"""
import os
import sys

import numpy as np

for run in sys.argv[1:]:
    a = np.loadtxt(os.path.join(run, 'edge_faces.txt'))
    t, ri, A, d, toff = a[:, 0], a[:, 3], a[:, 6], a[:, 7], a[:, 8]
    duc, dur, cs, P, Ps, ext, tqn, tqt = (a[:, k] for k in (9, 10, 11, 12, 13, 14, 15, 16))
    print(f"== {os.path.basename(run)}  faces={len(t)}")
    for w0 in sorted(set(np.round(t, 1))):
        wsel = np.abs(t - w0) < 0.1
        calls = len(set(a[wsel, 1]))
        print(f"  window t~{w0:.1f}  (force calls {calls})")
        print("    band            nface  med|du_rec|/cs  med|du_cell|/cs  rec/cell  med P*/P  frac P*>10P  ext  "
              "tq_normal/call  tq_tang/call")
        for lo, hi in ((0.0, 0.49), (0.49, 0.53), (0.53, 0.56)):
            s = wsel & (ri >= lo) & (ri < hi)
            if s.sum() == 0:
                continue
            print(f"    [{lo:4.2f},{hi:4.2f})  {s.sum():6d}   {np.median(np.abs(dur[s]) / cs[s]):9.3e}     "
                  f"{np.median(np.abs(duc[s]) / cs[s]):9.3e}   {np.median(np.abs(dur[s]) / (np.abs(duc[s]) + 1e-30)):6.3f}"
                  f"  {np.median(Ps[s] / P[s]):8.3f}   {np.mean(Ps[s] > 10 * P[s]):6.3f}   {int(ext[s].sum()):4d}"
                  f"   {tqn[s].sum() / calls:+.3e}   {tqt[s].sum() / calls:+.3e}")
        s = wsel & (ri < 0.53)
        if s.sum():
            q = toff[s] / d[s]
            for lo, hi in ((0, 0.25), (0.25, 0.75), (0.75, 99)):
                k = (q >= lo) & (q < hi)
                if k.sum():
                    print(f"      r<0.53, toff/d in [{lo},{hi}): n={k.sum():5d} med|du_rec|/cs={np.median(np.abs(dur[s][k]) / cs[s][k]):.3e}"
                          f"  med P*/P={np.median(Ps[s][k] / P[s][k]):.3f}  tq_normal/call={tqn[s][k].sum() / calls:+.3e}")
