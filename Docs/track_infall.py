"""Track particles that fall toward the central mass in khout.NNNNNN.dat snapshots.

For every particle whose final radius is below RMIN, print r, specific
angular momentum j = (x v_y - y v_x) relative to the centre, and j / j_K(r_0)
at each snapshot, plus the snapshot interval where j dropped fastest.
Usage: python track_infall.py RUNDIR [cx cy RMIN]
"""
import glob
import os
import struct
import sys

import numpy as np

RUN = sys.argv[1]
CX = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
CY = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
RMIN = float(sys.argv[4]) if len(sys.argv) > 4 else 0.35
EPS = 0.01


def load(path):
    with open(path, 'rb') as f:
        first = struct.unpack('i', f.read(4))[0]
        npart, psize, _ = struct.unpack('iii', f.read(12))
        a = np.frombuffer(f.read(npart * psize), dtype=np.uint8).reshape(npart, psize)
        tr = f.read()
    d = lambda o: a[:, o:o + 8].copy().view(np.float64).ravel()
    idx = (a[:, 8:16].copy().view(np.uint64).ravel() << np.uint64(8)) >> np.uint64(8)
    return struct.unpack('d', tr[:8])[0], idx, d(24), d(32), d(56), d(64)


snaps = sorted(glob.glob(os.path.join(RUN, 'khout.*.dat')))
T, R, J = [], {}, {}
for p in snaps:
    t, idx, x, y, vx, vy = load(p)
    rx, ry = x - CX, y - CY
    r = np.hypot(rx, ry)
    j = rx * vy - ry * vx
    T.append(t)
    for k, i in enumerate(idx):
        R.setdefault(int(i), []).append(r[k])
        J.setdefault(int(i), []).append(j[k])
T = np.array(T)
fall = [i for i in R if len(R[i]) == len(T) and R[i][-1] < RMIN]
print(f"# {len(snaps)} snapshots t={T[0]:.2f}..{T[-1]:.2f}; particles with final r < {RMIN}: {len(fall)}")
for i in fall:
    r, j = np.array(R[i]), np.array(J[i])
    jk0 = r[0] * r[0] * (r[0] ** 2 + EPS ** 2) ** -0.75
    dj = np.diff(j) / np.diff(T)
    k = int(np.argmin(dj))
    print(f"\nparticle {i}: r0={r[0]:.4f} j0/jK={j[0]/jk0:.4f}  fastest loss dj/dt={dj[k]:.3e} in t=[{T[k]:.2f},{T[k+1]:.2f}]")
    step = max(1, len(T) // 17)
    for m in list(range(0, len(T), step)) + [len(T) - 1]:
        print(f"   t={T[m]:6.2f} r={r[m]:.4f} j/jK(r0)={j[m]/jk0:.4f}")
