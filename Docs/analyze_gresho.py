"""L1 error of v_phi for the Gresho vortex in khout.NNNNNN.dat snapshots (box [0,1]^2, centre 0.5).

Usage: python analyze_gresho.py RUNDIR [RUNDIR ...]
"""
import glob
import os
import struct
import sys

import numpy as np


def load(path):
    with open(path, 'rb') as f:
        first = struct.unpack('i', f.read(4))[0]
        if first == -1:
            npart, psize, _stress = struct.unpack('iii', f.read(12))
        else:
            npart, psize = first, 360
        data = f.read(npart * psize)
        trailer = f.read()
    a = np.frombuffer(data, dtype=np.uint8).reshape(npart, psize)
    col = lambda o: a[:, o:o + 8].copy().view(np.float64).ravel()
    t = struct.unpack('d', trailer[0:8])[0]
    return t, col(24), col(32), col(56), col(64)


def vphi_exact(r):
    return np.where(r < 0.2, 5 * r, np.where(r < 0.4, 2 - 5 * r, 0.0))


for run in sys.argv[1:]:
    for p in sorted(glob.glob(os.path.join(run, 'khout.*.dat'))):
        t, x, y, vx, vy = load(p)
        rx, ry = x - 0.5, y - 0.5
        r = np.hypot(rx, ry) + 1e-30
        vphi = (-ry * vx + rx * vy) / r
        sel = r < 0.5
        L1 = np.mean(np.abs(vphi[sel] - vphi_exact(r[sel])))
        print(f"{os.path.basename(run):32s} t={t:6.3f} N={len(x)} L1(vphi, r<0.5)={L1:.4e}")
