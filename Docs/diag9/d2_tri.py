"""D2 for candidate R1: material-triangle areas (fixed initial connectivity) along the prescribed Gresho trajectory."""
import numpy as np
from scipy.spatial import Delaunay
def om(r): return np.where(r < 0.2, 5.0, np.where(r < 0.4, 2/r - 5, 0.0))
for n in (32, 64, 128):
    h = 1/n; i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
    x0 = np.stack([(i.ravel()+.5)*h, (j.ravel()+.5)*h], 1) + np.random.default_rng(3).uniform(-1e-6, 1e-6, (n*n, 2))*h
    tri = Delaunay(x0).simplices; X0 = x0 - 0.5; r = np.hypot(*X0.T); th0 = np.arctan2(X0[:, 1], X0[:, 0])
    def area(x): a, b, c = x[tri[:, 0]], x[tri[:, 1]], x[tri[:, 2]]; return 0.5*((b-a)[:, 0]*(c-a)[:, 1] - (b-a)[:, 1]*(c-a)[:, 0])
    A0 = area(x0); rc = r[tri].mean(1); msk = rc < 0.45
    out = []
    for t in (0.25, 0.5, 1.0):
        th = th0 + om(r)*t; x = 0.5 + np.stack([r*np.cos(th), r*np.sin(th)], 1); A = area(x)
        out.append('t=%.2f rms dA/A %.4f max %.3f inverted %d' % (t, np.sqrt(np.mean((A/A0 - 1)[msk]**2)), np.abs(A/A0 - 1)[msk].max(), int((A[msk]*A0[msk] <= 0).sum())))
    print('n=%d  ' % n + ' | '.join(out))
