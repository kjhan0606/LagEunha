"""Feasibility of candidate G1 on the prescribed Gresho trajectory: per-particle metric G_i = F_i^-T F_i^-1
(F_i exact deformation gradient), cells C_i = argmin_j (r-x_j)^T G_j (r-x_j), volumes by rasterization;
compared with plain Voronoi on the same positions. Material volume = initial volume (incompressible flow)."""
import sys, numpy as np
from scipy.spatial import cKDTree
def om(r): return np.where(r < 0.2, 5.0, np.where(r < 0.4, 2/r - 5, 0.0))
def omp(r): return np.where(r < 0.2, 0.0, np.where(r < 0.4, -2/r**2, 0.0))
n = int(sys.argv[1]); M = int(sys.argv[2]); K = 40
h = 1/n; i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
x0 = np.stack([(i.ravel()+.5)*h, (j.ravel()+.5)*h], 1); X0 = x0 - 0.5; r = np.hypot(*X0.T); th0 = np.arctan2(X0[:, 1], X0[:, 0])
g = (np.arange(M) + 0.5)/M; GX, GY = np.meshgrid(g, g, indexing='ij'); S = np.stack([GX.ravel(), GY.ravel()], 1); pa = 1/M**2
def volumes(x, Ginv_list):
    xs = np.concatenate([x + np.array(s) for s in [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)]]); ids = np.tile(np.arange(len(x)), 9)
    Gs = np.tile(Ginv_list, (9, 1, 1)); tree = cKDTree(xs); V = np.zeros(len(x)); far = 0
    for c0 in range(0, len(S), 400000):
        P = S[c0:c0+400000]; d, k = tree.query(P, K)
        D = P[:, None, :] - xs[k]; phi = np.einsum('pki,pkij,pkj->pk', D, Gs[k], D)
        a = phi.argmin(1); far += int((a == K - 1).sum())
        np.add.at(V, ids[k[np.arange(len(P)), a]], pa)
    return V, far
V0 = np.full(n*n, h*h); msk = r < 0.45
for t in (0.25, 0.5, 1.0):
    th = th0 + om(r)*t; x = 0.5 + np.stack([r*np.cos(th), r*np.sin(th)], 1)
    c, s = np.cos(om(r)*t), np.sin(om(r)*t); Rm = np.stack([np.stack([c, -s], -1), np.stack([s, c], -1)], -2)
    that = np.stack([-np.sin(th), np.cos(th)], 1); r0 = X0/np.maximum(r, 1e-12)[:, None]
    F = Rm + (r*omp(r)*t)[:, None, None]*that[:, :, None]*r0[:, None, :]
    Fi = np.linalg.inv(F); Gm = np.einsum('nji,njk->nik', Fi, Fi)
    Vg, fg = volumes(x, Gm); Vv, fv = volumes(x, np.tile(np.eye(2), (n*n, 1, 1)))
    print('n=%d raster %d t=%.2f  Voronoi rms dV/V %.4f | G1 metric rms dV/V %.4f max %.3f (cond F max %.1f, far-assign %d/%d)' %
          (n, M, t, np.sqrt(np.mean((Vv/V0 - 1)[msk]**2)), np.sqrt(np.mean((Vg/V0 - 1)[msk]**2)), np.abs(Vg/V0 - 1)[msk].max(), np.linalg.cond(F).max(), fg, fv), flush=True)
