"""TCF reversible force (C dump, Q=0, HG off) vs -grad U of the disc-clipped Voronoi thermodynamics (python FD)."""
import sys, numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce/EdgeModeStiff'); import arcgeom as g
from scipy.spatial import cKDTree
G = 5/3; d = 6/384; R = d
c = np.loadtxt('/gpfs/kjhan/LagForce/KepR24_revdump_bfd170d/rhs', comments='#')
X = np.stack([c[:, 1] - 3, c[:, 2] - 3], 1); P = c[:, 5]; V = c[:, 6]; m = c[:, 7]; F = c[:, 8:10]
K = P/(m/V)**G; r = np.hypot(*X.T); k = np.floor((r - 0.5)/d).astype(int)
tree = cKDTree(X)
def vol(a, Y):
    nb = [q for q in tree.query_ball_point(X[a], 3*R) if q != a]
    return g.area_moment(g.clipped_cell(Y[nb] - Y[a], R), R)[0]
Vpy = np.array([vol(a, X) for a in range(len(X))]) if False else None
def Uloc(cells, Y):
    return sum(K[a]*m[a]**G*vol(a, Y)**(1 - G)/(G - 1) for a in cells)
rng = np.random.default_rng(0); h = 1e-7
for ring in (0, 1, 2, 20):
    idx = rng.choice(np.where(k == ring)[0], 12, replace=False); errs = []; rel = []
    for p in idx:
        cells = [q for q in tree.query_ball_point(X[p], 3*R)]
        Vc = vol(p, X); 
        fd = np.zeros(2)
        for comp in (0, 1):
            Yp = X.copy(); Yp[p, comp] += h; Ym = X.copy(); Ym[p, comp] -= h
            fd[comp] = -(Uloc(cells, Yp) - Uloc(cells, Ym))/(2*h)
        er = np.array([X[p, 0], X[p, 1]])/r[p]
        errs.append(np.hypot(*(F[p] - fd))); rel.append(np.hypot(*(F[p] - fd))/max(np.hypot(*fd), 1e-300))
        if p == idx[0]: print('  ring %d sample: V_C %.6e V_py %.6e | F_C (%.3e, %.3e) | -gradU (%.3e, %.3e)' % (ring, V[p], Vc, F[p, 0], F[p, 1], fd[0], fd[1]))
    print('ring %2d: |F_C + grad U| / |grad U| median %.2e max %.2e' % (ring, np.median(rel), np.max(rel)), flush=True)
