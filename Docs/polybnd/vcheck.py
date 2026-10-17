import sys, glob, numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce/PolyBnd'); sys.path.insert(0, '/tmp/lageunha-push/Docs/diag9')
import arcgeom as g, lagsnap
from scipy.spatial import cKDTree
n = 192; d = 6/n; R = d
fs = sorted(glob.glob('/gpfs/kjhan/LagForce/Kepler_n192_ringIC_cv_bfd170d/khout.*.dat'))
for f in fs:
    s = lagsnap.read(f); t = s['t']
    if round(t) % 6 or t < 5 or t > 44: continue
    X = np.stack([s['x'] - 3, s['y'] - 3], 1); tree = cKDTree(X); Vc = s['volume']
    Vp = np.zeros(len(X))
    for a in range(len(X)):
        nb = [q for q in tree.query_ball_point(X[a], 2*R + 1e-9) if q != a]
        Vp[a] = g.area_moment(g.clipped_cell(X[nb] - X[a], R), R)[0]
    e = np.abs(Vc - Vp)/Vp; r = np.hypot(*X.T); bad = np.where(e > 1e-8)[0]
    print('t=%5.1f  max rel |V_C - V_exact| %.2e  cells > 1e-8: %d  (radii of worst: %s)' % (t, e.max(), len(bad), np.round(np.sort(r[bad])[:5], 3)), flush=True)
