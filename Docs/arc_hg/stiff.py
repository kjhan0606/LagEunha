"""Energy curvature of a prescribed azimuthal mode delta_theta = eps cos(m theta) applied to two adjacent rings,
for the disc-clipped Voronoi thermodynamics U = sum K m_i^g V_i^(1-g)/(g-1) (matched masses, uniform P).
d2U/deps2 > 0: restoring (stable); < 0: the reversible force drives the mode. Edge rings (0,1) vs interior (20,21).
Relative phase of the two rings is scanned (shear slides them)."""
import sys, numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce/EdgeModeStiff'); import arcgeom as g
from scipy.spatial import cKDTree
G = 5/3; d = 6/384; R = d; RIN = 0.5
def rings(kmax):
    P = []; lab = []
    for k in range(kmax):
        r = RIN + (k + 0.5)*d; n = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(n) + 0.5)*2*np.pi/n
        P.append(np.stack([r*np.cos(th), r*np.sin(th)], 1)); lab += [k]*n
    return np.concatenate(P), np.array(lab)
def energy(X, cells, m_i):
    tree = cKDTree(X); U = 0.0
    for a in cells:
        nb = [q for q in tree.query_ball_point(X[a], 2*R + 1e-12) if q != a]
        V, _ = g.area_moment(g.clipped_cell(X[nb] - X[a], R), R)
        U += m_i[a]*(m_i[a]/V)**(G - 1)/(G - 1)        # K = 1, rho0 = 1 -> P0 = 1
    return U
X0, lab = rings(26)
tree = cKDTree(X0)
m0 = np.zeros(len(X0))
for a in range(len(X0)):
    nb = [q for q in tree.query_ball_point(X0[a], 2*R + 1e-12) if q != a]
    m0[a] = g.area_moment(g.clipped_cell(X0[nb] - X0[a], R), R)[0]     # matched masses (rho = 1)
th0 = np.arctan2(X0[:, 1], X0[:, 0]); r0 = np.hypot(*X0.T)
for (ka, kb) in ((0, 1), (20, 21)):
    cells = np.where((lab >= max(ka - 2, 0)) & (lab <= kb + 2))[0]
    for mm in (45, 13):
        res = []
        for phase in (0.0, 0.25, 0.5, 0.75):
            X = X0.copy(); sel_b = lab == kb
            # relative phase: rotate ring kb by phase*spacing
            nb_ = int(round(2*np.pi*(RIN + (kb + 0.5)*d)/d)); dphi = phase*2*np.pi/nb_
            tb = th0[sel_b] + dphi; X[sel_b] = np.stack([r0[sel_b]*np.cos(tb), r0[sel_b]*np.sin(tb)], 1)
            thc = np.arctan2(X[:, 1], X[:, 0]); eps0 = 0.01*2*np.pi/len(np.where(lab == ka)[0])
            Us = []
            for e in (-eps0, 0.0, eps0):
                Y = X.copy(); s = (lab == ka) | (lab == kb); t = thc[s] + e*np.cos(mm*thc[s]); rr = np.hypot(*X[s].T)
                Y[s] = np.stack([rr*np.cos(t), rr*np.sin(t)], 1); Us.append(energy(Y, cells, m0))
            curv = (Us[0] - 2*Us[1] + Us[2])/eps0**2
            res.append('phase %.2f: %+.3e' % (phase, curv))
        print('rings %d,%d  m=%d  d2U/deps2: ' % (ka, kb, mm) + ' | '.join(res), flush=True)
