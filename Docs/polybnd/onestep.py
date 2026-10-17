"""C vs python one-step comparison from the same snapshot (n192, reversible: Q=0, HG off, fixed dt).
Python: exact disc-clipped Voronoi (arcgeom), F_k = sum_a P_a dV_a/dx_k (faces A/d (c-x), arcs rigid),
P_a = K_a (m_a/V_a)^g with K_a from the snapshot (fixed entropy), central mass GM=1 eps=0.01, RK4."""
import sys, numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce/PolyBnd'); import arcgeom as g, lagsnap
from scipy.spatial import cKDTree
from multiprocessing import Pool
G_ = 5/3; EPS = 0.01; n = 192; d = 6/n; Rc = d; C0 = 3.0
B = '/gpfs/kjhan/LagForce/KepOneStep_bfd170d/'
def cellgrad(args):
    idx, X, nbl = args; out = []
    for a, nb in zip(idx, nbl):
        D = X[nb] - X[a]; pcs = g.clipped_cell(D, Rc); V, M = g.area_moment(pcs, Rc); gr = {a: np.zeros(2)}
        for pc in pcs:
            if pc[0] == 'seg':
                p, q, L = pc[1], pc[2], pc[5]
                if isinstance(L, str): raise RuntimeError('box edge')
                A = np.linalg.norm(q - p); c = 0.5*(p + q); dj = np.linalg.norm(D[L]); j = nb[L]
                gr[a] += A/dj*c; gr[j] = gr.get(j, np.zeros(2)) - A/dj*(c - D[L])
            else:
                th0, th1 = pc[1], pc[2]; gr[a] += Rc*np.array([np.sin(th1) - np.sin(th0), np.cos(th0) - np.cos(th1)])
        out.append((a, V, gr))
    return out
pool = None
def acc(X, m, K):
    tree = cKDTree(X); N = len(X)
    nbl = [[q for q in tree.query_ball_point(X[a], 2*Rc + 1e-12) if q != a] for a in range(N)]
    chunks = np.array_split(np.arange(N), 64)
    res = pool.map(cellgrad, [(ch, X, [nbl[a] for a in ch]) for ch in chunks])
    V = np.zeros(N); F = np.zeros_like(X); grs = [None]*N
    for r in res:
        for a, v, gr in r: V[a] = v; grs[a] = gr
    P = K*(m/V)**G_
    for a in range(N):
        for k, v in grs[a].items(): F[k] += P[a]*v
    r2 = (X*X).sum(1) + EPS*EPS
    return F/m[:, None], -X/r2[:, None]**1.5, V
if __name__ == '__main__':
    pool = Pool(16)
    s0 = lagsnap.read(B + 'khout.000541.dat'); o = np.argsort(s0['id'])
    X = np.stack([s0['x'][o] - C0, s0['y'][o] - C0], 1); Vv = np.stack([s0['vx'][o], s0['vy'][o]], 1); m = s0['mass'][o]
    Pc = (G_ - 1)*s0['ie'][o]/s0['volume'][o]; K = Pc/(m/s0['volume'][o])**G_
    R0 = np.hypot(*X.T); k = np.floor((R0 - 0.5)/d).astype(int)
    dt = 0.02
    ah, ag, V = acc(X, m, K)
    print('t0 %.6f  max rel |V_py - V_C| %.2e' % (s0['t'], np.max(np.abs(V - s0['volume'][o])/V)), flush=True)
    for step in range(1, 4):
        def tot(X_): a_, g_, _ = acc(X_, m, K); return a_ + g_
        k1x = Vv; k1v = tot(X); k2x = Vv + 0.5*dt*k1v; k2v = tot(X + 0.5*dt*k1x)
        k3x = Vv + 0.5*dt*k2v; k3v = tot(X + 0.5*dt*k2x); k4x = Vv + dt*k3v; k4v = tot(X + dt*k3x)
        X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        sc = lagsnap.read(B + 'khout.%06d.dat' % (541 + step)); oc = np.argsort(sc['id'])
        Xc = np.stack([sc['x'][oc] - C0, sc['y'][oc] - C0], 1); Vc = np.stack([sc['vx'][oc], sc['vy'][oc]], 1)
        ah, _, Vp = acc(X, m, K)
        out = []
        for ring in (0, 1, 2, 10, 40):
            sel = k == ring
            dv = np.hypot(*(Vc[sel] - Vv[sel]).T); dx = np.hypot(*(Xc[sel] - X[sel]).T); ahn = np.hypot(*ah[sel].T)
            out.append('k%d |dv| %.1e (|a_hyd|dt %.1e) |dx| %.1e' % (ring, dv.max(), ahn.max()*dt, dx.max()))
        Pcs = (G_ - 1)*sc['ie'][oc]/sc['volume'][oc]; Kc = Pcs/(m/sc['volume'][oc])**G_
        print('step %d t %.4f | %s | entropy C drift max %.1e' % (step, sc['t'], ' | '.join(out), np.max(np.abs(Kc/K - 1))), flush=True)
