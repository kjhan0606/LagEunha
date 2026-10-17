"""Parallel version of dyn.py (disc closure only), for N ~ 1e4: same equations (fixed entropy, RK4, GM=1, eps=0.01)."""
import sys, time, numpy as np
from multiprocessing import Pool
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
import shapely
G_ = 5/3; P0 = 1e-6; EPS = 0.01; RIN, ROUT = 0.5, 2.0
def J(v): return np.array([-v[1], v[0]])
def cell_job(args):
    idx, XA, NBs, Rc, nseg = args; out = []
    circ0 = np.array([[np.cos(th), np.sin(th)] for th in np.linspace(0, 2*np.pi, nseg, endpoint=False)])*Rc
    rmin = Rc*np.cos(np.pi/nseg) - 1e-9
    for a, xa, (nbid, XN) in zip(idx, XA, NBs):
        B = 8*Rc + 1.0
        cell = Polygon([(xa[0]-B, xa[1]-B), (xa[0]+B, xa[1]-B), (xa[0]+B, xa[1]+B), (xa[0]-B, xa[1]+B)])
        for xj in XN:
            dv = xj - xa; mid = 0.5*(xj + xa); t = J(dv)/np.linalg.norm(dv)*50; n = dv/np.linalg.norm(dv)*50
            cell = cell.intersection(Polygon([mid + t, mid - t, mid - t - n, mid + t - n]))
        cell = cell.intersection(Polygon(xa + circ0))
        V = cell.area; g = {}
        polys = [cell] if isinstance(cell, Polygon) else [q for q in cell.geoms if isinstance(q, Polygon)]
        for pgn in polys:
            cs = np.array(shapely.geometry.polygon.orient(pgn, 1.0).exterior.coords)
            for k in range(len(cs) - 1):
                p, q = cs[k], cs[k+1]; L = np.linalg.norm(q - p)
                if L < 1e-14: continue
                done = False
                for j, xj in zip(nbid, XN):
                    dv = xj - xa; dn = np.linalg.norm(dv); e = dv/dn; mid = 0.5*(xj + xa)
                    if abs((p - mid) @ e) < 1e-9 and abs((q - mid) @ e) < 1e-9:
                        c = 0.5*(p + q); g[a] = g.get(a, 0) + L/dn*(c - xa); g[j] = g.get(j, 0) - L/dn*(c - xj); done = True; break
                if not done:
                    if np.linalg.norm(p - xa) >= rmin and np.linalg.norm(q - xa) >= rmin and np.linalg.norm(0.5*(p + q) - xa) >= rmin:
                        g[a] = g.get(a, 0) - J(q - p)
                    else: raise RuntimeError('unclassified segment, cell %d' % a)
        out.append((a, V, g))
    return out
if __name__ == '__main__':
    n = int(sys.argv[1]); dt = float(sys.argv[2]); T = float(sys.argv[3]); nproc = int(sys.argv[4]); nseg = int(sys.argv[5])
    d = 6.0/n; Rc = d
    P = []; lab = []
    for k in range(int(round((ROUT - RIN)/d))):
        r = RIN + (k + 0.5)*d; m_ = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(m_) + 0.5)*2*np.pi/m_
        P.append(np.stack([r*np.cos(th), r*np.sin(th)], 1)); lab += [k]*m_
    X = np.concatenate(P); lab = np.array(lab); N = len(X)
    pool = Pool(nproc)
    def geom(X):
        tri = Delaunay(X); ip, ix = tri.vertex_neighbor_vertices
        chunks = np.array_split(np.arange(N), nproc*4); jobs = []
        for ch in chunks:
            jobs.append((ch, X[ch], [(ix[ip[a]:ip[a+1]], X[ix[ip[a]:ip[a+1]]]) for a in ch], Rc, nseg))
        V = np.zeros(N); Gs = [None]*N
        for res in pool.map(cell_job, jobs):
            for a, v, g in res: V[a] = v; Gs[a] = g
        return V, Gs
    V0, _ = geom(X); m = V0.copy(); Kent = P0
    R = np.hypot(*X.T); Om = (R*R + EPS*EPS)**-0.75; Vv = np.stack([-Om*X[:, 1], Om*X[:, 0]], 1)
    t_start = 0.0
    if len(sys.argv) > 6:
        # start from a C snapshot (same particle ordering by ring construction is NOT assumed: map by id order)
        sys.path.insert(0, '/gpfs/kjhan/LagForce/PolyBnd'); import lagsnap
        sc = lagsnap.read(sys.argv[6]); o = np.argsort(sc['id'])
        X = np.stack([sc['x'][o] - 3, sc['y'][o] - 3], 1); Vv = np.stack([sc['vx'][o], sc['vy'][o]], 1); m = sc['mass'][o]
        Kent = (G_ - 1)*sc['ie'][o]/sc['volume'][o]/(m/sc['volume'][o])**G_; t_start = sc['t']
        th0c = None
    def acc(X):
        V, Gs = geom(X); Pc = Kent*(m/V)**G_; F = np.zeros_like(X)
        for a in range(N):
            for k, v in Gs[a].items(): F[k] += Pc[a]*v
        r2 = (X*X).sum(1) + EPS*EPS
        return F/m[:, None] - X/r2[:, None]**1.5
    th0 = np.arctan2(X[:, 1], X[:, 0])
    if len(sys.argv) > 6:
        # reference angles: the C IC (ring layout) mapped by id
        s0 = lagsnap.read(sys.argv[7]); o0 = np.argsort(s0['id']); th0 = np.arctan2(s0['y'][o0] - 3, s0['x'][o0] - 3)
    print('n %d d %.4f N %d rings0/1 %d %d nseg %d' % (n, d, N, (lab == 0).sum(), (lab == 1).sum(), nseg), flush=True)
    t = t_start; t0 = time.time(); every = int(round(2.0/dt))
    for s in range(1, int(round(T/dt)) + 1):
        k1x = Vv; k1v = acc(X); k2x = Vv + 0.5*dt*k1v; k2v = acc(X + 0.5*dt*k1x)
        k3x = Vv + 0.5*dt*k2v; k3v = acc(X + 0.5*dt*k2x); k4x = Vv + dt*k3v; k4v = acc(X + dt*k3x)
        X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v); t += dt
        if s % every == 0:
            R = np.hypot(*X.T); vr = (X*Vv).sum(1)/R; vk = R*(R*R + EPS*EPS)**-0.75; th = np.arctan2(X[:, 1], X[:, 0]); out = []
            for j in (0, 1):
                sel = np.where(lab == j)[0]; o = sel[np.argsort(th0[sel])]; nn = len(o)
                dth = np.angle(np.exp(1j*(th[o] - th0[o]))); dth -= np.angle(np.mean(np.exp(1j*dth)))
                F = np.abs(np.fft.rfft(dth))*2/nn/(2*np.pi/nn); F[:2] = 0; q = np.argmax(F)
                out.append('k%d vr/vK %.1e m%d %.3f' % (j, np.sqrt(np.mean((vr[sel]/vk[sel])**2)), q, F[q]))
            print('t %5.1f | %s  [%.0fs]' % (t, ' | '.join(out), time.time() - t0), flush=True)
            np.save('/gpfs/kjhan/LagForce/PolyBnd/state_n%d%s.npy' % (n, '_fromC' if len(sys.argv) > 6 else ''), np.concatenate([X, Vv], 1))
