"""mech.py: mechanism tests for the reversible cell-scale azimuthal mode (from dynpar.py).
args: n dt T nproc nseg P0scale force(full|mid) seed_ring seed_m seed_amp(in units of d)
mid = midpoint force -A Pbar e only (tangential (Pj-Pi) t/d term dropped; not -grad U)."""
import sys, time, numpy as np
from multiprocessing import Pool
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
import shapely
G_ = 5/3; P0 = 1e-6; EPS = 0.01; RIN, ROUT = 0.5, 2.0
def J(v): return np.array([-v[1], v[0]])
def cell_job(args):
    idx, XA, NBs, Rc, nseg, fmode = args; out = []
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
                        c = 0.5*(p + q) if fmode == 'full' else mid
                        g[a] = g.get(a, 0) + L/dn*(c - xa); g[j] = g.get(j, 0) - L/dn*(c - xj); done = True; break
                if not done:
                    if np.linalg.norm(p - xa) >= rmin and np.linalg.norm(q - xa) >= rmin and np.linalg.norm(0.5*(p + q) - xa) >= rmin:
                        g[a] = g.get(a, 0) - J(q - p)
                    else: raise RuntimeError('unclassified segment, cell %d' % a)
        out.append((a, V, g))
    return out
if __name__ == '__main__':
    n = int(sys.argv[1]); dt = float(sys.argv[2]); T = float(sys.argv[3]); nproc = int(sys.argv[4]); nseg = int(sys.argv[5])
    P0 = P0*float(sys.argv[6]); fmode = sys.argv[7]; sring = int(sys.argv[8]); smode = int(sys.argv[9]); samp = float(sys.argv[10])
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
            jobs.append((ch, X[ch], [(ix[ip[a]:ip[a+1]], X[ix[ip[a]:ip[a+1]]]) for a in ch], Rc, nseg, fmode))
        V = np.zeros(N); Gs = [None]*N
        for res in pool.map(cell_job, jobs):
            for a, v, g in res: V[a] = v; Gs[a] = g
        return V, Gs
    V0, _ = geom(X); m = V0.copy(); Kent = P0
    R = np.hypot(*X.T); Om = (R*R + EPS*EPS)**-0.75; Vv = np.stack([-Om*X[:, 1], Om*X[:, 0]], 1)
    t_start = 0.0
    th0 = np.arctan2(X[:, 1], X[:, 0]); R0 = np.hypot(*X.T)
    def acc(X):
        V, Gs = geom(X); Pc = Kent*(m/V)**G_; F = np.zeros_like(X)
        for a in range(N):
            for k, v in Gs[a].items(): F[k] += Pc[a]*v
        r2 = (X*X).sum(1) + EPS*EPS
        return F/m[:, None] - (0 if PRESSURE_ONLY else X/r2[:, None]**1.5)
    PRESSURE_ONLY = True
    sel = np.where(lab == sring)[0]; Nr = len(sel); eps = 1e-4
    a0 = acc(X); Rr = R0[sel].mean(); cs2 = G_*P0; print('n %d ring %d N %d R %.3f cs/d %.4f  (1D Voronoi chain: w^2 = cs^2 sin^2(kd)/d^2; continuum cs^2 k^2)' % (n, sring, Nr, Rr, np.sqrt(cs2)/d), flush=True)
    er = np.stack([np.cos(th0[sel]), np.sin(th0[sel])], 1); et = np.stack([-np.sin(th0[sel]), np.cos(th0[sel])], 1)
    for m_ in [int(v) for v in sys.argv[12].split(',')]:
        for kind in ('tan', 'rad'):
            Xp = X.copy(); c_ = np.cos(m_*th0[sel])
            if kind == 'tan':
                thp = th0[sel] + eps*d/R0[sel]*c_; Xp[sel] = np.stack([R0[sel]*np.cos(thp), R0[sel]*np.sin(thp)], 1)
            else:
                Rp = R0[sel] + eps*d*c_; Xp[sel] = np.stack([Rp*np.cos(th0[sel]), Rp*np.sin(th0[sel])], 1)
            da = acc(Xp) - a0
            at = (da[sel]*et).sum(1); ar = (da[sel]*er).sum(1)
            proj = lambda f: 2*np.sum(f*c_)/Nr
            w2 = -proj(at if kind == 'tan' else ar)/(eps*d); kd = 2*np.pi*m_/Nr
            print('m %3d kd/pi %.2f %s: w^2/(cs^2/d^2) %+.4f | chain sin^2(kd) %.4f | continuum (kd)^2 %.3f | cross %+.3e' % (m_, kd/np.pi, kind, w2/(cs2/d/d), np.sin(kd)**2, kd*kd, -proj(ar if kind == 'tan' else at)/(eps*d)/(cs2/d/d)), flush=True)
