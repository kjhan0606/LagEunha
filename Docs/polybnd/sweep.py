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
    sel = np.where(lab == sring)[0]; Nr = len(sel); eps = 1e-4; cs2 = G_*P0
    er = np.stack([np.cos(th0[sel]), np.sin(th0[sel])], 1); et = np.stack([-np.sin(th0[sel]), np.cos(th0[sel])], 1)
    modes = [int(v) for v in sys.argv[12].split(',')]; phis = [float(v) for v in sys.argv[13].split(',')]
    print('phase sweep: rings %d and %d shifted by phi*spacing; ring %d N %d; w^2 in units cs^2/d^2' % (sring-1, sring+1, sring, Nr), flush=True)
    acc_all = {}
    for phi in phis:
        Xb = X.copy()
        for rr in (sring-1, sring+1):
            ss = np.where(lab == rr)[0]; dphi = phi*2*np.pi/len(ss); thp = th0[ss] + dphi
            Xb[ss] = np.stack([R0[ss]*np.cos(thp), R0[ss]*np.sin(thp)], 1)
        a0 = acc(Xb); row = []
        for m_ in modes:
            Xp = Xb.copy(); c_ = np.cos(m_*th0[sel]); thp = th0[sel] + eps*d/R0[sel]*c_
            Xp[sel] = np.stack([R0[sel]*np.cos(thp), R0[sel]*np.sin(thp)], 1)
            da = acc(Xp) - a0; at = (da[sel]*et).sum(1)
            w2 = -2*np.sum(at*c_)/Nr/(eps*d)/(cs2/d/d); row.append(w2); acc_all.setdefault(m_, []).append(w2)
        print('phi %.2f | ' % phi + ' '.join('m%d %+.4f' % (m_, w) for m_, w in zip(modes, row)), flush=True)
    print('MEAN   | ' + ' '.join('m%d %+.4f' % (m_, np.mean(acc_all[m_])) for m_ in modes), flush=True)
