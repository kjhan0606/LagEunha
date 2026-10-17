"""mechS.py: smooth-volume variational formulation test (approved 2026-10-08).
Thermodynamic volume Vt_i = sum_j K_ij V_j, K_ij = W(|x_i-x_j|/h)/S_i (cubic-spline kernel, position-only weights,
S_i = sum_j W_ij incl. j=i); U = sum_i u_i(Vt_i), fixed entropy; F = -grad U exactly:
  F_a = sum_j Q_j dV_j/dx_a  (Q_j = sum_i K_ij P~_i, back-smoothed pressure)  +  kernel term
  kernel term, per pair (i,j): f = P~_i (V_j - Vt_i) W'_ij/(h S_i) e_ij on i, -f on j.
V_j are the disc-clipped Voronoi volumes (same geometry as mech.py; later replaceable by Laguerre V_j(x,W)).
args: n dt T nproc nseg P0scale eta seed_ring seed_m seed_amp tag [GM|rigid]   (h = eta*d0; eta=0 -> plain Voronoi)"""
import sys, time, numpy as np
from multiprocessing import Pool
from scipy.spatial import Delaunay, cKDTree
from scipy.sparse import coo_matrix
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
def Wk(q):   # 2D cubic spline, unnormalised (normalisation cancels in K)
    q = np.asarray(q); return np.where(q < 0.5, 1 - 6*q*q + 6*q**3, np.where(q < 1, 2*(1 - q)**3, 0.0))
def dWk(q):
    q = np.asarray(q); return np.where(q < 0.5, -12*q + 18*q*q, np.where(q < 1, -6*(1 - q)**2, 0.0))
if __name__ == '__main__':
    n = int(sys.argv[1]); dt = float(sys.argv[2]); T = float(sys.argv[3]); nproc = int(sys.argv[4]); nseg = int(sys.argv[5])
    P0 = P0*float(sys.argv[6]); eta = float(sys.argv[7]); sring = int(sys.argv[8]); smode = int(sys.argv[9]); samp = float(sys.argv[10])
    tag = sys.argv[11]; grav = sys.argv[12] if len(sys.argv) > 12 else '1'; RIGID = (grav == 'rigid'); GM = 1.0 if RIGID else float(grav)
    d = 6.0/n; Rc = d; h = eta*d
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
    th0 = np.arctan2(X[:, 1], X[:, 0]); R0 = np.hypot(*X.T)
    sel = np.where(lab == sring)[0]; dth_seed = samp*d/R0[sel]*np.cos(smode*(th0[sel] - th0[sel][0]))
    X[sel] = np.stack([R0[sel]*np.cos(th0[sel] + dth_seed), R0[sel]*np.sin(th0[sel] + dth_seed)], 1)
    R = np.hypot(*X.T); Om = np.ones(N) if RIGID else np.sqrt(GM)*(R*R + EPS*EPS)**-0.75; Vv = np.stack([-Om*X[:, 1], Om*X[:, 0]], 1)
    def smooth(X, V):
        """returns Vt, and the sparse pair data (i, j, Wij, dWij/dr, e_ij) for the kernel force"""
        if h <= 0: return V.copy(), None
        tree = cKDTree(X); pairs = tree.query_pairs(2*h, output_type='ndarray')   # i<j within 2h
        i, j = pairs[:, 0], pairs[:, 1]; dx = X[j] - X[i]; r = np.hypot(*dx.T); e = dx/r[:, None]
        w = Wk(r/(2*h)); dw = dWk(r/(2*h))/(2*h)   # W as a function of q = r/(2h), support r < 2h
        S = np.full(N, Wk(0.0)) + np.bincount(i, w, N) + np.bincount(j, w, N)
        B = Wk(0.0)*V + np.bincount(i, w*V[j], N) + np.bincount(j, w*V[i], N)
        return B/S, (i, j, w, dw, e, S)
    def acc(X):
        V, Gs = geom(X); Vt, pr = smooth(X, V); Pt = Kent*(m/Vt)**G_
        if pr is None: Q = Pt
        else:
            i, j, w, dw, e, S = pr
            # Q_j = sum_i K_ij Pt_i = sum_i W_ij Pt_i / S_i  (incl. i=j)
            Q = Wk(0.0)*Pt/S + np.bincount(j, w*Pt[i]/S[i], N) + np.bincount(i, w*Pt[j]/S[j], N)
        F = np.zeros_like(X)
        for a in range(N):
            for k, v in Gs[a].items(): F[k] += Q[a]*v
        if pr is not None:
            # kernel term: pair (i,j): on i: Pt_i (V_j - Vt_i) dW/dr /S_i * e_ij ; on j the reaction; plus roles swapped
            fi = (Pt[i]*(V[j] - Vt[i])/S[i] + Pt[j]*(V[i] - Vt[j])/S[j])*dw   # along e_ij = (x_j - x_i)/r
            # dW_ij/dx_i = dW/dr * (x_i - x_j)/r = -dw e ; contribution to F_i = Pt_i (V_j - Vt_i)/S_i * (-dw e)
            fvec = -fi[:, None]*e
            np.add.at(F, i, fvec); np.add.at(F, j, -fvec)
        r2 = (X*X).sum(1) + EPS*EPS
        return F/m[:, None] - (X if RIGID else GM*X/r2[:, None]**1.5), Vt
    def Etot(X, Vv):
        _, Vt = acc(X) if False else (None, None)
        return None
    print('n %d d %.4f N %d rings0/1 %d %d P0 %.1e cs/d %.3f eta %.2f (kernel radius %.1f d) seed ring %d m %d amp %.3f d grav %s' % (
        n, d, N, (lab == 0).sum(), (lab == 1).sum(), P0, np.sqrt(G_*P0)/d, eta, 2*eta, sring, smode, samp, grav), flush=True)
    t = 0.0; t0 = time.time(); every = int(round(1.0/dt))
    # energy check at t=0: U(Vt) + KE + Phi
    V, Gs = geom(X); Vt0, _ = smooth(X, V); U0 = (Kent*m**G_*Vt0**(1 - G_)/(G_ - 1)).sum()
    def energy(X, Vv):
        V, _ = geom(X); Vt, _ = smooth(X, V); U = (Kent*m**G_*Vt**(1 - G_)/(G_ - 1)).sum()
        r2 = (X*X).sum(1) + EPS*EPS; Phi = (0.5*m*(X*X).sum(1)).sum() if RIGID else -(GM*m/np.sqrt(r2)).sum()
        return 0.5*(m*(Vv*Vv).sum(1)).sum() + U + Phi
    E0 = energy(X, Vv)
    for s in range(1, int(round(T/dt)) + 1):
        k1x = Vv; k1v, _ = acc(X); k2x = Vv + 0.5*dt*k1v; k2v, _ = acc(X + 0.5*dt*k1x)
        k3x = Vv + 0.5*dt*k2v; k3v, _ = acc(X + 0.5*dt*k2x); k4x = Vv + dt*k3v; k4v, _ = acc(X + dt*k3x)
        X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v); t += dt
        if s % every == 0:
            th = np.arctan2(X[:, 1], X[:, 0]); out = []
            for jr in (sring, 5):
                ss = np.where(lab == jr)[0]; o = ss[np.argsort(th0[ss])]; nn = len(o)
                dth = np.angle(np.exp(1j*(th[o] - th0[o]))); dth -= dth.mean()
                F = np.abs(np.fft.rfft(dth))*2/nn; F[:2] = 0; q = np.argmax(F)
                out.append('ring%d seed m%d R*dth/d %.2e | max m%d %.2e' % (jr, smode, np.hypot(*X[o].T).mean()*F[smode]/d if smode < len(F) else 0, q, np.hypot(*X[o].T).mean()*F[q]/d))
            print('t %5.1f | %s | dE/E0 %+.2e  [%.0fs]' % (t, ' | '.join(out), (energy(X, Vv) - E0)/abs(E0), time.time() - t0), flush=True)
            np.save('/gpfs/kjhan/LagForce/PolyBnd/mechS_%s.npy' % tag, np.concatenate([X, Vv], 1))
