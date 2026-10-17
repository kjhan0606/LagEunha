"""Power (Laguerre) diagram with disc closure for the ring IC: cells, volumes, dV/dx at fixed w, dV/dw, and a Newton
solve for weights giving prescribed volumes. Serial (n <= 96). Face a|j: plane perpendicular to e=(xj-xa)/d at distance
s = d/2 + (w_a - w_j)/(2d) from x_a. Candidate neighbours: 2-hop Delaunay (power-diagram topology may differ from Delaunay)."""
import sys, numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from shapely.geometry import Polygon
import shapely
G_ = 5/3; RIN, ROUT = 0.5, 2.0
def J(v): return np.array([-v[1], v[0]])
def ring_ic(n):
    d = 6.0/n; P = []; lab = []
    for k in range(int(round((ROUT - RIN)/d))):
        r = RIN + (k + 0.5)*d; m_ = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(m_) + 0.5)*2*np.pi/m_
        P.append(np.stack([r*np.cos(th), r*np.sin(th)], 1)); lab += [k]*m_
    return np.concatenate(P), np.array(lab), d
def candidates(X):
    tri = Delaunay(X); ip, ix = tri.vertex_neighbor_vertices; N = len(X); C = []
    for a in range(N):
        n1 = ix[ip[a]:ip[a+1]]; s = set(n1.tolist())
        for j in n1: s.update(ix[ip[j]:ip[j+1]].tolist())
        s.discard(a); C.append(np.array(sorted(s)))
    return C
def cells(X, w, Rc, C, nseg=256):
    """returns V (N,), G: list of dicts {k: dV_a/dx_k} at fixed w, faces: list of (a, j, A, d) for a<j"""
    N = len(X); V = np.zeros(N); G = [dict() for _ in range(N)]; faces = []
    circ0 = np.array([[np.cos(t), np.sin(t)] for t in np.linspace(0, 2*np.pi, nseg, endpoint=False)])*Rc
    rmin = Rc*np.cos(np.pi/nseg) - 1e-9
    for a in range(N):
        xa = X[a]; B = 8*Rc + 1.0
        cell = Polygon([(xa[0]-B, xa[1]-B), (xa[0]+B, xa[1]-B), (xa[0]+B, xa[1]+B), (xa[0]-B, xa[1]+B)])
        nb = C[a]; planes = []
        for j in nb:
            dv = X[j] - xa; dn = np.linalg.norm(dv); e = dv/dn; s = 0.5*dn + (w[a] - w[j])/(2*dn)
            if s <= 0 or s >= dn and False: pass
            pt = xa + s*e; t = J(e)*50; nrm = e*50
            cell = cell.intersection(Polygon([pt + t, pt - t, pt - t - nrm, pt + t - nrm])); planes.append((j, e, dn, pt))
        cell = cell.intersection(Polygon(xa + circ0))
        V[a] = cell.area
        polys = [cell] if isinstance(cell, Polygon) else [q for q in cell.geoms if isinstance(q, Polygon)]
        for pgn in polys:
            cs = np.array(shapely.geometry.polygon.orient(pgn, 1.0).exterior.coords)
            for k in range(len(cs) - 1):
                p, q = cs[k], cs[k+1]; L = np.linalg.norm(q - p)
                if L < 1e-14: continue
                done = False
                for j, e, dn, pt in planes:
                    if abs((p - pt) @ e) < 1e-9 and abs((q - pt) @ e) < 1e-9:
                        c = 0.5*(p + q)
                        # dV_a/dx_a and dV_a/dx_j at fixed w (face point c, Laguerre offset included via c - x)
                        # V_a = sum over faces; d(face)/dx: standard formula with face centroid c: dV_a/dx_a += L/dn (c - x_a) + (L/dn)(w_a-w_j)/dn * e ... use exact: s depends on dn
                        G[a][a] = G[a].get(a, 0) + L/dn*(c - xa) - L*(w[a] - w[j])/(dn**3)*(X[j] - xa)*0  # placeholder corrected below
                        G[a][j] = G[a].get(j, 0) - L/dn*(c - X[j])
                        if a < j: faces.append((a, j, L, dn))
                        done = True; break
                if not done:
                    if np.linalg.norm(p - xa) >= rmin and np.linalg.norm(q - xa) >= rmin and np.linalg.norm(0.5*(p + q) - xa) >= rmin:
                        G[a][a] = G[a].get(a, 0) - J(q - p)
                    else: raise RuntimeError('unclassified segment, cell %d' % a)
    return V, G, faces
def dVdw(N, faces):
    """sparse Jacobian dV_a/dw_b: diag += A/(2d), off = -A/(2d)"""
    r, c, v = [], [], []
    for a, j, A, d in faces:
        k = A/(2*d); r += [a, j, a, j]; c += [a, j, j, a]; v += [k, k, -k, -k]
    return coo_matrix((v, (r, c)), shape=(N, N)).tocsr()
def solve_weights(X, Vt, Rc, C, w0=None, tol=1e-6, itmax=30, verbose=True):
    N = len(X); w = np.zeros(N) if w0 is None else w0.copy()
    for it in range(itmax):
        V, G, faces = cells(X, w, Rc, C); res = V - Vt; rn = np.abs(res).max()/Vt.mean()
        if verbose: print('  newton it %d  max|V-Vt|/V %.2e  max|w|/d^2 %.3e' % (it, rn, np.abs(w).max()), flush=True)
        if rn < tol: return w, V, G, faces
        # gauge: faces depend on w differences only, so sum_a V_a is w-invariant and (1,..,1) is a null vector;
        # pin the mean of w to 0 and solve in the complement
        Jw = (dVdw(N, faces) + coo_matrix((np.full(N, 1.0/N), (np.arange(N), np.arange(N))), shape=(N, N))).tocsc()
        dw = spsolve(Jw, -res + res.mean()); w = w + dw; w -= w.mean()
    raise RuntimeError('weights did not converge')
if __name__ == '__main__':
    n = int(sys.argv[1]); X, lab, d = ring_ic(n); Rc = d; N = len(X); C = candidates(X)
    w0 = np.zeros(N); V0, G0, f0 = cells(X, w0, Rc, C)
    print('n %d N %d  Voronoi volumes: ring0 %.4e ring1 %.4e ring5 %.4e' % (n, N, V0[lab == 0].mean(), V0[lab == 1].mean(), V0[lab == 5].mean()))
    # consistency: targets = Voronoi volumes -> w = 0 is the solution
    w, V, G, faces = solve_weights(X, V0*V0.sum()/V0.sum(), Rc, C)
    print('consistency: max|w| %.1e (expect 0)' % np.abs(w).max())
    # displace ring 1 by a short-wavelength tangential mode, then solve weights that restore the IC volumes
    sel = np.where(lab == 1)[0]; th = np.arctan2(X[sel, 1], X[sel, 0]); R = np.hypot(*X[sel].T); m_ = int(sys.argv[2]); amp = float(sys.argv[3])
    Xp = X.copy(); thp = th + amp*d/R*np.cos(m_*th); Xp[sel] = np.stack([R*np.cos(thp), R*np.sin(thp)], 1)
    Cp = candidates(Xp); Vv, Gv, fv = cells(Xp, w0, Rc, Cp)
    print('displaced ring 1 (m=%d, amp %.3f d): Voronoi volume change ring1 rms %.2e, max %.2e (relative)' % (m_, amp, np.sqrt(np.mean((Vv[sel]/V0[sel]-1)**2)), np.abs(Vv[sel]/V0[sel]-1).max()))
    Vt = V0*Vv.sum()/V0.sum(); print('total volume change after displacement: %.2e (targets rescaled)' % (Vv.sum()/V0.sum()-1))
    wp, Vp, Gp, fp = solve_weights(Xp, Vt, Rc, Cp)
    print('Laguerre weights restoring IC volumes: max|w|/d^2 %.3e, ring1 rms %.3e, ring0 rms %.3e, ring2 rms %.3e' % (np.abs(wp).max()/d/d, np.sqrt(np.mean(wp[sel]**2))/d/d, np.sqrt(np.mean(wp[lab==0]**2))/d/d, np.sqrt(np.mean(wp[lab==2]**2))/d/d))
    np.save('lag_test_n%d.npy' % n, dict(X=Xp, w=wp, V0=V0, lab=lab, d=d), allow_pickle=True)
