"""mech.py: mechanism tests for the reversible cell-scale azimuthal mode (from dynpar.py).
args: n dt T nproc nseg P0scale force(full|mid) seed_ring seed_m seed_amp(in units of d) tag [GM | 'rigid']
GM scales the point mass (shear rate); 'rigid' = uniform rotation Omega=1 with a = -x (no shear, no sliding).
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
    tag = sys.argv[11] if len(sys.argv) > 11 else 'run'; grav = sys.argv[12] if len(sys.argv) > 12 else '1'
    RIGID = (grav == 'rigid'); GM = 1.0 if RIGID else float(grav)
    th0 = np.arctan2(X[:, 1], X[:, 0]); R0 = np.hypot(*X.T)
    R = R0; Om = np.ones(N) if RIGID else np.sqrt(GM)*(R*R + EPS*EPS)**-0.75; Vv = np.stack([-Om*X[:, 1], Om*X[:, 0]], 1)
    def acc(X):
        V, Gs = geom(X); Pc = Kent*(m/V)**G_; F = np.zeros_like(X)
        for a in range(N):
            for k, v in Gs[a].items(): F[k] += Pc[a]*v
        r2 = (X*X).sum(1) + EPS*EPS
        return F/m[:, None] - (X if RIGID else GM*X/r2[:, None]**1.5)
    print('n %d d %.4f N %d rings0/1 %d %d nseg %d P0 %.1e cs %.2e cs/d %.3f force %s seed ring %d m %d amp %.3f d grav %s' % (n, d, N, (lab == 0).sum(), (lab == 1).sum(), nseg, P0, np.sqrt(G_*P0), np.sqrt(G_*P0)/d, fmode, sring, smode, samp, grav), flush=True)
    # tangent-linear perturbation: dx = seed on ring sring (mode smode, tangential), dv = 0; |dx| = eta*d
    sel = np.where(lab == sring)[0]; eta = 1e-6
    dX = np.zeros_like(X); c_ = np.cos(smode*th0[sel]); et = np.stack([-np.sin(th0[sel]), np.cos(th0[sel])], 1)
    dX[sel] = eta*d*c_[:, None]*et; dV = np.zeros_like(X); lognorm = 0.0
    def jv(Xb, dx):
        h = 1e-3*d/ (np.abs(dx).max() + 1e-300)   # FD step so that max displacement is 1e-3 d
        return (acc(Xb + h*dx) - acc(Xb - h*dx))/(2*h)
    print('tangent-linear: seed ring %d m %d, eta %.0e d, FD central step 1e-3 d' % (sring, smode, eta), flush=True)
    t = t_start; t0 = time.time(); every = int(round(1.0/dt))
    for s in range(1, int(round(T/dt)) + 1):
        k1x = Vv; k1v = acc(X); d1x = dV; d1v = jv(X, dX)
        X2 = X + 0.5*dt*k1x; k2x = Vv + 0.5*dt*k1v; k2v = acc(X2); d2x = dV + 0.5*dt*d1v; d2v = jv(X2, dX + 0.5*dt*d1x)
        X3 = X + 0.5*dt*k2x; k3x = Vv + 0.5*dt*k2v; k3v = acc(X3); d3x = dV + 0.5*dt*d2v; d3v = jv(X3, dX + 0.5*dt*d2x)
        X4 = X + dt*k3x; k4x = Vv + dt*k3v; k4v = acc(X4); d4x = dV + dt*d3v; d4v = jv(X4, dX + dt*d3x)
        X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
        dX = dX + dt/6*(d1x + 2*d2x + 2*d3x + d4x); dV = dV + dt/6*(d1v + 2*d2v + 2*d3v + d4v); t += dt
        nrm = np.sqrt((dX*dX).sum())/np.sqrt(len(sel))/d
        if nrm > 1e-4: f = 1e-6/nrm; dX *= f; dV *= f; lognorm += np.log(1/f)
        if s % every == 0:
            th = np.arctan2(X[:, 1], X[:, 0]); out = []
            for j in (sring, sring-1, sring+1, 5):
                ss = np.where(lab == j)[0]; o = ss[np.argsort(th0[ss])]; nn = len(o)
                dt_ = (dX[o]*np.stack([-np.sin(th[o]), np.cos(th[o])], 1)).sum(1); dr_ = (dX[o]*np.stack([np.cos(th[o]), np.sin(th[o])], 1)).sum(1)
                F = np.abs(np.fft.rfft(dt_))*2/nn; q = np.argmax(F[2:]) + 2
                out.append('r%d tan rms %.1e rad rms %.1e peak m%d' % (j, np.sqrt(np.mean(dt_**2))/d, np.sqrt(np.mean(dr_**2))/d, q))
            print('t %5.1f | log|dx|/|dx0| %+.3f | %s  [%.0fs]' % (t, lognorm + np.log(nrm/1e-6), ' | '.join(out), time.time() - t0), flush=True)
