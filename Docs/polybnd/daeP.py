"""daeP.py: Laguerre-DAE (or plain Voronoi, mode=vor) in a PERIODIC box with the non-vacuum Kepler setup (disc rho=1 on
0.5<=r<=2, background rho=BG, uniform P, spline-softened point mass rc=0.5, R>=2.5 slaved/frozen). No disc clipping.
args: n dt T nproc P seed_R seed_m seed_amp tag divmode(face|flux2|vor) BG
Derived from dae.py (approved 2026-10-08).
State: X, V (velocity), VT (thermodynamic volume), U (internal energy); W algebraic.
Each RHS evaluation at (X, V, VT, U):
  1. project W: V(X, W) = VT (Newton, dV/dW = 1/2 L_{A/d}, zero-mean gauge, residual mean removed and logged)
  2. P_a = (g-1) U_a / VT_a
  3. force F_a = sum_b P_b dV_b/dx_a |_W   (exact for power diagrams: dV_a/dx_a = A/d (c - x_a), dV_a/dx_j = -A/d (c - x_j))
  4. VT'_a = sum_faces A (u_a+u_j)/2 . e  + arcs A u_a.n   (material/midpoint face kinematics)
  5. Wdot from the constraint tangent: 1/2 L Wdot = VT' - dV/dt|_W
  6. U'_a = -P_a VT'_a + 1/2 sum_faces D_f,  D_f = A (P_a - P_j) (Wdot_a - Wdot_j)/(2d)   (W-work, split 1/2-1/2) -> total energy exact
Diagnostics: seeded-mode growth, max|dW|/d^2, saturated faces (|dW|/d^2 > 0.9), E_tot residual, negative D_f share, K drift.
args: n dt T nproc nseg P0scale seed_ring seed_m seed_amp tag [edgew0]
divmode=flux3: like flux2 with CUBIC reconstruction (9 coefficients) and 2-point Gauss quadrature along each face/arc segment.
divmode=flux2: VT'_a = sum_faces A u_f.e with u_f = 1/2 (v~_a(c_f) + v~_j(c_f)) from the cells' quadratic LSQ reconstructions at the face centroid (conservative, high order); arcs use v~_a at the segment midpoint.
divmode=lsq (rejected, non-conservative: total VT drifts): VT'_a = VT_a * tr(G_a) with G_a the quadratic least-squares velocity gradient (2-hop stencil) instead of midpoint face kinematics.
tauv>0 (argv[14]): VT' += (V_geo - VT)/tauv — low-pass relaxation of the material volume toward the geometric volume (scheme H Phase A tau_V)
edgew0=1: arc (vacuum-edge) cells are slaved: W_a = 0 pinned, VT_a := V_a (no constraint row), VT'_a = geometric rate incl. neighbour Wdot (energy exact)."""
import sys, time, numpy as np
from multiprocessing import Pool
from scipy.spatial import Delaunay, cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from shapely.geometry import Polygon
import shapely
G_ = 5/3; P0 = 1e-6; EPS = 0.01; RIN, ROUT = 0.5, 2.0
def J(v): return np.array([-v[1], v[0]])
def cell_job(args):
    idx, XA, WA, NBs, Rc, nseg = args; out = []
    for a, xa, wa, (nbid, XN, WN) in zip(idx, XA, WA, NBs):
        B = 8*Rc + 1.0
        cell = Polygon([(xa[0]-B, xa[1]-B), (xa[0]+B, xa[1]-B), (xa[0]+B, xa[1]+B), (xa[0]-B, xa[1]+B)])
        planes = []
        for j, xj, wj in zip(nbid, XN, WN):
            dv = xj - xa; dn = np.linalg.norm(dv); e = dv/dn; s = 0.5*dn + (wa - wj)/(2*dn); pt = xa + s*e
            t = J(e)*50; nrm = e*50
            cell = cell.intersection(Polygon([pt + t, pt - t, pt - t - nrm, pt + t - nrm])); planes.append((j, e, dn, pt, xj))
        V = cell.area; g = {}; fl = []
        polys = [cell] if isinstance(cell, Polygon) else [q for q in cell.geoms if isinstance(q, Polygon)]
        for pgn in polys:
            cs = np.array(shapely.geometry.polygon.orient(pgn, 1.0).exterior.coords)
            for k in range(len(cs) - 1):
                p, q = cs[k], cs[k+1]; L = np.linalg.norm(q - p)
                if L < 1e-14: continue
                done = False
                for j, e, dn, pt, xj in planes:
                    if abs((p - pt) @ e) < 1e-9 and abs((q - pt) @ e) < 1e-9:
                        c = 0.5*(p + q); g[a] = g.get(a, 0) + L/dn*(c - xa); g[j] = g.get(j, 0) - L/dn*(c - xj)
                        fl.append((j, L, e, dn, c, p, q)); done = True; break
                if not done: raise RuntimeError('unclassified segment, cell %d (candidate radius too small?)' % a)
        out.append((a, V, g, fl))
    return out
if __name__ == '__main__':
    n = int(sys.argv[1]); dt = float(sys.argv[2]); T = float(sys.argv[3]); nproc = int(sys.argv[4]); nseg = 0
    P0 = float(sys.argv[5]); sR = float(sys.argv[6]); smode = int(sys.argv[7]); samp = float(sys.argv[8]); tag = sys.argv[9]
    DIVMODE = sys.argv[10]; BG = float(sys.argv[11]); EDGEW0 = 1; DBG = (len(sys.argv) > 12 and sys.argv[12] == 'dbg'); TAUV = 0.0
    L = 6.0; C = 3.0; d = L/n; Rc = d; RC = 0.5; RFRZ = 2.5; RFRZ_IN = 0.25   # inner core (background only, within the spline core) frozen too
    RT1, RT2 = 2.5, 2.9   # taper: rotation and gravity -> 0 between RT1 and RT2 (static, periodic-safe far background)
    def taper(r):
        r = np.asarray(r, float); x = np.clip((r - RT1)/(RT2 - RT1), 0, 1); return 0.5*(1 + np.cos(np.pi*x))
    # uniform rings, golden phase, periodic-safe (as ICgraded/mkrings_uniform.py)
    xs = []; ys = []; k = 1
    while k*d < 3*np.sqrt(2) + d:
        r = k*d; m_ = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(m_) + (0.6180339887*k) % 1.0)*2*np.pi/m_
        xs.append(C + r*np.cos(th)); ys.append(C + r*np.sin(th)); k += 1
    X = np.stack([np.concatenate(xs), np.concatenate(ys)], 1); inb = (X[:, 0] > 0) & (X[:, 0] < L) & (X[:, 1] > 0) & (X[:, 1] < L); X = X[inb]
    Tq = cKDTree(X, boxsize=[L, L]); drop = set(b for a_, b in Tq.query_pairs(0.7*d)); keep = np.ones(len(X), bool); keep[list(drop)] = False; X = X[keep]
    N = len(X); R0 = np.hypot(X[:, 0] - C, X[:, 1] - C); lab = np.round(R0/d).astype(int)
    def Mspl(u):
        u = np.asarray(u, float)
        return np.where(u >= 1, 1.0, np.where(u <= 0.5, u**3*(32/3 - 192/5*u*u + 32*u**3), -1/15 + u**3*(64/3 - 48*u + 192/5*u*u - 32/3*u**3)))
    def gacc(X):
        dx = X - C; r = np.hypot(*dx.T); g = -Mspl(r/RC)/r**3*taper(r)**2; return g[:, None]*dx   # g = v^2/r with v = v_K*taper: exact circular balance everywhere
    _rg = np.linspace(1e-4, 5.0, 20001); _gg = Mspl(_rg/RC)/_rg**2*taper(_rg)**2; _ph = np.concatenate([[0.0], np.cumsum(0.5*(_gg[1:] + _gg[:-1])*np.diff(_rg))]); _ph -= _ph[-1]   # Phi(r) = -int_r^inf g dr, g=0 beyond RT2
    def Phi(X):
        dx = X - C; r = np.hypot(*dx.T); return np.interp(r, _rg, _ph)
    Sig = np.where((R0 >= 0.5) & (R0 <= 2.0), 1.0, BG)
    FREEZE = int(sys.argv[13]) if len(sys.argv) > 13 else 1
    frz = ((R0 >= RFRZ) | (R0 < RFRZ_IN)) if FREEZE else np.zeros(N, bool)   # FREEZE=0: fully periodic, everything evolves (no artificial interfaces)
    EDGEW0 = 1 if frz.any() else 0   # without pinned (frozen) cells the Laplacian needs the zero-mean gauge
    pool = Pool(nproc)
    def geom(X, W):
        # exact periodic Voronoi neighbours: Delaunay of the points plus image copies of boundary-near points
        Xw = np.mod(X, L); marg = 5*d; src = [np.arange(N)]; pts = [Xw]
        for sx in (-1, 0, 1):
            for sy in (-1, 0, 1):
                if sx == 0 and sy == 0: continue
                sel_ = np.where(((sx == 0) | ((Xw[:, 0] < marg) if sx > 0 else (Xw[:, 0] > L - marg))) & ((sy == 0) | ((Xw[:, 1] < marg) if sy > 0 else (Xw[:, 1] > L - marg))))[0]
                if len(sel_): pts.append(Xw[sel_] + np.array([sx*L, sy*L])); src.append(sel_)
        Q = np.concatenate(pts); srcidx = np.concatenate(src); tri = Delaunay(Q); ip, ix = tri.vertex_neighbor_vertices
        nb = []; nbx = []
        for a in range(N):
            js = ix[ip[a]:ip[a+1]]; nb.append(srcidx[js]); nbx.append(Q[js])
        # LSQ stencil (2-hop by real index, minimum image handled in lsq_coef)
        geom.nb = [np.unique(np.concatenate([nb[a]] + [nb[j] for j in nb[a]]))[lambda u: u != a] if False else np.array([j for j in np.unique(np.concatenate([nb[a]] + [nb[j] for j in nb[a]])) if j != a]) for a in range(N)]
        chunks = np.array_split(np.arange(N), nproc*4); jobs = []
        for ch in chunks:
            jobs.append((ch, Xw[ch], W[ch], [(nb[a], nbx[a], W[nb[a]]) for a in ch], Rc, nseg))
        V = np.zeros(N); Gs = [None]*N; FL = [None]*N
        for res in pool.map(cell_job, jobs):
            for a, v, g, fl in res: V[a] = v; Gs[a] = g; FL[a] = fl
        return V, Gs, FL
    def arcflag(FL): return frz.copy()   # slaved (frozen) far background instead of arc cells
    def lsq_coef(X, Vv):
        # per cell: c (5,2): [Gx, Gy, Hxx, Hxy, Hyy] for (vx, vy); v~(x+y) = v + c0 yx + c1 yy + 1/2 c2 yx^2 + c3 yx yy + 1/2 c4 yy^2
        nb = geom.nb; Cc = np.zeros((N, 5, 2))
        for a in range(N):
            ys = X[nb[a]] - X[a]; ys -= L*np.round(ys/L); w = 1.0/(ys*ys).sum(1); sw = np.sqrt(w)
            A = np.stack([ys[:, 0], ys[:, 1], 0.5*ys[:, 0]**2, ys[:, 0]*ys[:, 1], 0.5*ys[:, 1]**2], 1)*sw[:, None]
            b = (Vv[nb[a]] - Vv[a])*sw[:, None]
            Cc[a], *_ = np.linalg.lstsq(A, b, rcond=None)
        return Cc
    def lsq_coef3(X, Vv):
        # cubic LSQ: basis [x, y, x2/2, xy, y2/2, x3/6, x2y/2, xy2/2, y3/6], weights 1/|y|^2
        nb = geom.nb; C = np.zeros((N, 9, 2))
        for a in range(N):
            ys = X[nb[a]] - X[a]; w = 1.0/(ys*ys).sum(1); sw = np.sqrt(w); x_, y_ = ys[:, 0], ys[:, 1]
            A = np.stack([x_, y_, 0.5*x_*x_, x_*y_, 0.5*y_*y_, x_**3/6, 0.5*x_*x_*y_, 0.5*x_*y_*y_, y_**3/6], 1)*sw[:, None]
            b = (Vv[nb[a]] - Vv[a])*sw[:, None]
            C[a], *_ = np.linalg.lstsq(A, b, rcond=None)
        return C
    def vrec3(a, y, Vv, C):
        c = C[a]; x_, y_ = y[0], y[1]
        return Vv[a] + c[0]*x_ + c[1]*y_ + 0.5*c[2]*x_*x_ + c[3]*x_*y_ + 0.5*c[4]*y_*y_ + c[5]*x_**3/6 + 0.5*c[6]*x_*x_*y_ + 0.5*c[7]*x_*y_*y_ + c[8]*y_**3/6
    def vrec(a, y, Vv, C):
        c = C[a]; return Vv[a] + c[0]*y[0] + c[1]*y[1] + 0.5*c[2]*y[0]**2 + c[3]*y[0]*y[1] + 0.5*c[4]*y[1]**2
    def lsq_div(X, Vv):
        # quadratic LSQ fit of (vx, vy) on the 2-hop stencil, weights 1/|y|^2; returns tr G
        nb = geom.nb; div = np.zeros(N)
        for a in range(N):
            ys = X[nb[a]] - X[a]; w = 1.0/(ys*ys).sum(1); sw = np.sqrt(w)
            A = np.stack([ys[:, 0], ys[:, 1], 0.5*ys[:, 0]**2, ys[:, 0]*ys[:, 1], 0.5*ys[:, 1]**2], 1)*sw[:, None]
            b = (Vv[nb[a]] - Vv[a])*sw[:, None]
            c, *_ = np.linalg.lstsq(A, b, rcond=None)
            div[a] = c[0, 0] + c[1, 1]
        return div
    def Lmat(FL):
        r, cc, v = [], [], []
        for a in range(N):
            for j, A, e, dn, c_, p_, q_ in FL[a]:
                if j == a: continue
                k = A/(2*dn); r += [a, a]; cc += [a, j]; v += [k, -k]
        return coo_matrix((v, (r, cc)), shape=(N, N)).tocsr()
    Pin = coo_matrix((np.full(N, 1.0/N), (np.arange(N), np.arange(N))), shape=(N, N)).tocsr()   # gauge pin (mean)
    stat = dict(newton=0, resid=0.0, meanres=0.0)
    def project(X, W, VT, tol=1e-7, itmax=12):
        for it in range(itmax):
            V, Gs, FL = geom(X, W); res = V - VT
            arc = arcflag(FL) if EDGEW0 else np.zeros(N, bool)
            if EDGEW0 and it == 0: VT = VT.copy(); VT[arc] = V[arc]; res = V - VT   # slave once; keep the target fixed during Newton
            # with pinned arc cells the system is nonsingular (arc cells are the volume reservoir): no gauge, no mean removal
            mres = 0.0 if EDGEW0 else res.mean(); res0 = res - mres; res0[arc] = 0
            rn = np.abs(res0).max()/VT.mean(); stat['newton'] += 1
            if DBG: print('   newton it %d rn %.2e mres %.2e |W|max/d2 %.3e narc %d' % (it, rn, mres/VT.mean(), np.abs(W).max()/d/d, arc.sum()), flush=True)
            if rn < tol:
                if EDGEW0: VT = VT.copy(); VT[arc] = V[arc]
                stat['resid'] = rn; stat['meanres'] = mres/VT.mean(); return W, V, Gs, FL, VT, arc
            Lw = Lmat(FL).tolil()
            if EDGEW0:
                for a in np.where(arc)[0]: Lw.rows[a] = [a]; Lw.data[a] = [1.0]
            Lw = (Lw.tocsr() + (0*Pin if EDGEW0 else Pin)).tocsc()
            if DBG: print('   it %d rn %.2e |res0|max %.2e nan(res) %d nan(Lw) %d diag0 %d VT min %.2e V min %.2e nonfrozen-zero-rows %d' % (it, rn, np.abs(res0).max(), int(np.isnan(res0).sum()), int(np.isnan(Lw.data).sum()), int((Lw.diagonal()==0).sum()), VT.min(), V.min(), int(((Lw.diagonal()==0) & ~arc).sum())), flush=True)
            W = W + spsolve(Lw, -res0); W[arc] = 0
            if not EDGEW0: W -= W.mean()
        if rn < 1e-5:   # arc-polygonisation floor: accept and count
            stat['loose'] = stat.get('loose', 0) + 1; stat['resid'] = rn; stat['meanres'] = mres/VT.mean()
            if EDGEW0: VT = VT.copy(); VT[arc] = V[arc]
            return W, V, Gs, FL, VT, arc
        raise RuntimeError('projection did not converge: %.2e' % rn)
    def rhs(X, Vv, VT, U, W):
        if (VT[~frz] <= 0).any(): raise RuntimeError('VT <= 0 in a stage (min %.3e): reduce dt' % VT[~frz].min())
        if DIVMODE == 'vor':
            V, Gs, FL = geom(X, W*0); VT = V.copy(); arc = frz.copy(); W = W*0
        else: W, V, Gs, FL, VT, arc = project(X, W, VT)
        Pc = (G_ - 1)*U/VT; F = np.zeros_like(X); dVgeo = np.zeros(N); dVT = np.zeros(N)
        for a in range(N):
            for k, g in Gs[a].items(): F[k] += Pc[a]*g; dVgeo[a] += g @ Vv[k]
            for j, A, e, dn, c, p_, q_ in FL[a]: dVT[a] += A*0.5*((Vv[a] + Vv[j]) @ e)
        if DIVMODE == 'flux2':
            C = lsq_coef(X, Vv); dVT = np.zeros(N)
            for a in range(N):
                for j, A, e, dn, c, p_, q_ in FL[a]:
                    ya = c - np.mod(X[a], L); yj = c - np.mod(X[j], L); yj -= L*np.round(yj/L)
                    u = 0.5*(vrec(a, ya, Vv, C) + vrec(j, yj, Vv, C))
                    dVT[a] += A*(u @ e)
        if TAUV > 0:
            Vvor, _, _ = geom(X, np.zeros(N))   # relaxation target = W=0 Voronoi volume (V(x,W) == VT by construction)
            dVT = dVT + (Vvor - VT)/TAUV
        if DBG: print('   dVT*dt/V: rms %.2e max %.2e (ring of max: %d)  dVgeo*dt/V rms %.2e' % (np.sqrt(np.mean((dVT*dt/V)**2)), np.abs(dVT*dt/V).max(), lab[np.argmax(np.abs(dVT/V))], np.sqrt(np.mean((dVgeo*dt/V)**2))), flush=True)
        if DIVMODE == 'vor': dVT = dVgeo.copy()
        stat['FL'] = FL; stat['W'] = W
        Lm = Lmat(FL); rr = dVT - dVgeo; rr0 = rr - (0.0 if EDGEW0 else rr.mean()); rr0[arc] = 0
        Lw = Lm.tolil()
        if EDGEW0:
            for a in np.where(arc)[0]: Lw.rows[a] = [a]; Lw.data[a] = [1.0]
        Lw = (Lw.tocsr() + (0*Pin if EDGEW0 else Pin)).tocsc(); Wd = spsolve(Lw, rr0); Wd[arc] = 0
        if EDGEW0: dVT = dVT.copy(); dVT[arc] = (dVgeo + Lm @ Wd)[arc]   # slaved cells: true geometric rate
        dU = -Pc*dVT; Dneg = 0.0; Dpos = 0.0
        for a in range(N):
            for j, A, e, dn, c, p_, q_ in FL[a]:
                if j == a: continue
                Df = A*(Pc[a] - Pc[j])*(Wd[a] - Wd[j])/(2*dn); dU[a] += 0.5*Df
                if Df < 0: Dneg += 0.5*Df
                else: Dpos += 0.5*Df
        acc = F/m[:, None] + gacc(X)
        dW = np.array([abs(W[a] - W[j])/(dn*dn) for a in range(N) for j, A, e, dn, c, p_, q_ in FL[a] if j != a])
        return acc, dVT, dU, W, dict(maxdW=dW.max(), sat=np.mean(dW > 0.9), Dneg=Dneg, Dpos=Dpos, dVmis=np.sqrt(np.mean((dVT - dVgeo)**2))*dt/V.mean(), VT=VT)
    W = np.zeros(N); V0, G0, FL0 = geom(X, W); m = Sig*V0; VT = V0.copy(); U = P0*VT/(G_ - 1); K0 = P0/(m/VT)**G_
    th0 = np.arctan2(X[:, 1] - C, X[:, 0] - C)
    sel = np.where(np.abs(R0 - sR) < 0.5*d)[0]; dth = samp*d/R0[sel]*np.cos(smode*(th0[sel] - th0[sel][0]))
    X[sel] = np.stack([C + R0[sel]*np.cos(th0[sel] + dth), C + R0[sel]*np.sin(th0[sel] + dth)], 1)
    R = np.hypot(X[:, 0] - C, X[:, 1] - C); vc = np.sqrt(Mspl(R/RC)/R)*taper(R); Vv = np.stack([-vc*(X[:, 1] - C)/R, vc*(X[:, 0] - C)/R], 1)
    if FREEZE: Vv[frz] = 0.0   # frozen regions: at rest, slaved (W=0, VT=V)
    sring = -1
    def Etot(X, Vv, U):
        ev = ~frz; return 0.5*(m*(Vv*Vv).sum(1))[ev].sum() + U[ev].sum() + (m*Phi(X))[ev].sum()
    E0 = Etot(X, Vv, U)
    print('daeP n %d N %d (frozen %d) P0 %.2e cs_disc/d %.3f cs_bg/d %.3f BG %.0e seed R %.2f (N %d) m %d amp %.3f d  E0 %.6e divmode %s' % (n, N, frz.sum(), P0, np.sqrt(G_*P0)/d, np.sqrt(G_*P0/BG)/d, BG, sR, len(sel), smode, samp, E0, DIVMODE), flush=True)
    t = 0.0; t0 = time.time(); every = int(round(1.0/dt)); Dneg_c = 0.0; Dpos_c = 0.0
    for s in range(1, int(round(T/dt)) + 1):
        a1, t1, u1, W, d1 = rhs(X, Vv, VT, U, W); VT = d1['VT']; k1x = Vv
        a2, t2, u2, W, d2 = rhs(X + 0.5*dt*k1x, Vv + 0.5*dt*a1, VT + 0.5*dt*t1, U + 0.5*dt*u1, W); k2x = Vv + 0.5*dt*a1
        a3, t3, u3, W, d3 = rhs(X + 0.5*dt*k2x, Vv + 0.5*dt*a2, VT + 0.5*dt*t2, U + 0.5*dt*u2, W); k3x = Vv + 0.5*dt*a2
        a4, t4, u4, W, d4 = rhs(X + dt*k3x, Vv + dt*a3, VT + dt*t3, U + dt*u3, W); k4x = Vv + dt*a3
        a1[frz] = 0; a2[frz] = 0; a3[frz] = 0; a4[frz] = 0; k1x[frz] = 0; k2x[frz] = 0; k3x[frz] = 0; k4x[frz] = 0
        X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(a1 + 2*a2 + 2*a3 + a4)
        VT = VT + dt/6*(t1 + 2*t2 + 2*t3 + t4); U = U + dt/6*(u1 + 2*u2 + 2*u3 + u4); t += dt
        Dneg_c += dt*(d1['Dneg'] + 2*d2['Dneg'] + 2*d3['Dneg'] + d4['Dneg'])/6; Dpos_c += dt*(d1['Dpos'] + 2*d2['Dpos'] + 2*d3['Dpos'] + d4['Dpos'])/6
        if s % every == 0:
            th = np.arctan2(X[:, 1] - C, X[:, 0] - C); out = []
            for Rr in (sR, 0.41, 0.84, 1.25, 2.1):
                kk = np.where(np.abs(R0 - Rr) < 0.5*d)[0]
                if len(kk) < 10: continue
                o = kk[np.argsort(th0[kk])]; nn = len(o); dth_ = np.angle(np.exp(1j*(th[o] - th0[o]))); dth_ -= dth_.mean()
                Fm = np.abs(np.fft.rfft(dth_))*2/nn; Fm[:2] = 0; q = np.argmax(Fm)
                out.append('R%.2f m%d %.2e' % (Rr, q, R0[o].mean()*Fm[q]/d))
            K = (G_ - 1)*U/VT/(m/VT)**G_; ev = ~frz
            FLl = stat.get('FL'); Wl = stat.get('W', W)
            def zmax(lo, hi): return (np.array([abs(Wl[a] - Wl[j])/(dn*dn) for a in np.where(ev & (R0 >= lo) & (R0 < hi))[0] for j, A, e, dn, c, p_, q_ in FLl[a] if j != a]).max() if (DIVMODE != 'vor' and FLl is not None) else 0.0)
            dWf = np.array([zmax(0.0, 0.5), zmax(0.5, 2.0), zmax(2.0, 2.5)])
            dWg = np.array([abs(Wl[a] - Wl[j])/(dn*dn) for a in np.where(ev)[0] for j, A, e, dn, c, p_, q_ in FLl[a] if j != a]) if (DIVMODE != 'vor' and FLl is not None) else np.zeros(1)
            print('t %5.1f | %s | max|dW|/d2 zones(<0.5, disc, 2-2.5) %.3f %.3f %.3f (global %.3f) sat %.3f | dE/E0 %+.2e | Dneg %.2e | K rms %.1e | newton %d res %.1e loose %d  [%.0fs]' % (
                t, ' | '.join(out), dWf[0], dWf[1], dWf[2], dWg.max(), np.mean(dWg > 0.9), (Etot(X, Vv, U) - E0)/abs(E0), Dneg_c/abs(E0), np.sqrt(np.mean((K[ev]/K0[ev] - 1)**2)), stat['newton'], stat['resid'], stat.get('loose', 0), time.time() - t0), flush=True)
            stat['newton'] = 0
            np.save('/gpfs/kjhan/LagForce/PolyBnd/daeP_%s.npy' % tag, dict(X=X, V=Vv, VT=VT, U=U, W=W, R0=R0), allow_pickle=True)
