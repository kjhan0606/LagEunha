"""dae.py: minimal consistent Laguerre-DAE on the ring IC (particle = generator). Approved 2026-10-08.
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
from scipy.spatial import Delaunay
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from shapely.geometry import Polygon
import shapely
G_ = 5/3; P0 = 1e-6; EPS = 0.01; RIN, ROUT = 0.5, 2.0
def J(v): return np.array([-v[1], v[0]])
def cell_job(args):
    idx, XA, WA, NBs, Rc, nseg = args; out = []
    circ0 = np.array([[np.cos(th), np.sin(th)] for th in np.linspace(0, 2*np.pi, nseg, endpoint=False)])*Rc
    rmin = Rc*np.cos(np.pi/nseg) - 1e-9
    for a, xa, wa, (nbid, XN, WN) in zip(idx, XA, WA, NBs):
        B = 8*Rc + 1.0
        cell = Polygon([(xa[0]-B, xa[1]-B), (xa[0]+B, xa[1]-B), (xa[0]+B, xa[1]+B), (xa[0]-B, xa[1]+B)])
        planes = []
        for j, xj, wj in zip(nbid, XN, WN):
            dv = xj - xa; dn = np.linalg.norm(dv); e = dv/dn; s = 0.5*dn + (wa - wj)/(2*dn); pt = xa + s*e
            t = J(e)*50; nrm = e*50
            cell = cell.intersection(Polygon([pt + t, pt - t, pt - t - nrm, pt + t - nrm])); planes.append((j, e, dn, pt, xj))
        cell = cell.intersection(Polygon(xa + circ0))
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
                if not done:
                    if np.linalg.norm(p - xa) >= rmin and np.linalg.norm(q - xa) >= rmin and np.linalg.norm(0.5*(p + q) - xa) >= rmin:
                        g[a] = g.get(a, 0) - J(q - p); fl.append((a, L, -J(q - p)/L, 0.0, 0.5*(p + q), p, q))
                    else: raise RuntimeError('unclassified segment, cell %d' % a)
        out.append((a, V, g, fl))
    return out
if __name__ == '__main__':
    n = int(sys.argv[1]); dt = float(sys.argv[2]); T = float(sys.argv[3]); nproc = int(sys.argv[4]); nseg = int(sys.argv[5])
    P0 = P0*float(sys.argv[6]); sring = int(sys.argv[7]); smode = int(sys.argv[8]); samp = float(sys.argv[9]); tag = sys.argv[10]
    EDGEW0 = int(sys.argv[11]) if len(sys.argv) > 11 else 0; DIVMODE = sys.argv[12] if len(sys.argv) > 12 else 'face'; DBG = len(sys.argv) > 13 and sys.argv[13] == 'dbg'; TAUV = float(sys.argv[14]) if len(sys.argv) > 14 else 0.0
    d = 6.0/n; Rc = d
    P = []; lab = []
    for k in range(int(round((ROUT - RIN)/d))):
        r = RIN + (k + 0.5)*d; m_ = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(m_) + 0.5)*2*np.pi/m_
        P.append(np.stack([r*np.cos(th), r*np.sin(th)], 1)); lab += [k]*m_
    X = np.concatenate(P); lab = np.array(lab); N = len(X)
    pool = Pool(nproc)
    def geom(X, W):
        tri = Delaunay(X); ip, ix = tri.vertex_neighbor_vertices
        # 2-hop candidates (power-diagram topology may differ from Delaunay)
        nb = []
        for a in range(N):
            n1 = ix[ip[a]:ip[a+1]]; s = set(n1.tolist())
            for j in n1: s.update(ix[ip[j]:ip[j+1]].tolist())
            s.discard(a); nb.append(np.fromiter(s, dtype=np.int64))
        geom.nb = nb
        chunks = np.array_split(np.arange(N), nproc*4); jobs = []
        for ch in chunks:
            jobs.append((ch, X[ch], W[ch], [(nb[a], X[nb[a]], W[nb[a]]) for a in ch], Rc, nseg))
        V = np.zeros(N); Gs = [None]*N; FL = [None]*N
        for res in pool.map(cell_job, jobs):
            for a, v, g, fl in res: V[a] = v; Gs[a] = g; FL[a] = fl
        return V, Gs, FL
    def arcflag(FL): return np.array([any(j == a for j, A, e, dn, c, p_, q_ in FL[a]) for a in range(N)])
    def lsq_coef(X, Vv):
        # per cell: c (5,2): [Gx, Gy, Hxx, Hxy, Hyy] for (vx, vy); v~(x+y) = v + c0 yx + c1 yy + 1/2 c2 yx^2 + c3 yx yy + 1/2 c4 yy^2
        nb = geom.nb; C = np.zeros((N, 5, 2))
        for a in range(N):
            ys = X[nb[a]] - X[a]; w = 1.0/(ys*ys).sum(1); sw = np.sqrt(w)
            A = np.stack([ys[:, 0], ys[:, 1], 0.5*ys[:, 0]**2, ys[:, 0]*ys[:, 1], 0.5*ys[:, 1]**2], 1)*sw[:, None]
            b = (Vv[nb[a]] - Vv[a])*sw[:, None]
            C[a], *_ = np.linalg.lstsq(A, b, rcond=None)
        return C
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
            Lw = (Lw.tocsr() + (0*Pin if EDGEW0 else Pin)).tocsc(); W = W + spsolve(Lw, -res0); W[arc] = 0
            if not EDGEW0: W -= W.mean()
        if rn < 1e-5:   # arc-polygonisation floor: accept and count
            stat['loose'] = stat.get('loose', 0) + 1; stat['resid'] = rn; stat['meanres'] = mres/VT.mean()
            if EDGEW0: VT = VT.copy(); VT[arc] = V[arc]
            return W, V, Gs, FL, VT, arc
        raise RuntimeError('projection did not converge: %.2e' % rn)
    def rhs(X, Vv, VT, U, W):
        W, V, Gs, FL, VT, arc = project(X, W, VT)
        Pc = (G_ - 1)*U/VT; F = np.zeros_like(X); dVgeo = np.zeros(N); dVT = np.zeros(N)
        for a in range(N):
            for k, g in Gs[a].items(): F[k] += Pc[a]*g; dVgeo[a] += g @ Vv[k]
            for j, A, e, dn, c, p_, q_ in FL[a]: dVT[a] += A*0.5*((Vv[a] + Vv[j]) @ e)
        if DIVMODE == 'lsq': dVT = VT*lsq_div(X, Vv)
        if DIVMODE == 'flux3':
            C = lsq_coef3(X, Vv); dVT = np.zeros(N); gq = 0.5/np.sqrt(3.0)
            for a in range(N):
                for j, A, e, dn, c, p_, q_ in FL[a]:
                    f = 0.0
                    for s_ in (-gq, gq):
                        pt = c + s_*(q_ - p_)
                        if j == a: u = vrec3(a, pt - X[a], Vv, C)
                        else: u = 0.5*(vrec3(a, pt - X[a], Vv, C) + vrec3(j, pt - X[j], Vv, C))
                        f += 0.5*(u @ e)
                    dVT[a] += A*f
        if DIVMODE == 'flux2':
            C = lsq_coef(X, Vv); dVT = np.zeros(N)
            for a in range(N):
                for j, A, e, dn, c, p_, q_ in FL[a]:
                    if j == a: u = vrec(a, c - X[a], Vv, C)
                    else: u = 0.5*(vrec(a, c - X[a], Vv, C) + vrec(j, c - X[j], Vv, C))
                    dVT[a] += A*(u @ e)
        if TAUV > 0:
            Vvor, _, _ = geom(X, np.zeros(N))   # relaxation target = W=0 Voronoi volume (V(x,W) == VT by construction)
            dVT = dVT + (Vvor - VT)/TAUV
        if DBG: print('   dVT*dt/V: rms %.2e max %.2e (ring of max: %d)  dVgeo*dt/V rms %.2e' % (np.sqrt(np.mean((dVT*dt/V)**2)), np.abs(dVT*dt/V).max(), lab[np.argmax(np.abs(dVT/V))], np.sqrt(np.mean((dVgeo*dt/V)**2))), flush=True)
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
        r2 = (X*X).sum(1) + EPS*EPS
        acc = F/m[:, None] - X/r2[:, None]**1.5
        dW = np.array([abs(W[a] - W[j])/(dn*dn) for a in range(N) for j, A, e, dn, c, p_, q_ in FL[a] if j != a])
        return acc, dVT, dU, W, dict(maxdW=dW.max(), sat=np.mean(dW > 0.9), Dneg=Dneg, Dpos=Dpos, dVmis=np.sqrt(np.mean((dVT - dVgeo)**2))*dt/V.mean(), VT=VT)
    W = np.zeros(N); V0, G0, FL0 = geom(X, W); m = V0.copy(); VT = V0.copy(); U = P0*VT/(G_ - 1); K0 = P0/(m/VT)**G_
    th0 = np.arctan2(X[:, 1], X[:, 0]); R0 = np.hypot(*X.T)
    sel = np.where(lab == sring)[0]; dth = samp*d/R0[sel]*np.cos(smode*(th0[sel] - th0[sel][0]))
    X[sel] = np.stack([R0[sel]*np.cos(th0[sel] + dth), R0[sel]*np.sin(th0[sel] + dth)], 1)
    R = np.hypot(*X.T); Om = (R*R + EPS*EPS)**-0.75; Vv = np.stack([-Om*X[:, 1], Om*X[:, 0]], 1)
    def Etot(X, Vv, U):
        r2 = (X*X).sum(1) + EPS*EPS
        return 0.5*(m*(Vv*Vv).sum(1)).sum() + U.sum() - (m/np.sqrt(r2)).sum()
    E0 = Etot(X, Vv, U)
    print('DAE n %d N %d ring%d N %d P0 %.1e cs/d %.3f seed m %d amp %.3f d  E0 %.6e edgew0 %d divmode %s tauv %.2f' % (n, N, sring, len(sel), P0, np.sqrt(G_*P0)/d, smode, samp, E0, EDGEW0, DIVMODE, TAUV), flush=True)
    t = 0.0; t0 = time.time(); every = int(round(1.0/dt)); Dneg_c = 0.0; Dpos_c = 0.0
    for s in range(1, int(round(T/dt)) + 1):
        a1, t1, u1, W, d1 = rhs(X, Vv, VT, U, W); VT = d1['VT']; k1x = Vv
        a2, t2, u2, W, d2 = rhs(X + 0.5*dt*k1x, Vv + 0.5*dt*a1, VT + 0.5*dt*t1, U + 0.5*dt*u1, W); k2x = Vv + 0.5*dt*a1
        a3, t3, u3, W, d3 = rhs(X + 0.5*dt*k2x, Vv + 0.5*dt*a2, VT + 0.5*dt*t2, U + 0.5*dt*u2, W); k3x = Vv + 0.5*dt*a2
        a4, t4, u4, W, d4 = rhs(X + dt*k3x, Vv + dt*a3, VT + dt*t3, U + dt*u3, W); k4x = Vv + dt*a3
        X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(a1 + 2*a2 + 2*a3 + a4)
        VT = VT + dt/6*(t1 + 2*t2 + 2*t3 + t4); U = U + dt/6*(u1 + 2*u2 + 2*u3 + u4); t += dt
        Dneg_c += dt*(d1['Dneg'] + 2*d2['Dneg'] + 2*d3['Dneg'] + d4['Dneg'])/6; Dpos_c += dt*(d1['Dpos'] + 2*d2['Dpos'] + 2*d3['Dpos'] + d4['Dpos'])/6
        if s % every == 0:
            th = np.arctan2(X[:, 1], X[:, 0]); o = sel[np.argsort(th0[sel])]; nn = len(o)
            dth = np.angle(np.exp(1j*(th[o] - th0[o]))); dth -= dth.mean(); Fm = np.abs(np.fft.rfft(dth))*2/nn; Fm[:2] = 0; q = np.argmax(Fm)
            K = (G_ - 1)*U/VT/(m/VT)**G_
            print('t %5.1f | seed m%d %.2e max m%d %.2e | max|dW|/d2 %.3f sat %.3f | dE/E0 %+.2e | Dneg %.2e Dpos %.2e | K rms %.1e | VT/Vvor mis/step %.1e | newton %d res %.1e mean %.1e loose %d  [%.0fs]' % (
                t, smode, np.hypot(*X[o].T).mean()*Fm[smode]/d, q, np.hypot(*X[o].T).mean()*Fm[q]/d, d4['maxdW'], d4['sat'], (Etot(X, Vv, U) - E0)/abs(E0),
                Dneg_c/abs(E0), Dpos_c/abs(E0), np.sqrt(np.mean((K/K0 - 1)**2)), d4['dVmis'], stat['newton'], stat['resid'], stat['meanres'], stat.get('loose', 0), time.time() - t0), flush=True)
            stat['newton'] = 0
            np.save('/gpfs/kjhan/LagForce/PolyBnd/dae_%s.npy' % tag, dict(X=X, V=Vv, VT=VT, U=U, W=W, lab=lab), allow_pickle=True)
