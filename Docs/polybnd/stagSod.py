#!/usr/bin/env python3
"""2D Sod shock tube with the staggered (zone-centred) formulation of stagP.py + zone artificial viscosity.
Zones = Delaunay triangles of the particles; m_T, K_T per zone; P_T = K_T (m_T/A_T)^gamma;
F_a = sum_T (P_T + Q_T) dA_T/dx_a ; dU_T/dt = -(P_T+Q_T) dA_T/dt  (exact, since U_T = K_T m_T^g A_T^(1-g)/(g-1) and K_T
is advanced by the Q heating only: dK_T/dt = (g-1) (-Q_T dA_T/dt) / (m_T^g A_T^(1-g)) ).
Q_T (Wilkins/Caramana form, compression only, affine-limited):
   div_T = (dA_T/dt)/A_T (exact for the linear velocity field inside the triangle)
   Q_T = rho_T * psi_T * (c2 * l_T^2 * div_T^2 + c1 * c_s,T * l_T * |div_T|)  if div_T < 0 else 0,  l_T = sqrt(A_T)
   psi_T = Caramana-type limiter: 1 - min over neighbour zones of (minmod(div_T, div_N)/div_T) ... (0 for an affine field)
Usage: stagSod.py n dt T tag [c1 c2 psi(0/1) remap(0/1)]
Box [0,2]x[0,1] periodic in y, reflecting walls in x realised by frozen particles at x<0.1 and x>1.9.
Left state rho=1,P=1 ; right rho=0.125,P=0.1 ; interface at x=1 ; gamma=1.4 ; equal-mass particles.
"""
import sys, time, numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely.strtree import STRtree
n = int(sys.argv[1]); dt0 = float(sys.argv[2]); T = float(sys.argv[3]); tag = sys.argv[4]
c1 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5; c2 = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
PSI = int(sys.argv[7]) if len(sys.argv) > 7 else 1; REMAP = int(sys.argv[8]) if len(sys.argv) > 8 else 1
G_ = 1.4; Lx, Ly = 2.0, 1.0
# ---- particles: equal mass, hexagonal-ish lattice, spacing d on the left, d*sqrt(8) on the right (rho ratio 8)
def lattice(x0, x1, dd, ny):
    dy = Ly/ny; nx = int(round((x1 - x0)/(dd*np.sqrt(3)/2))); dx = (x1 - x0)/nx
    P = []
    for i in range(nx):
        for j in range(ny):
            P.append((x0 + (i + 0.5)*dx, (j + 0.5 + 0.5*(i % 2))*dy))
    return np.array(P)
d = Ly/n; XL = lattice(0.0, 1.0, d, n); nyR = max(4, int(round(n/np.sqrt(8)))); XR = lattice(1.0, 2.0, Ly/nyR, nyR)
X = np.vstack([XL, XR]); N = len(X); rng = np.random.default_rng(1); X += 1e-3*d*rng.standard_normal(X.shape)
frozen = (X[:, 0] < 0.1) | (X[:, 0] > 1.9)
def rho0(x): return np.where(x < 1.0, 1.0, 0.125)
def P_0(x): return np.where(x < 1.0, 1.0, 0.1)
# periodic images in y for the triangulation
def triangulate(X):
    Xi = np.vstack([X, X + [0, Ly], X - [0, Ly]]); src = np.concatenate([np.arange(N)]*3); sh = np.concatenate([np.zeros(N, int), np.ones(N, int), -np.ones(N, int)])
    tri = Delaunay(Xi); zones = {}
    for simp in tri.simplices:
        c = Xi[simp].mean(0)
        if not (0 <= c[1] < Ly): continue
        if frozen[src[simp]].all(): continue   # triangles inside the static wall layers are not zones
        key = tuple(sorted((int(src[i]), int(sh[i])) for i in simp)); P = np.array([Xi[i] for i in simp]); P = np.array([X[i_] + [0, Ly*s_] for i_, s_ in key])
        e1 = P[1] - P[0]; e2 = P[2] - P[0]
        if e1[0]*e2[1] - e1[1]*e2[0] < 0: key = (key[0], key[2], key[1])
        zones[key] = None   # orientation fixed at creation (CCW)
    return zones
def zone_geom(key, X):
    ids = [k_[0] for k_ in key]; P = np.array([X[i] + [0, Ly*s_] for i, s_ in key])
    e1 = P[1] - P[0]; e2 = P[2] - P[0]; A2 = e1[0]*e2[1] - e1[1]*e2[0]
    dA = 0.5*np.array([[P[1][1] - P[2][1], P[2][0] - P[1][0]], [P[2][1] - P[0][1], P[0][0] - P[2][0]], [P[0][1] - P[1][1], P[1][0] - P[0][0]]])
    return ids, P, 0.5*A2, dA
def all_geom(keys, X):
    ids = np.zeros((len(keys), 3), int); areas = np.zeros(len(keys)); dAs = np.zeros((len(keys), 3, 2)); polys = []
    for q, key in enumerate(keys):
        i_, P, A, dA = zone_geom(key, X); ids[q] = i_; areas[q] = A; dAs[q] = dA; polys.append(P)
    return ids, areas, dAs, polys
keys = list(triangulate(X).keys()); ids, areas, dAs, polys = all_geom(keys, X)
cen = np.array([p.mean(0) for p in polys]); mT = rho0(cen[:, 0])*areas; KT = P_0(cen[:, 0])/(mT/areas)**G_
mA = np.zeros(N); np.add.at(mA, ids.ravel(), np.repeat(mT/3, 3))
zstate = {k_: (mT[q], KT[q]) for q, k_ in enumerate(keys)}
Vv = np.zeros((N, 2))
# zone neighbours (sharing an edge) for the limiter, rebuilt with the triangulation
def neighbours(keys):
    edge = {}
    for q, k_ in enumerate(keys):
        for a_ in range(3):
            u_, w_ = sorted((k_[a_], k_[(a_ + 1) % 3])); e = (u_[0], w_[0], tuple(np.subtract(w_[1], u_[1]).ravel())); edge.setdefault(e, []).append(q)
    nb = [[] for _ in keys]
    for e, qs in edge.items():
        if len(qs) == 2: nb[qs[0]].append(qs[1]); nb[qs[1]].append(qs[0])
    return nb
def forcesQ(ids, areas, dAs, mTz, KTz, Vs, nb):
    rho = mTz/areas; P = KTz*rho**G_; cs = np.sqrt(G_*P/rho)
    dAdt = np.einsum('qvk,qvk->q', dAs, Vs[ids]); div = dAdt/areas; l = np.sqrt(areas)
    Q = np.zeros_like(P); comp = div < 0
    if PSI:
        psi = np.ones_like(P)
        for q in np.where(comp)[0]:
            r = [min(max(div[m]/div[q], 0.0), 1.0) if div[m] < 0 else 0.0 for m in nb[q]]   # minmod ratio of neighbour compression
            psi[q] = 1.0 - min(r) if r else 1.0
    else: psi = np.ones_like(P)
    Q[comp] = rho[comp]*psi[comp]*(c2*l[comp]**2*div[comp]**2 + c1*cs[comp]*l[comp]*np.abs(div[comp]))
    F = np.zeros((N, 2))
    for v in range(3): np.add.at(F, ids[:, v], (P + Q)[:, None]*dAs[:, v])
    heat = -Q*dAdt   # dU_T/dt from viscosity (>=0)
    return F, heat, Q
def energy(Vv, areas, mTz, KTz): return 0.5*(mA*(Vv*Vv).sum(1)).sum() + (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1)).sum()
def step(X, Vv, keys, zstate, dt, nb):
    mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
    def rhs(Xs, Vs, Ks):
        ids_, areas_, dAs_, _ = all_geom(keys, Xs)
        if (areas_ <= 0).any(): raise RuntimeError('inverted zone')
        F, heat, Q = forcesQ(ids_, areas_, dAs_, mTz, Ks, Vs, nb)
        a = F/mA[:, None]; a[frozen] = 0
        dK = (G_ - 1)*heat/(mTz**G_*areas_**(1 - G_))
        return Vs*(~frozen)[:, None], a, dK
    k1 = rhs(X, Vv, KTz); k2 = rhs(X + 0.5*dt*k1[0], Vv + 0.5*dt*k1[1], KTz + 0.5*dt*k1[2])
    k3 = rhs(X + 0.5*dt*k2[0], Vv + 0.5*dt*k2[1], KTz + 0.5*dt*k2[2]); k4 = rhs(X + dt*k3[0], Vv + dt*k3[1], KTz + dt*k3[2])
    X = X + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); Vv = Vv + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]); KTz = KTz + dt/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
    zs = {k_: (mTz[q], KTz[q]) for q, k_ in enumerate(keys)}
    ids, areas, dAs, polys = all_geom(keys, X)
    if (areas <= 0).any(): raise RuntimeError('inverted zone at step end')
    Eb = energy(Vv, areas, mTz, KTz); dS = 0.0; nf = 0
    if REMAP:
        newz = triangulate(X); newkeys = list(newz.keys()); fresh = [k_ for k_ in newkeys if k_ not in zs]
        if fresh:
            oldpolys = [Polygon(p) for p in polys]; tree = STRtree(oldpolys); rho_old = mTz/areas; e_old = (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1))/areas
            gone = [k_ for k_ in keys if k_ not in newz]; S_old = sum(zs[k_][0]*np.log(zs[k_][1]) for k_ in gone); m_old = sum(zs[k_][0] for k_ in gone)
            S_new = 0.0; m_sum = 0.0
            for k_ in fresh:
                _, Pn, An, _ = zone_geom(k_, X); pn = Polygon(Pn); m_new = 0.0; u_new = 0.0
                for sy in (-1, 0, 1):
                    pq = Polygon(Pn + [0, sy*Ly])
                    for qi in tree.query(pq):
                        ov = pq.intersection(oldpolys[qi]).area
                        if ov > 0: m_new += ov*rho_old[qi]; u_new += ov*e_old[qi]
                if m_new <= 0: raise RuntimeError('remap: no overlap')
                K_new = (G_ - 1)*u_new/An/(m_new/An)**G_; zs[k_] = (m_new, K_new); S_new += m_new*np.log(K_new); m_sum += m_new
            if abs(m_sum - m_old) > 1e-8*max(m_old, 1e-30): raise RuntimeError('remap: mass mismatch %.3e vs %.3e' % (m_sum, m_old))
            dS = S_new - S_old; nf = len(fresh)
            for k_ in gone: del zs[k_]
        keys = newkeys
    return X, Vv, keys, zs, Eb, nf, dS
def exact_sod(x, t):
    # standard Sod exact solution (gamma=1.4), interface at x=1
    gl, gr = (1.0, 1.0, 0.0), (0.125, 0.1, 0.0); g = G_
    rl, pl, ul = gl; rr, pr, ur = gr; cl = np.sqrt(g*pl/rl); cr = np.sqrt(g*pr/rr)
    def f(p, r, pk, ck):
        if p > pk: return (p - pk)*np.sqrt(2/((g + 1)*r)/(p + (g - 1)/(g + 1)*pk))
        return 2*ck/(g - 1)*((p/pk)**((g - 1)/(2*g)) - 1)
    from scipy.optimize import brentq
    ps = brentq(lambda p: f(p, rl, pl, cl) + f(p, rr, pr, cr) + ur - ul, 1e-6, 1.0); us = 0.5*(ul + ur) + 0.5*(f(ps, rr, pr, cr) - f(ps, rl, pl, cl))
    rsl = rl*(ps/pl)**(1/g); rsr = rr*(ps/pr + (g - 1)/(g + 1))/((g - 1)/(g + 1)*ps/pr + 1)
    S = ur + cr*np.sqrt((g + 1)/(2*g)*ps/pr + (g - 1)/(2*g)); csl = np.sqrt(g*ps/rsl)
    xi = (x - 1.0)/t; rho = np.empty_like(x); P = np.empty_like(x); u = np.empty_like(x)
    r1 = xi < ul - cl; r2 = (xi >= ul - cl) & (xi < us - csl); r3 = (xi >= us - csl) & (xi < us); r4 = (xi >= us) & (xi < S); r5 = xi >= S
    rho[r1] = rl; P[r1] = pl; u[r1] = ul
    uf = 2/(g + 1)*(cl + (g - 1)/2*ul + xi[r2]); cf = cl - (g - 1)/2*(uf - ul); rho[r2] = rl*(cf/cl)**(2/(g - 1)); P[r2] = pl*(cf/cl)**(2*g/(g - 1)); u[r2] = uf
    rho[r3] = rsl; P[r3] = ps; u[r3] = us; rho[r4] = rsr; P[r4] = ps; u[r4] = us; rho[r5] = rr; P[r5] = pr; u[r5] = ur
    return rho, P, u
nb = neighbours(keys); mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
E0 = energy(Vv, areas, mTz, KTz); S_tot = 0.0; nrej = 0; nremap = 0; Eremap = 0.0; t = 0.0; t0 = time.time(); s = 0
print('stagSod n %d N %d zones %d c1 %.2f c2 %.2f psi %d remap %d E0 %.6e' % (n, N, len(keys), c1, c2, PSI, REMAP, E0), flush=True)
every = max(1, int(round(0.05/dt0)))
while t < T - 1e-12:
    s += 1; dtn = min(dt0, T - t); dts = dtn; nsub = 1; done = 0; saved = (X.copy(), Vv.copy(), list(keys), dict(zstate), S_tot, nremap, Eremap); ok = True
    while done < nsub:
        try:
            X, Vv, keys, zstate, Eb, nf, dS = step(X, Vv, keys, zstate, dts, nb); done += 1; nremap += nf; S_tot += dS
            if nf: nb = neighbours(keys)
            ids, areas, dAs, polys = all_geom(keys, X); mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
            Eremap += energy(Vv, areas, mTz, KTz) - Eb
        except RuntimeError as ex:
            nrej += 1
            if nsub >= 16: print('ABORT t %.3f: %s' % (t, ex), flush=True); ok = False; break
            X, Vv, keys, zstate, S_tot, nremap, Eremap = saved; nb = neighbours(keys); nsub *= 2; dts = dtn/nsub; done = 0
    if not ok: break
    t += dtn
    if s % every == 0:
        Ea = energy(Vv, areas, mTz, KTz); cen = np.array([p.mean(0) for p in polys]); rho = mTz/areas; P = KTz*rho**G_
        inner = (cen[:, 0] > 0.15) & (cen[:, 0] < 1.85); re, pe, ue = exact_sod(cen[inner, 0], t)
        L1r = np.mean(np.abs(rho[inner] - re)); L1p = np.mean(np.abs(P[inner] - pe))
        vz = Vv[ids].mean(1)[:, 0]; L1u = np.mean(np.abs(vz[inner] - ue))
        print('t %.3f | L1 rho %.4f P %.4f u %.4f | dE/E0 %+.2e (remap %+.2e) | S_remap %+.3e | <K> jump: max K/K0 %.3f | remapped %d rej %d | zones %d  [%.0fs]' % (
            t, L1r, L1p, L1u, (Ea - E0)/abs(E0), Eremap/abs(E0), S_tot, (KTz/KTz.min()).max(), nremap, nrej, len(keys), time.time() - t0), flush=True)
# final profile dump
cen = np.array([p.mean(0) for p in polys]); rho = mTz/areas; P = KTz*rho**G_; vz = Vv[ids].mean(1)[:, 0]
np.savetxt(tag + '_prof.txt', np.c_[cen[:, 0], cen[:, 1], rho, P, vz, areas], header='x y rho P u area')
