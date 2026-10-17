#!/usr/bin/env python3
"""Gresho vortex with the staggered (zone-centred) formulation + zone Q of stagSod.py.
Periodic box [0,1]^2, gamma=5/3, rho=1, standard Gresho profile (P0=5); equal-mass particles on a jittered hex lattice.
Measures: L1 of v_phi vs exact, kinetic energy decay, viscous heating ledger (sum of Q work), remap entropy.
Usage: stagGresho.py n dt T tag [c1 c2 psi(0/1)]
"""
import sys, time, numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely.strtree import STRtree
n = int(sys.argv[1]); dt0 = float(sys.argv[2]); T = float(sys.argv[3]); tag = sys.argv[4]
c1 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5; c2 = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
PSI = int(sys.argv[7]) if len(sys.argv) > 7 else 1
G_ = 5.0/3.0; L = 1.0; d = L/n
X = np.array([((i + 0.5 + 0.25*(j % 2)*0)*d, (j + 0.5)*d) for i in range(n) for j in range(n)]); N = len(X)
X[:, 0] += 0.5*d*(np.arange(N) // n % 2)*0   # square lattice + jitter (hex offset breaks periodicity for odd n)
rng = np.random.default_rng(1); X += 0.05*d*rng.standard_normal(X.shape); X = np.mod(X, L)
def gresho(r):
    v = np.where(r < 0.2, 5*r, np.where(r < 0.4, 2 - 5*r, 0.0))
    P = np.where(r < 0.2, 5 + 12.5*r*r, np.where(r < 0.4, 9 + 12.5*r*r - 20*r + 4*np.log(r/0.2), 3 + 4*np.log(2)))
    return v, P
def triangulate(X):
    Xi = np.vstack([X + [sx*L, sy*L] for sx in (-1, 0, 1) for sy in (-1, 0, 1)]); src = np.tile(np.arange(N), 9)
    sh = np.repeat(np.array([(sx, sy) for sx in (-1, 0, 1) for sy in (-1, 0, 1)]), N, axis=0)
    tri = Delaunay(Xi); zones = {}
    for simp in tri.simplices:
        c = Xi[simp].mean(0)
        if not (0 <= c[0] < L and 0 <= c[1] < L): continue
        key = tuple(sorted((int(src[i]), (int(sh[i][0]), int(sh[i][1]))) for i in simp)); P = np.array([X[i_] + L*np.array(s_) for i_, s_ in key])
        e1 = P[1] - P[0]; e2 = P[2] - P[0]
        if e1[0]*e2[1] - e1[1]*e2[0] < 0: key = (key[0], key[2], key[1])
        zones[canon(key)] = None   # orientation fixed at creation (CCW), canonical cyclic order
    return zones
def canon(key):
    i = min(range(3), key=lambda q: key[q]); return (key[i], key[(i + 1) % 3], key[(i + 2) % 3])
def zone_geom(key, X):
    ids = [k_[0] for k_ in key]; P = np.array([X[i] + L*np.array(s_) for i, s_ in key])
    e1 = P[1] - P[0]; e2 = P[2] - P[0]; A2 = e1[0]*e2[1] - e1[1]*e2[0]
    dA = 0.5*np.array([[P[1][1] - P[2][1], P[2][0] - P[1][0]], [P[2][1] - P[0][1], P[0][0] - P[2][0]], [P[0][1] - P[1][1], P[1][0] - P[0][0]]])
    return ids, P, 0.5*A2, dA
def all_geom(keys, X):
    ids = np.zeros((len(keys), 3), int); areas = np.zeros(len(keys)); dAs = np.zeros((len(keys), 3, 2)); polys = []
    for q, key in enumerate(keys):
        i_, P, A, dA = zone_geom(key, X); ids[q] = i_; areas[q] = A; dAs[q] = dA; polys.append(P)
    return ids, areas, dAs, polys
keys = list(triangulate(X).keys()); ids, areas, dAs, polys = all_geom(keys, X)
cen = np.array([p.mean(0) for p in polys]); rc = np.hypot(cen[:, 0] - 0.5, cen[:, 1] - 0.5); _, Pz0 = gresho(rc)
mT = areas*1.0; KT = Pz0/(mT/areas)**G_
mA = np.zeros(N); np.add.at(mA, ids.ravel(), np.repeat(mT/3, 3))
zstate = {k_: (mT[q], KT[q]) for q, k_ in enumerate(keys)}
rp = np.hypot(X[:, 0] - 0.5, X[:, 1] - 0.5); vp, _ = gresho(rp); Vv = np.c_[-vp*(X[:, 1] - 0.5)/np.maximum(rp, 1e-12), vp*(X[:, 0] - 0.5)/np.maximum(rp, 1e-12)]
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
    if PSI and comp.any():
        psi = np.ones_like(P)
        for q in np.where(comp)[0]:
            r = [min(max(div[m]/div[q], 0.0), 1.0) if div[m] < 0 else 0.0 for m in nb[q]]
            psi[q] = 1.0 - min(r) if r else 1.0
    else: psi = np.ones_like(P)
    Q[comp] = rho[comp]*psi[comp]*(c2*l[comp]**2*div[comp]**2 + c1*cs[comp]*l[comp]*np.abs(div[comp]))
    F = np.zeros((N, 2))
    for v in range(3): np.add.at(F, ids[:, v], (P + Q)[:, None]*dAs[:, v])
    return F, -Q*dAdt
def energy(Vv, areas, mTz, KTz): return 0.5*(mA*(Vv*Vv).sum(1)).sum() + (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1)).sum()
Qwork = 0.0
def step(X, Vv, keys, zstate, dt, nb):
    global Qwork
    mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
    def rhs(Xs, Vs, Ks):
        ids_, areas_, dAs_, _ = all_geom(keys, Xs)
        if (areas_ <= 0).any(): raise RuntimeError('inverted zone')
        F, heat = forcesQ(ids_, areas_, dAs_, mTz, Ks, Vs, nb)
        return Vs, F/mA[:, None], (G_ - 1)*heat/(mTz**G_*areas_**(1 - G_)), heat.sum()
    k1 = rhs(X, Vv, KTz); k2 = rhs(X + 0.5*dt*k1[0], Vv + 0.5*dt*k1[1], KTz + 0.5*dt*k1[2])
    k3 = rhs(X + 0.5*dt*k2[0], Vv + 0.5*dt*k2[1], KTz + 0.5*dt*k2[2]); k4 = rhs(X + dt*k3[0], Vv + dt*k3[1], KTz + dt*k3[2])
    X = X + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); Vv = Vv + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]); KTz = KTz + dt/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
    Qwork += dt/6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3])
    zs = {k_: (mTz[q], KTz[q]) for q, k_ in enumerate(keys)}
    ids, areas, dAs, polys = all_geom(keys, X)
    if (areas <= 0).any(): raise RuntimeError('inverted zone at step end')
    Eb = energy(Vv, areas, mTz, KTz); dS = 0.0; nf = 0
    Xw = np.mod(X, L); shift = np.round((X - Xw)/L).astype(int)   # wrap particles; zone keys shift accordingly
    newz = triangulate(Xw); newkeys = list(newz.keys())
    # re-express old keys in the wrapped frame
    zs2 = {}
    for k_, v in zs.items():
        kk = tuple((i, (s_[0] + shift[i][0], s_[1] + shift[i][1])) for i, s_ in k_)
        # canonical frame: centroid in the box
        P = np.array([Xw[i] + L*np.array(s_) for i, s_ in kk]); c = P.mean(0); off = np.floor(c/L).astype(int)
        kk = canon(tuple((i, (s_[0] - off[0], s_[1] - off[1])) for i, s_ in kk)); zs2[kk] = v
    zs = zs2; X = Xw; fresh = [k_ for k_ in newkeys if k_ not in zs]
    if fresh:
        oldkeys = list(zs.keys()); _, oareas, _, opolys = all_geom(oldkeys, X); om = np.array([zs[k_][0] for k_ in oldkeys]); oK = np.array([zs[k_][1] for k_ in oldkeys])
        oldpolys = [Polygon(p) for p in opolys]; tree = STRtree(oldpolys); rho_old = om/oareas; e_old = (oK*om**G_*oareas**(1 - G_)/(G_ - 1))/oareas
        gone = [k_ for k_ in oldkeys if k_ not in newz]; S_old = sum(zs[k_][0]*np.log(zs[k_][1]) for k_ in gone); m_old = sum(zs[k_][0] for k_ in gone)
        S_new = 0.0; m_sum = 0.0
        for k_ in fresh:
            _, Pn, An, _ = zone_geom(k_, X); m_new = 0.0; u_new = 0.0
            for sx in (-1, 0, 1):
                for sy in (-1, 0, 1):
                    pq = Polygon(Pn + [sx*L, sy*L])
                    for qi in tree.query(pq):
                        ov = pq.intersection(oldpolys[qi]).area
                        if ov > 0: m_new += ov*rho_old[qi]; u_new += ov*e_old[qi]
            if m_new <= 0: raise RuntimeError('remap: no overlap')
            K_new = (G_ - 1)*u_new/An/(m_new/An)**G_; zs[k_] = (m_new, K_new); S_new += m_new*np.log(K_new); m_sum += m_new
        if abs(m_sum - m_old) > 1e-8*max(m_old, 1e-30): raise RuntimeError('remap: mass mismatch %.3e vs %.3e' % (m_sum, m_old))
        dS = S_new - S_old; nf = len(fresh)
        for k_ in gone: del zs[k_]
    return X, Vv, newkeys, zs, Eb, nf, dS
nb = neighbours(keys); mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
E0 = energy(Vv, areas, mTz, KTz); KE0 = 0.5*(mA*(Vv*Vv).sum(1)).sum(); S_tot = 0.0; nrej = 0; nremap = 0; Eremap = 0.0; t = 0.0; t0 = time.time(); s = 0
print('stagGresho n %d N %d zones %d c1 %.2f c2 %.2f psi %d E0 %.6e KE0 %.6e' % (n, N, len(keys), c1, c2, PSI, E0, KE0), flush=True)
every = max(1, int(round(0.25/dt0)))
while t < T - 1e-12:
    s += 1; dtn = min(dt0, T - t); dts = dtn; nsub = 1; done = 0; saved = (X.copy(), Vv.copy(), list(keys), dict(zstate), Qwork, S_tot, nremap, Eremap); ok = True
    while done < nsub:
        try:
            X, Vv, keys, zstate, Eb, nf, dS = step(X, Vv, keys, zstate, dts, nb); done += 1; nremap += nf; S_tot += dS
            nb = neighbours(keys)
            ids, areas, dAs, polys = all_geom(keys, X); mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
            Eremap += energy(Vv, areas, mTz, KTz) - Eb
        except RuntimeError as ex:
            nrej += 1
            if nsub >= 16: print('ABORT t %.3f: %s' % (t, ex), flush=True); ok = False; break
            X, Vv, keys, zstate, Qwork, S_tot, nremap, Eremap = saved; nb = neighbours(keys); nsub *= 2; dts = dtn/nsub; done = 0
    if not ok: break
    t += dtn
    if s % every == 0:
        Ea = energy(Vv, areas, mTz, KTz); KE = 0.5*(mA*(Vv*Vv).sum(1)).sum()
        rp = np.hypot(X[:, 0] - 0.5, X[:, 1] - 0.5); vpe, _ = gresho(rp); vphi = (-Vv[:, 0]*(X[:, 1] - 0.5) + Vv[:, 1]*(X[:, 0] - 0.5))/np.maximum(rp, 1e-12)
        print('t %.2f | L1 vphi %.4f | KE/KE0 %.5f | Qwork/KE0 %.3e | dE/E0 %+.2e (remap %+.2e) | S_remap %+.3e | remapped %d rej %d  [%.0fs]' % (
            t, np.mean(np.abs(vphi - vpe)), KE/KE0, Qwork/KE0, (Ea - E0)/abs(E0), Eremap/abs(E0), S_tot, nremap, nrej, time.time() - t0), flush=True)
