"""stagShear.py: planar Couette shear box for the staggered (zone) formulation — Astra 27 §1 discriminating test.
Rows of particles (hex-like lattice, spacing d) in [0,Lx) x [0,Ly], periodic in x. Uniform rho=1, P0; velocity
v_x = S (y - Ly/2) (affine shear). Top and bottom rows are material boundaries in contact with a pressure reservoir
P_ext = P0 (as stagP2): F_a = sum_T (P_T - P0) dA_T/dx_a, E = KE + sum_T(u_T + P0 A_T). No gravity.
A longitudinal mode (displacement along x, wavenumber m on the row at y = y_seed) is seeded with amplitude amp*d.
Linear response: run with and without the seed (same everything else) and subtract (mode amplitude of the difference).
Diagnostics: seeded-row mode amplitude (seeded run minus unseeded run if REF given), remap counts, energy, S ledger.
args: nx ny dt T S P0 m amp tag   env: REMAP2 (0/1), REF (log of the unseeded run to subtract; optional), FIXCONN (0/1: never retriangulate; control)"""
import sys, os, time, numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely.strtree import STRtree
G_ = 5/3
nx = int(sys.argv[1]); ny = int(sys.argv[2]); dt = float(sys.argv[3]); T = float(sys.argv[4]); S_ = float(sys.argv[5]); P0 = float(sys.argv[6])
smode = int(sys.argv[7]); samp = float(sys.argv[8]); tag = sys.argv[9]
REMAP2 = int(os.environ.get('REMAP2', '0')); FIXCONN = int(os.environ.get('FIXCONN', '0'))
d = 1.0; Lx = nx*d; Ly = (ny - 1)*d*np.sqrt(3)/2 + 1e-9   # rows at y_j = j*d*sqrt(3)/2, half-offset alternate rows
X = np.array([(((i + 0.5*(j % 2)) % nx)*d, j*d*np.sqrt(3)/2) for j in range(ny) for i in range(nx)]); N = len(X)
row = np.repeat(np.arange(ny), nx); X0 = X.copy()
TOP = row == ny - 1; BOT = row == 0; OUTR = TOP | BOT
jseed = ny//2; sel = np.where(row == jseed)[0]
X[sel, 0] += samp*d*np.cos(2*np.pi*smode*X0[sel, 0]/Lx)
Vv = np.stack([S_*(X[:, 1] - Ly/2), np.zeros(N)], 1)
cs = np.sqrt(G_*P0); print('stagShear nx %d ny %d N %d S %.3f P0 %.3e c_s %.4f (c_s/d %.3f, c_s/(S d) %.3f) m %d amp %.3f REMAP2 %d FIXCONN %d' % (nx, ny, N, S_, P0, cs, cs/d, cs/(S_*d) if S_ else np.inf, smode, samp, REMAP2, FIXCONN), flush=True)
# ---- periodic-x Delaunay with image columns; zone keys (id, shift) oriented CCW at creation; exterior slivers (3 boundary-row vertices) excluded
def canon(key):
    i = min(range(3), key=lambda q: key[q]); return (key[i], key[(i + 1) % 3], key[(i + 2) % 3])
def triangulate(X):
    Xi = np.vstack([X + [s*Lx, 0] for s in (-1, 0, 1)]); src = np.tile(np.arange(N), 3); sh = np.repeat([-1, 0, 1], N)
    tri = Delaunay(Xi); zones = {}
    for simp in tri.simplices:
        c = Xi[simp].mean(0)
        if not (0 <= c[0] < Lx): continue
        i3 = [int(src[i]) for i in simp]
        if OUTR[i3].all() and (row[i3] == row[i3[0]]).all(): continue   # sliver along a boundary row: exterior
        key = [(int(src[i]), int(sh[i])) for i in simp]; P = np.array([X[i_] + [s_*Lx, 0] for i_, s_ in key])
        if (P[1][0] - P[0][0])*(P[2][1] - P[0][1]) - (P[2][0] - P[0][0])*(P[1][1] - P[0][1]) < 0: key = [key[0], key[2], key[1]]
        zones[canon(tuple(key))] = None
    return zones
def zone_geom(key, X):
    a, b, c = [X[i] + [s_*Lx, 0] for i, s_ in key]; A = 0.5*((b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1]))
    def J(v): return np.array([-v[1], v[0]])
    return np.array([a, b, c]), A, [0.5*J(c - b), 0.5*J(a - c), 0.5*J(b - a)]
def all_geom(keys, X):
    ids = np.array([[k_[0][0], k_[1][0], k_[2][0]] for k_ in keys], int); areas = np.zeros(len(keys)); dAs = np.zeros((len(keys), 3, 2)); polys = []
    for q, key in enumerate(keys):
        P, A, dA = zone_geom(key, X); areas[q] = A; dAs[q] = dA; polys.append(P)
    return ids, areas, dAs, polys
keys = list(triangulate(X).keys()); ids, areas, dAs, polys = all_geom(keys, X)
mT = 1.0*areas; KT = P0/(mT/areas)**G_
mA = np.zeros(N); np.add.at(mA, ids.ravel(), np.repeat(mT/3, 3))
zstate = {k_: (mT[q], KT[q]) for q, k_ in enumerate(keys)}
def forces(ids, areas, dAs, mTz, KTz):
    P = KTz*(mTz/areas)**G_ - P0; F = np.zeros((N, 2))
    for v in range(3): np.add.at(F, ids[:, v], P[:, None]*dAs[:, v])
    return F
def energy(Vv, areas, mTz, KTz): return 0.5*(mA*(Vv*Vv).sum(1)).sum() + (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1) + P0*areas).sum()
def wrap(X):
    Xw = X.copy(); Xw[:, 0] = np.mod(Xw[:, 0], Lx); return Xw
def rekey(keys, zs, X, Xw):
    shift = np.round((X[:, 0] - Xw[:, 0])/Lx).astype(int); zs2 = {}; keys2 = []
    for k_ in keys:
        kk = tuple((i, s_ + shift[i]) for i, s_ in k_); P = np.array([Xw[i] + [s_*Lx, 0] for i, s_ in kk]); off = int(np.floor(P.mean(0)[0]/Lx))
        kk = canon(tuple((i, s_ - off) for i, s_ in kk)); zs2[kk] = zs[k_]; keys2.append(kk)
    return keys2, zs2
def step(X, Vv, keys, zstate, dt):
    mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
    def acc(Xs):
        ids_, areas_, dAs_, _ = all_geom(keys, Xs)
        if (areas_ <= 0).any(): raise RuntimeError('inverted zone within a step (min signed area/d2 %.2e)' % (areas_.min()/d**2))
        a = forces(ids_, areas_, dAs_, mTz, KTz)/mA[:, None]; return a
    k1x = Vv; k1v = acc(X); k2x = Vv + 0.5*dt*k1v; k2v = acc(X + 0.5*dt*k1x)
    k3x = Vv + 0.5*dt*k2v; k3v = acc(X + 0.5*dt*k2x); k4x = Vv + dt*k3v; k4v = acc(X + dt*k3x)
    X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
    ids, areas, dAs, polys = all_geom(keys, X)
    if (areas <= 0).any(): raise RuntimeError('inverted zone at step end (min signed area/d2 %.2e)' % (areas.min()/d**2))
    Eb = energy(Vv, areas, mTz, KTz)
    Xw = wrap(X); keys, zs = rekey(keys, zstate, X, Xw); X = Xw
    if FIXCONN: return X, Vv, keys, zs, Eb, 0, 0.0
    ids, areas, dAs, polys = all_geom(keys, X)
    newz = triangulate(X); newkeys = list(newz.keys()); fresh = [k_ for k_ in newkeys if k_ not in zs]; dS = 0.0
    if fresh:
        gone = [k_ for k_ in keys if k_ not in newz]; gidx = {k_: q for q, k_ in enumerate(keys)}
        if REMAP2:
            cz_all = np.array([p.mean(0) for p in polys]); rho_all = mTz/areas; eps_all = (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1))/areas
            edge = {}
            for q, k_ in enumerate(keys):
                for a_ in range(3):
                    u_, w_ = sorted((k_[a_], k_[(a_ + 1) % 3])); e = (u_[0], w_[0], w_[1] - u_[1]); edge.setdefault(e, []).append(q)
            nbr = [[] for _ in keys]
            for e, qs in edge.items():
                if len(qs) == 2: nbr[qs[0]].append(qs[1]); nbr[qs[1]].append(qs[0])
            grad = {}
            for k_ in gone:
                i = gidx[k_]; nb = nbr[i]
                if len(nb) < 2: grad[k_] = (np.zeros(2), np.zeros(2)); continue
                D = np.array([cz_all[j] - cz_all[i] for j in nb]); D[:, 0] -= Lx*np.round(D[:, 0]/Lx)
                out = []
                for qv in (rho_all, eps_all):
                    g = np.linalg.lstsq(D, qv[nb] - qv[i], rcond=None)[0]; qmin = min(qv[i], qv[nb].min()); qmax = max(qv[i], qv[nb].max()); phi = 1.0
                    for v in polys[i]:
                        dq = g @ (v - cz_all[i])
                        if dq > 1e-300: phi = min(phi, (qmax - qv[i])/dq)
                        elif dq < -1e-300: phi = min(phi, (qmin - qv[i])/dq)
                    out.append(max(phi, 0.0)*g)
                grad[k_] = tuple(out)
        oldpolys = [Polygon(polys[gidx[k_]]) for k_ in gone]; tree = STRtree(oldpolys)
        rho_old = np.array([mTz[gidx[k_]]/areas[gidx[k_]] for k_ in gone]); e_old = np.array([(KTz[gidx[k_]]*mTz[gidx[k_]]**G_*areas[gidx[k_]]**(1 - G_)/(G_ - 1))/areas[gidx[k_]] for k_ in gone])
        A_old = np.array([areas[gidx[k_]] for k_ in gone]); cov_old = np.zeros(len(gone))
        S_old = sum(zs[k_][0]*np.log(zs[k_][1]) for k_ in gone); m_old = sum(zs[k_][0] for k_ in gone); S_new = 0.0; m_sum = 0.0; mH = {}; uH = {}
        for k_ in fresh:
            Pn, An, _ = zone_geom(k_, X); m_new = 0.0; u_new = 0.0; cov = 0.0; mh = 0.0; uh = 0.0
            for sx in (-1, 0, 1):
                pq = Polygon(Pn + [sx*Lx, 0])
                for qi in tree.query(pq):
                    op = pq.intersection(oldpolys[qi]); ov = op.area
                    if ov > 0:
                        m_new += ov*rho_old[qi]; u_new += ov*e_old[qi]; cov += ov; cov_old[qi] += ov
                        if REMAP2:
                            gr, ge = grad[gone[qi]]; dc = np.array(op.centroid.coords[0]) - polys[gidx[gone[qi]]].mean(0)
                            mh += ov*(rho_old[qi] + gr @ dc); uh += ov*(e_old[qi] + ge @ dc)
            unc = An - cov
            if unc > 1e-9*An and not OUTR[[k_[0][0], k_[1][0], k_[2][0]]].any(): raise RuntimeError('remap: interior new zone not covered')
            if m_new <= 0: raise RuntimeError('remap: new zone without overlap')
            K_new = (G_ - 1)*u_new/An/(m_new/An)**G_; zs[k_] = (m_new, K_new); S_new += m_new*np.log(K_new); m_sum += m_new; mH[k_] = mh; uH[k_] = uh
        lost = A_old - cov_old
        if (lost > 1e-9*A_old).any() and not OUTR[np.array([[k_[0][0], k_[1][0], k_[2][0]] for k_ in gone])[lost > 1e-9*A_old]].any(axis=1).all(): raise RuntimeError('remap: old interior zone not covered')
        if abs(m_sum - m_old) > 1e-9*max(m_old, 1e-30): raise RuntimeError('remap: mass mismatch %.3e vs %.3e' % (m_sum, m_old))
        dS = S_new - S_old
        if REMAP2:
            AA = {k_: zone_geom(k_, X)[1] for k_ in fresh}
            def state(theta):
                st = {}; S = 0.0
                for k_ in fresh:
                    mL, KL = zs[k_]; An = AA[k_]; uL = KL*mL**G_*An**(1 - G_)/(G_ - 1); m = mL + theta*(mH[k_] - mL); u = uL + theta*(uH[k_] - uL)
                    if m <= 0 or u <= 0: return None, -np.inf
                    K = (G_ - 1)*u/An/(m/An)**G_; st[k_] = (m, K); S += m*np.log(K)
                return st, S
            theta = 1.0
            while theta > 1e-3:
                st, S = state(theta)
                if st is not None and S >= S_old - 1e-12*abs(S_old): zs.update(st); dS = S - S_old; break
                theta *= 0.5
        for k_ in gone: del zs[k_]
    return X, Vv, newkeys, zs, Eb, len(fresh), dS
# reference (unseeded) run positions for the linear-response difference
REF = os.environ.get('REF', ''); ref = np.load(REF) if REF else None
mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys]); E0 = energy(Vv, areas, mTz, KTz)
S_tot = 0.0; Eremap = 0.0; nremap = 0; nrej = 0; s = 0; dt0 = dt; t = 0.0; t0 = time.time(); every = max(1, int(round(1.0/dt))); hist = []
def mode_amp(X):
    x = X[sel, 0]; dx = x - np.mod(X0[sel, 0] + S_*(X0[sel, 1] - Ly/2)*t, Lx); dx -= Lx*np.round(dx/Lx)
    if ref is not None and len(hist) < len(ref): dx = dx - ref[len(hist)]
    return np.abs(np.sum(dx*np.exp(-2j*np.pi*smode*X0[sel, 0]/Lx)))*2/len(sel)/d, dx
while t < T - 1e-12:
    s += 1; dtn = min(dt0, T - t); dts = dtn; nsub = 1; done = 0; ok = True; saved = (X.copy(), Vv.copy(), list(keys), dict(zstate), S_tot, Eremap, nremap)
    while done < nsub:
        try:
            X, Vv, keys, zstate, Eb, nf, dS = step(X, Vv, keys, zstate, dts); done += 1; nremap += nf; S_tot += dS
            ids, areas, dAs, polys = all_geom(keys, X); mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
            Ea = energy(Vv, areas, mTz, KTz); Eremap += Ea - Eb
        except RuntimeError as ex:
            nrej += 1
            if nsub >= 16: print('ABORT t %.3f: %s' % (t, ex), flush=True); ok = False; break
            X, Vv, keys, zstate, S_tot, Eremap, nremap = saved; nsub *= 2; dts = dtn/nsub; done = 0
    if not ok: break
    t += dtn
    if s % every == 0:
        amp, dx = mode_amp(X); hist.append(dx)
        rho = mTz/areas; Pz = KTz*rho**G_
        print('t %6.2f | seed m%d amp %.3e | P/P0 rms %.2e | dE/E0 %+.2e (remap %+.2e) | S_remap %+.3e | remapped %d rej %d | min area/d2 %.3f  [%.0fs]' % (
            t, smode, amp, np.sqrt(np.mean((Pz/P0 - 1)**2)), (Ea - E0)/abs(E0), Eremap/abs(E0), S_tot, nremap, nrej, areas.min()/d**2, time.time() - t0), flush=True)
np.save(tag + '_dx.npy', np.array(hist))
