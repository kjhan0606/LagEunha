"""stagP2.py: staggered (zone-centred) Kepler testbed, revision after Astra 27 (2026-10-10).
Changes vs stagP.py: (1) signed-area inversion check against the orientation fixed when the zone was created;
(2) no frozen wall and no rotation taper: the gas (disc + co-rotating background, Kepler gravity everywhere) fills a
disc of radius RMAX whose outer ring is a tracked material boundary in contact with a PRESSURE RESERVOIR P_ext = P0
(no mass exchange): F_a = sum_T (P_T - P0) dA_T/dx_a, conserved energy E = KE + sum_T(u_T + P0 A_T) + Phi;
triangles whose three vertices are all outer-ring particles are exterior (not zones); (3) no fallback of any kind:
a fresh zone without overlap rejects the step; both overlap closures (old covered, new covered) are checked, the
uncovered new area is allowed only for zones touching the boundary ring and is booked as boundary area change;
(4) rejected sub-steps restore every ledger; (5) remap diagnostic = change of nodal force at the remap (old zones vs
new zones at the same positions), rms per region and projection on the seeded mode.
args: n dt T P0 seed_R seed_m seed_amp tag BG   env: RMAX (3.6), GMF (1), RIGID (0)"""
import sys, os, time, numpy as np
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely.strtree import STRtree
G_ = 5/3
n = int(sys.argv[1]); dt = float(sys.argv[2]); T = float(sys.argv[3]); P0 = float(sys.argv[4]); sR = float(sys.argv[5])
smode = int(sys.argv[6]); samp = float(sys.argv[7]); tag = sys.argv[8]; BG = float(sys.argv[9])
DK = float(os.environ.get('DK', '0')); REMAP2 = int(os.environ.get('REMAP2', '0')); RMAX = float(os.environ.get('RMAX', '3.6')); L = 6.0; C = 3.0; d = L/n; RC = 0.5
def Mspl(u):
    u = np.asarray(u, float)
    return np.where(u >= 1, 1.0, np.where(u <= 0.5, u**3*(32/3 - 192/5*u*u + 32*u**3), -1/15 + u**3*(64/3 - 48*u + 192/5*u*u - 32/3*u**3)))
GMF = float(os.environ.get('GMF', '1')); _M0 = Mspl; Mspl = lambda u: GMF*_M0(u)
if os.environ.get('RIGID', '0') == '1':   # rigid rotation control: harmonic potential, Omega matched to Kepler at R=2.5
    _W0 = np.sqrt(_M0(2.5/RC)/2.5**3); Mspl = lambda u: _W0**2*(np.asarray(u, float)*RC)**3
def gacc(X):
    dx = X - C; r = np.hypot(*dx.T); g = -Mspl(r/RC)/r**3; return g[:, None]*dx
_rg = np.linspace(1e-4, 8.0, 32001); _gg = Mspl(_rg/RC)/_rg**2
_ph = np.concatenate([[0.0], np.cumsum(0.5*(_gg[1:] + _gg[:-1])*np.diff(_rg))]); _ph -= _ph[-1]
def Phi(X): dx = X - C; r = np.hypot(*dx.T); return np.interp(r, _rg, _ph)
# ---- IC: uniform rings (golden-ratio phases) out to RMAX; the last ring is the material boundary ----
xs = []; ys = []; k = 1; kmax = int(round(RMAX/d))
for k in range(1, kmax + 1):
    r = k*d; m_ = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(m_) + (0.6180339887*k) % 1.0)*2*np.pi/m_
    xs.append(C + r*np.cos(th)); ys.append(C + r*np.sin(th))
X = np.stack([np.concatenate(xs), np.concatenate(ys)], 1); N = len(X)
R0 = np.hypot(X[:, 0] - C, X[:, 1] - C); th0 = np.arctan2(X[:, 1] - C, X[:, 0] - C); OUT = R0 > (kmax - 0.5)*d
def Sig(p): r = np.hypot(p[:, 0] - C, p[:, 1] - C); return np.where((r >= 0.5) & (r <= 2.0), 1.0, BG)
sel = np.where(np.abs(R0 - sR) < 0.5*d)[0]; dth = samp*d/R0[sel]*np.cos(smode*(th0[sel] - th0[sel][0]))
X[sel] = np.stack([C + R0[sel]*np.cos(th0[sel] + dth), C + R0[sel]*np.sin(th0[sel] + dth)], 1)
R = np.hypot(X[:, 0] - C, X[:, 1] - C); vc = np.sqrt(Mspl(R/RC)/R); Vv = np.stack([-vc*(X[:, 1] - C)/R, vc*(X[:, 0] - C)/R], 1)
# ---- zones: Delaunay triangles, oriented CCW at creation; exterior slivers (3 boundary vertices) excluded ----
def triangulate(X):
    tri = Delaunay(X); zones = {}
    for simp in tri.simplices:
        i3 = [int(i) for i in simp]
        if OUT[i3].all(): continue
        a, b, c = X[i3]; A2 = (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])
        if A2 < 0: i3 = [i3[0], i3[2], i3[1]]
        zones[tuple(i3)] = None   # key = CCW-ordered vertex ids (orientation fixed at creation)
    return zones
def zone_geom(key, X):
    a, b, c = X[list(key)]; A = 0.5*((b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1]))   # signed w.r.t. creation orientation
    def J(v): return np.array([-v[1], v[0]])
    return np.array([a, b, c]), A, [0.5*J(c - b), 0.5*J(a - c), 0.5*J(b - a)]
from p2_geometry_backend import all_geom

def canon(key):   # same triangle, any cyclic order -> canonical (min id first); mirror order = different (inverted) zone
    i = int(np.argmin(key)); return (key[i], key[(i + 1) % 3], key[(i + 2) % 3])
zones = {canon(k_): None for k_ in triangulate(X)}; keys = list(zones.keys()); ids, areas, dAs, polys = all_geom(keys, X)
cen = np.array([p.mean(0) for p in polys]); mT = Sig(cen)*areas; KT = P0/(mT/areas)**G_
mA = np.zeros(N); np.add.at(mA, ids.ravel(), np.repeat(mT/3, 3))   # particle inertia, fixed
zstate = {key: (mT[q], KT[q]) for q, key in enumerate(keys)}
K0 = P0/1.0**G_
def shapeE(ids, areas, X):
    # E_D = kappa sum_T A_T (1 - q_T),  q_T = 4 sqrt3 A_T / L2_T,  L2_T = sum of squared edge lengths; returns (E_D, dE_D/dx)
    if DK == 0: return 0.0, np.zeros((N, 2))
    kap = DK*P0; P3 = X[ids]; e = np.stack([P3[:, 2] - P3[:, 1], P3[:, 0] - P3[:, 2], P3[:, 1] - P3[:, 0]], 1)   # edge opposite vertex v
    L2 = (e*e).sum((1, 2)); q = 4*np.sqrt(3)*areas/L2; E = kap*(areas*(1 - q)).sum()
    # dE/dx_a = kap[(1-q) dA/dx_a - A dq/dx_a],  dq/dx_a = 4sqrt3 (dA/dx_a / L2 - A dL2/dx_a / L2^2),  dL2/dx_a = 2*(e_b... ) : vertex a appears in two edges
    def J(v): return np.stack([-v[..., 1], v[..., 0]], -1)
    dA = 0.5*J(e)                       # dA/dx_v for vertex v (CCW)
    dL2 = np.zeros_like(dA)
    dL2[:, 0] = 2*(e[:, 1]*(-1) + e[:, 2]*(1))*0 + 2*((P3[:, 0] - P3[:, 2]) + (P3[:, 0] - P3[:, 1]))
    dL2[:, 1] = 2*((P3[:, 1] - P3[:, 0]) + (P3[:, 1] - P3[:, 2]))
    dL2[:, 2] = 2*((P3[:, 2] - P3[:, 1]) + (P3[:, 2] - P3[:, 0]))
    dq = 4*np.sqrt(3)*(dA/L2[:, None, None] - (areas/L2**2)[:, None, None]*dL2)
    dE = kap*((1 - q)[:, None, None]*dA - areas[:, None, None]*dq); G = np.zeros((N, 2))
    for v in range(3): np.add.at(G, ids[:, v], dE[:, v])
    return E, G
def forces(ids, areas, dAs, mTz, KTz, X=None):
    P = KTz*(mTz/areas)**G_ - P0; F = np.zeros((N, 2))   # pressure relative to the reservoir
    for v in range(3): np.add.at(F, ids[:, v], P[:, None]*dAs[:, v])
    if DK and X is not None: F -= shapeE(ids, areas, X)[1]
    return F
def energy(X, Vv, areas, mTz, KTz, ids=None):
    ED = shapeE(ids, areas, X)[0] if (DK and ids is not None) else 0.0
    return 0.5*(mA*(Vv*Vv).sum(1)).sum() + (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1) + P0*areas).sum() + (mA*Phi(X)).sum() + ED
mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
E0 = energy(X, Vv, areas, mTz, KTz, ids); A0 = areas.sum()
cs_bg = np.sqrt(G_*P0/BG); cs_d = np.sqrt(G_*P0)
print('stagP2 DK %.2f REMAP2 %d n %d N %d (boundary ring %d) zones %d P0 %.2e BG %.0e c_s bg %.4f disc %.4f d %.4f (c_s/d bg %.3f disc %.3f) RMAX %.2f seed R %.2f m %d amp %.3f d  E0 %.6e' % (
    DK, REMAP2, n, N, OUT.sum(), len(keys), P0, BG, cs_bg, cs_d, d, cs_bg/d, cs_d/d, RMAX, sR, smode, samp, E0), flush=True)
# seeded-mode pattern for the force-projection diagnostic
tsel = np.stack([-np.sin(th0[sel]), np.cos(th0[sel])], 1)*np.cos(smode*(th0[sel] - th0[sel][0]))[:, None]
def step(X, Vv, keys, zstate, dt):
    mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
    def acc(Xs):
        ids_, areas_, dAs_, _ = all_geom(keys, Xs)
        if (areas_ <= 0).any(): raise RuntimeError('inverted zone within a step (min signed area/d2 %.2e)' % (areas_.min()/d**2))
        return forces(ids_, areas_, dAs_, mTz, KTz, Xs)/mA[:, None] + gacc(Xs)
    k1x = Vv; k1v = acc(X); k2x = Vv + 0.5*dt*k1v; k2v = acc(X + 0.5*dt*k1x)
    k3x = Vv + 0.5*dt*k2v; k3v = acc(X + 0.5*dt*k2x); k4x = Vv + dt*k3v; k4v = acc(X + dt*k3x)
    X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
    ids, areas, dAs, polys = all_geom(keys, X)
    if (areas <= 0).any(): raise RuntimeError('inverted zone at step end (min signed area/d2 %.2e)' % (areas.min()/d**2))
    Eb = energy(X, Vv, areas, mTz, KTz, ids); Fb = forces(ids, areas, dAs, mTz, KTz, X)
    newz = {canon(k_): None for k_ in triangulate(X)}; newkeys = list(newz.keys()); zs = dict(zstate); fresh = [k_ for k_ in newkeys if k_ not in zs]
    dS = 0.0; dAb = 0.0; dF = None
    if fresh:
        gone = [k_ for k_ in keys if k_ not in newz]; gidx = {k_: q for q, k_ in enumerate(keys)}
        if REMAP2:
            # gradients of q = (rho, eps) for every old zone from edge neighbours (least squares), limited (BJ)
            cz_all = np.array([p.mean(0) for p in polys]); rho_all = mTz/areas; eps_all = (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1))/areas
            edge = {}
            for q, k_ in enumerate(keys):
                for a_ in range(3):
                    e = tuple(sorted((k_[a_], k_[(a_ + 1) % 3]))); edge.setdefault(e, []).append(q)
            nbr = [[] for _ in keys]
            for e, qs in edge.items():
                if len(qs) == 2: nbr[qs[0]].append(qs[1]); nbr[qs[1]].append(qs[0])
            grad = {}
            for k_ in gone:
                i = gidx[k_]; nb = nbr[i]
                if len(nb) < 2: grad[k_] = (np.zeros(2), np.zeros(2)); continue
                D = cz_all[nb] - cz_all[i]
                gr = np.linalg.lstsq(D, rho_all[nb] - rho_all[i], rcond=None)[0]; ge = np.linalg.lstsq(D, eps_all[nb] - eps_all[i], rcond=None)[0]
                out = []
                for g, qv in ((gr, rho_all), (ge, eps_all)):
                    qmin = min(qv[i], qv[nb].min()); qmax = max(qv[i], qv[nb].max()); phi = 1.0
                    for v in polys[i]:
                        dq = g @ (v - cz_all[i])
                        if dq > 1e-300: phi = min(phi, (qmax - qv[i])/dq)
                        elif dq < -1e-300: phi = min(phi, (qmin - qv[i])/dq)
                    out.append(max(phi, 0.0)*g)
                grad[k_] = (out[0], out[1])
            mH = {}; uH = {}
        oldpolys = [Polygon(polys[gidx[k_]]) for k_ in gone]; tree = STRtree(oldpolys)
        rho_old = np.array([mTz[gidx[k_]]/areas[gidx[k_]] for k_ in gone]); e_old = np.array([(KTz[gidx[k_]]*mTz[gidx[k_]]**G_*areas[gidx[k_]]**(1 - G_)/(G_ - 1))/areas[gidx[k_]] for k_ in gone])
        A_old = np.array([areas[gidx[k_]] for k_ in gone]); cov_old = np.zeros(len(gone))
        S_old = sum(zs[k_][0]*np.log(zs[k_][1]) for k_ in gone); m_old = sum(zs[k_][0] for k_ in gone)
        S_new = 0.0; m_sum = 0.0
        for k_ in fresh:
            Pn, An, _ = zone_geom(k_, X); pn = Polygon(Pn); m_new = 0.0; u_new = 0.0; cov = 0.0
            mh = 0.0; uh = 0.0
            for qi in tree.query(pn):
                op = pn.intersection(oldpolys[qi]); ov = op.area
                if ov > 0:
                    m_new += ov*rho_old[qi]; u_new += ov*e_old[qi]; cov += ov; cov_old[qi] += ov
                    if REMAP2:
                        gr, ge = grad[gone[qi]]; dc = np.array(op.centroid.coords[0]) - polys[gidx[gone[qi]]].mean(0)
                        mh += ov*(rho_old[qi] + gr @ dc); uh += ov*(e_old[qi] + ge @ dc)
            if REMAP2: mH[k_] = mh; uH[k_] = uh
            unc = An - cov
            if unc < -1e-9*An: raise RuntimeError('remap: new zone over-covered (duplicate old zones) %.2e' % (-unc/An))
            if unc > 1e-9*An:
                if OUT[list(k_)].any(): dAb += unc   # boundary polygon gained area from the reservoir (no mass): booked
                else: raise RuntimeError('remap: interior new zone not covered (uncovered %.2e of A)' % (unc/An))
            if m_new <= 0: raise RuntimeError('remap: new zone without overlap')
            K_new = (G_ - 1)*u_new/An/(m_new/An)**G_; zs[k_] = (m_new, K_new); S_new += m_new*np.log(K_new); m_sum += m_new
        lost = A_old - cov_old
        if (lost < -1e-9*A_old).any(): raise RuntimeError('remap: old zone over-covered (overlapping new zones) %.2e' % (-(lost/A_old).min()))
        if (lost > 1e-9*A_old).any():
            if not OUT[np.array([list(k_) for k_ in gone])[lost > 1e-9*A_old]].any(axis=1).all(): raise RuntimeError('remap: old interior zone not fully covered (max %.2e of A)' % (lost/A_old).max())
            dAb -= lost.sum()
        if abs(m_sum - m_old) > 1e-9*max(m_old, 1e-30): raise RuntimeError('remap: mass mismatch %.3e vs %.3e' % (m_sum, m_old))
        dS = S_new - S_old
        if REMAP2 and fresh:
            global n_theta, n_blend
            AA = {k_: zone_geom(k_, X)[1] for k_ in fresh}
            def state(theta):
                st = {}; S = 0.0
                for k_ in fresh:
                    mL, KL = zs[k_]; An = AA[k_]; uL = KL*mL**G_*An**(1 - G_)/(G_ - 1)
                    m = mL + theta*(mH[k_] - mL); u = uL + theta*(uH[k_] - uL)
                    if m <= 0 or u <= 0: return None, -np.inf
                    K = (G_ - 1)*u/An/(m/An)**G_; st[k_] = (m, K); S += m*np.log(K)
                return st, S
            theta = 1.0; chosen = None
            while theta > 1e-3:
                st, S = state(theta)
                if st is not None and S >= S_old - 1e-12*abs(S_old): chosen = (st, S, theta); break
                theta *= 0.5
            n_blend += 1
            if chosen is not None:
                st, S, theta = chosen
                mB = sum(v[0] for v in st.values()); uB = sum(v[1]*v[0]**G_*AA[k_]**(1 - G_)/(G_ - 1) for k_, v in st.items())
                uL_tot = sum(zs[k_][1]*zs[k_][0]**G_*AA[k_]**(1 - G_)/(G_ - 1) for k_ in fresh)
                if abs(mB - m_sum) > 1e-9*m_sum or abs(uB - uL_tot) > 1e-9*uL_tot: raise RuntimeError('remap2: blended state not conservative (dm %.2e du %.2e)' % ((mB - m_sum)/m_sum, (uB - uL_tot)/uL_tot))
                zs.update(st); dS = S - S_old
                if theta < 1: n_theta += 1
            else: n_theta += 1   # fall back to first order (theta = 0)
        for k_ in gone: del zs[k_]
        ids2, areas2, dAs2, _ = all_geom(newkeys, X); m2 = np.array([zs[k_][0] for k_ in newkeys]); K2 = np.array([zs[k_][1] for k_ in newkeys])
        dF = forces(ids2, areas2, dAs2, m2, K2, X) - Fb
    return X, Vv, newkeys, zs, Eb, len(fresh), dS, dAb, dF
n_theta = 0; n_blend = 0; S_tot = 0.0; Eremap = 0.0; nremap_tot = 0; Ab_tot = 0.0; nrej = 0; s = 0; dt0 = dt; t = 0.0; t0 = time.time(); every = int(round(1.0/dt))
DF = []   # (rms |dF| inner bg, rms disc, projection on seeded mode) per remap
while t < T - 1e-12:
    s += 1; dtn = min(dt0, T - t); dts = dtn; nsub = 1; done = 0; ok = True
    saved = (X.copy(), Vv.copy(), list(keys), dict(zstate), S_tot, Eremap, nremap_tot, Ab_tot, len(DF), n_theta, n_blend)
    while done < nsub:
        try:
            X, Vv, keys, zstate, Eb, nf, dS, dAb, dF = step(X, Vv, keys, zstate, dts); done += 1
            nremap_tot += nf; S_tot += dS; Ab_tot += dAb
            if dS < -1e-12*max(abs(S_tot), 1e-30): print('   WARNING: remap dS < 0: %.3e' % dS, flush=True)
            ids, areas, dAs, polys = all_geom(keys, X)
            mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
            Ea = energy(X, Vv, areas, mTz, KTz, ids); Eremap += Ea - Eb
            if dF is not None:
                fn = np.hypot(*dF.T)/(P0*d); bgi = R0 < 0.6; dsi = (R0 > 0.6) & (R0 < 1.9)
                DF.append((np.sqrt(np.mean(fn[bgi]**2)), np.sqrt(np.mean(fn[dsi]**2)), (dF[sel]*tsel).sum()/(P0*d*len(sel))))
        except RuntimeError as ex:
            nrej += 1
            if nsub >= 16: print('ABORT at t %.3f: %s (dt/16 failed)' % (t + done*dts, ex), flush=True); ok = False; break
            X, Vv, keys, zstate, S_tot, Eremap, nremap_tot, Ab_tot, ndf, n_theta, n_blend = saved; del DF[ndf:]
            nsub *= 2; dts = dtn/nsub; done = 0
            print('   reject t %.3f (%s) -> substeps %d' % (t, str(ex)[:70], nsub), flush=True)
            continue
    if not ok: break
    t += dtn
    if s % every == 0:
        th = np.arctan2(X[:, 1] - C, X[:, 0] - C); out = []
        for Rr in (sR, 0.41, 0.84, 1.25):
            kk = np.where(np.abs(R0 - Rr) < 0.5*d)[0]
            if len(kk) < 10: continue
            o = kk[np.argsort(th0[kk])]; nn = len(o); dth_ = np.angle(np.exp(1j*(th[o] - th0[o]))); dth_ -= dth_.mean()
            Fm = np.abs(np.fft.rfft(dth_))*2/nn; Fm[:2] = 0; q = np.argmax(Fm); out.append('R%.2f m%d %.2e' % (Rr, q, R0[o].mean()*Fm[q]/d))
        rho = mTz/areas; Pz = KTz*rho**G_; cz = np.array([p.mean(0) for p in polys]); Rz = np.hypot(*(cz - C).T); disc = (Rz > 0.6) & (Rz < 1.9)
        qP = np.argmax(Pz/P0); qA = np.argmin(areas/d**2); Rb = np.hypot(*(X[OUT] - C).T)
        print('   diag: max P/P0 %.2e at R %.3f (area/d2 %.3f) | min area/d2 %.3e at R %.3f | slivers<0.05d2 %d | boundary R %.4f-%.4f A_dom/A0-1 %+.2e remap-dA(cum) %+.2e' % (
            Pz[qP]/P0, Rz[qP], areas[qP]/d**2, areas[qA]/d**2, Rz[qA], int((areas < 0.05*d*d).sum()), Rb.min(), Rb.max(), areas.sum()/A0 - 1, Ab_tot), flush=True)
        if DF:
            df = np.array(DF); DF.clear()
            print('   remap dF/(P0 d): inner bg rms %.2e | disc rms %.2e | seeded-mode projection mean %+.2e rms %.2e (n %d)' % (df[:, 0].mean(), df[:, 1].mean(), df[:, 2].mean(), np.sqrt(np.mean(df[:, 2]**2)), len(df)), flush=True)
        Kd = (KTz[disc]*mTz[disc]).sum()/mTz[disc].sum()/K0; Ud = (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1))[disc].sum()
        print('t %5.1f | %s | zones %d remapped(cum) %d | dE/E0 %+.2e (remap %+.2e) | P/P0 disc rms %.2e | <K>/K0 disc %.6f U_disc %.5e | S_remap %+.3e rej %d | theta<1 %d/%d | mass %.6f  [%.0fs]' % (
            t, ' | '.join(out), len(keys), nremap_tot, (Ea - E0)/abs(E0), Eremap/abs(E0), np.sqrt(np.mean((Pz[disc]/P0 - 1)**2)), Kd, Ud, S_tot, nrej, n_theta, n_blend, mTz.sum(), time.time() - t0), flush=True)
