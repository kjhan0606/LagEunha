# Checked research copy; original sha256 b29927bdf22993f7cd28f45702ac7c38959334bfbd802962f32e4c588b494982
# Original: /gpfs/kjhan/LagForce/PolyBnd/stagP.py; baseline IC unchanged.
"""stagP.py: STAGGERED (zone-centred) formulation test, option B (approved 2026-10-09: B -> D).
Zones = Delaunay triangles of the particles (periodic box via image points). Each zone carries mass m_T and entropy K_T;
zone volume = triangle area A_T (exact material volume while the connectivity is fixed); P_T = K_T (m_T/A_T)^g.
U = sum_T u_T(A_T); force on particle a: F_a = sum_T P_T dA_T/dx_a (exact -grad U). Particle inertia m_a fixed (1/3 of
adjacent zone masses at t=0). The triangulation is frozen within an RK4 step; after the step the Delaunay is rebuilt and
NEW zones get mass and internal energy by exact area-overlap remap from the old zones (conservative, first order).
Setup as daeP.py: disc rho=1 (0.5<=r<=2), background BG, uniform P0, spline-softened point mass rc=0.5, rotation and gravity
tapered to zero over R=2.5..2.9. Diagnostics: seeded ring modes, energy (KE+U+Phi) with remap jumps booked separately,
number of remapped zones per step, K drift.
args: n dt T P0 seed_R seed_m seed_amp tag BG"""
import sys, time, numpy as np
from stag_checks import signed_area, check_orientation, overlap_remap, atomic_substeps
from scipy.spatial import Delaunay
from shapely.geometry import Polygon
from shapely.strtree import STRtree
G_ = 5/3
n = int(sys.argv[1]); dt = float(sys.argv[2]); T = float(sys.argv[3]); P0 = float(sys.argv[4]); sR = float(sys.argv[5])
smode = int(sys.argv[6]); samp = float(sys.argv[7]); tag = sys.argv[8]; BG = float(sys.argv[9])
import os
L = float(os.environ.get('BOXL', '6')); C = L/2; d = L/n; RC = 0.5; RT1, RT2 = 2.5, float(os.environ.get('RT2', '2.9'))   # BOXL=8 RT2=3.6: wide taper (shear ~ Kepler's), far field frozen at RT2
def taper(r):
    r = np.asarray(r, float); x = np.clip((r - RT1)/(RT2 - RT1), 0, 1); return 0.5*(1 + np.cos(np.pi*x))
def Mspl(u):
    u = np.asarray(u, float)
    return np.where(u >= 1, 1.0, np.where(u <= 0.5, u**3*(32/3 - 192/5*u*u + 32*u**3), -1/15 + u**3*(64/3 - 48*u + 192/5*u*u - 32/3*u**3)))
import os
GMF = float(os.environ.get('GMF', '1'))   # gravity (and circular speed) scaling for shear-rate tests
_M0 = Mspl; Mspl = lambda u: GMF*_M0(u)
if os.environ.get('RIGID', '0') == '1':   # rigid rotation control (no neighbour sliding): harmonic potential with Omega matched at R=0.41
    _W0 = np.sqrt(_M0(2.5/RC)/2.5**3)   # Omega matched to Kepler at R=2.5 so the taper/wall shear equals the baseline; interior slides nowhere
    Mspl = lambda u: _W0**2*(np.asarray(u, float)*RC)**3
def gacc(X):
    dx = X - C; r = np.hypot(*dx.T); g = -Mspl(r/RC)/r**3*taper(r)**2; return g[:, None]*dx
_rg = np.linspace(1e-4, 5.0, 20001); _gg = Mspl(_rg/RC)/_rg**2*taper(_rg)**2
_ph = np.concatenate([[0.0], np.cumsum(0.5*(_gg[1:] + _gg[:-1])*np.diff(_rg))]); _ph -= _ph[-1]
def Phi(X): dx = X - C; r = np.hypot(*dx.T); return np.interp(r, _rg, _ph)
# ---- IC (as daeP) ----
xs = []; ys = []; k = 1
while k*d < 3*np.sqrt(2) + d:
    r = k*d; m_ = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(m_) + (0.6180339887*k) % 1.0)*2*np.pi/m_
    xs.append(C + r*np.cos(th)); ys.append(C + r*np.sin(th)); k += 1
X = np.stack([np.concatenate(xs), np.concatenate(ys)], 1); inb = (X[:, 0] > 0) & (X[:, 0] < L) & (X[:, 1] > 0) & (X[:, 1] < L); X = X[inb]
from scipy.spatial import cKDTree
# (non-periodic: the far field R>2.9 is static, so no periodic images are needed and no points are dropped)
N = len(X); R0 = np.hypot(X[:, 0] - C, X[:, 1] - C); th0 = np.arctan2(X[:, 1] - C, X[:, 0] - C)
def Sig(p): r = np.hypot(p[:, 0] - C, p[:, 1] - C); return np.where((r >= 0.5) & (r <= 2.0), 1.0, BG)
sel = np.where(np.abs(R0 - sR) < 0.5*d)[0]; dth = samp*d/R0[sel]*np.cos(smode*(th0[sel] - th0[sel][0]))
X[sel] = np.stack([C + R0[sel]*np.cos(th0[sel] + dth), C + R0[sel]*np.sin(th0[sel] + dth)], 1)
R = np.hypot(X[:, 0] - C, X[:, 1] - C); vc = np.sqrt(Mspl(R/RC)/R)*taper(R); Vv = np.stack([-vc*(X[:, 1] - C)/R, vc*(X[:, 0] - C)/R], 1)
# ---- periodic Delaunay zones ----
def triangulate(X):
    tri = Delaunay(X); zones = {}
    for simp in tri.simplices:
        key = tuple(sorted((int(i), (0, 0)) for i in simp)); zones[key] = None
    return zones
def zone_geom(key, X, orientation=None):
    """vertex positions (representative frame) and ids; area; dA/dx for the 3 vertices"""
    ids = [k_[0] for k_ in key]; P = np.array([X[i] for i, s_ in key])
    if orientation is not None: check_orientation(P, orientation[key])
    a, b, c = P; A = 0.5*((b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1]))
    if A < 0: P = P[::-1]; ids = ids[::-1]; a, b, c = P; A = -A
    def J(v): return np.array([-v[1], v[0]])
    dA = [0.5*J(c - b), 0.5*J(a - c), 0.5*J(b - a)]   # dA/dx_a etc. (CCW)
    return ids, P, A, dA
def all_geom(zones, X, orientation=None):
    keys = list(zones.keys()); ids = np.zeros((len(keys), 3), int); areas = np.zeros(len(keys)); dAs = np.zeros((len(keys), 3, 2)); polys = []
    for q, key in enumerate(keys):
        i3, P, A, dA = zone_geom(key, X, orientation); ids[q] = i3; areas[q] = A; dAs[q] = dA; polys.append(P)
    return keys, ids, areas, dAs, polys
# ---- initial zone state ----
zones = triangulate(X); keys, ids, areas, dAs, polys = all_geom(zones, X)
cen = np.array([p.mean(0) for p in polys]); mT = Sig(cen)*areas; KT = P0/(mT/areas)**G_; UT = P0*areas/(G_ - 1)
mA = np.zeros(N); np.add.at(mA, ids.ravel(), np.repeat(mT/3, 3))   # particle inertia, fixed
zstate = {key: (mT[q], KT[q]) for q, key in enumerate(keys)}
K0 = P0/1.0**G_; K0bg = P0/BG**G_   # reference entropies (disc, background)
def forces(keys, ids, areas, dAs, mTz, KTz):
    P = KTz*(mTz/areas)**G_; F = np.zeros((N, 2))
    for v in range(3): np.add.at(F, ids[:, v], P[:, None]*dAs[:, v])
    return F
def energy(X, Vv, keys, ids, areas, mTz, KTz):
    U = (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1)).sum(); ev = np.ones(N, bool)
    return 0.5*(mA*(Vv*Vv).sum(1)).sum() + U + (mA*Phi(X)).sum()
mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
E0 = energy(X, Vv, keys, ids, areas, mTz, KTz); Eremap = 0.0; nremap_tot = 0
print('stagP n %d N %d zones %d P0 %.2e BG %.0e seed R %.2f m %d amp %.3f d  E0 %.6e' % (n, N, len(keys), P0, BG, sR, smode, samp, E0), flush=True)
t = 0.0; t0 = time.time(); every = int(round(1.0/dt))
FRZ = R0 >= RT2   # far field frozen (at rest, uniform P, no gravity): fixed convex hull, as the C setup freezes R>=2.5
Vv[FRZ] = 0
def step(X, Vv, keys, zstate, dt):
    mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
    orientation = {key: signed_area([X[i] for i, _ in key]) for key in keys}
    def acc(Xs):
        _, ids_, areas_, dAs_, _ = all_geom(dict.fromkeys(keys), Xs, orientation)
        if (areas_ <= 0).any(): raise RuntimeError('inverted zone within a step (min area %.2e)' % areas_.min())
        a_ = forces(keys, ids_, areas_, dAs_, mTz, KTz)/mA[:, None] + gacc(Xs); a_[FRZ] = 0
        return a_
    k1x = Vv; k1v = acc(X); k2x = Vv + 0.5*dt*k1v; k2v = acc(X + 0.5*dt*k1x)
    k3x = Vv + 0.5*dt*k2v; k3v = acc(X + 0.5*dt*k2x); k4x = Vv + dt*k3v; k4v = acc(X + dt*k3x)
    X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v)
    # end of step: energy before remap, rebuild triangulation, remap new zones by exact overlap (no fallback)
    _, ids, areas, dAs, polys = all_geom(dict.fromkeys(keys), X, orientation)
    if (areas <= 0).any(): raise RuntimeError('inverted zone at step end (min area %.2e)' % areas.min())
    Eb = energy(X, Vv, keys, ids, areas, mTz, KTz)
    newz = triangulate(X); newkeys = list(newz.keys()); zs = dict(zstate); fresh = [k_ for k_ in newkeys if k_ not in zs]
    dS = 0.0; jumps = np.empty((0, 2)); closure = np.zeros(4)
    if fresh:
        gone = [k_ for k_ in keys if k_ not in newz]
        index = {k_: q for q, k_ in enumerate(keys)}
        oi = [index[k_] for k_ in gone]
        oldpolys = [polys[q] for q in oi]
        newpolys = [zone_geom(k_, X)[1] for k_ in fresh]
        oldmass = mTz[oi]
        oldU = KTz[oi]*oldmass**G_*areas[oi]**(1-G_)/(G_-1)
        mn, un, closure = overlap_remap(oldpolys, newpolys, oldmass, oldU)
        S_old = sum(zs[k_][0]*np.log(zs[k_][1]) for k_ in gone)
        S_new = 0.0; jumps = []
        for j, k_ in enumerate(fresh):
            _, pn, an, _ = zone_geom(k_, X)
            knew = (G_-1)*un[j]/an/(mn[j]/an)**G_
            zs[k_] = (mn[j], knew); S_new += mn[j]*np.log(knew)
            # P reconstructed from remapped U must equal the overlap-volume
            # weighted donor P; no density weighting in this consistency check.
            pref = (G_-1)*un[j]/an
            pnew = knew*(mn[j]/an)**G_
            jumps.append((pnew/pref-1, np.hypot(*(pn.mean(0)-C))))
        dS = S_new-S_old; jumps = np.asarray(jumps)
        for k_ in gone: del zs[k_]
    return X, Vv, newkeys, zs, Eb, len(fresh), dS, jumps, closure

JUMP = []
S_tot = 0.0   # exact remap entropy ledger: sum over flips of (sum m ln K)_new - (sum m ln K)_old
nrej = 0; s = 0; dt0 = dt
mass0 = mTz.sum()
Lz0 = np.sum(mA*((X[:,0]-C)*Vv[:,1]-(X[:,1]-C)*Vv[:,0]))
closure_max = np.zeros(4)
def advance(state, dts):
    x, v, kk, zz = state
    xn, vn, kn, zn, eb, nf, ds, jumps, closure = step(x, v, kk, zz, dts)
    _, ii, aa, _, _ = all_geom(dict.fromkeys(kn), xn)
    mm = np.array([zn[k][0] for k in kn]); entropy = np.array([zn[k][1] for k in kn])
    ea = energy(xn, vn, kn, ii, aa, mm, entropy)
    return (xn, vn, kn, zn), (nf, ds, ea-eb, jumps, closure, ea)
while t < T - 1e-12:
    s += 1; dtn = min(dt0, T-t)
    try:
        (X, Vv, keys, zstate), records, rejected = atomic_substeps(
            (X, Vv, keys, zstate), dtn, advance)
    except RuntimeError as exc:
        print('ABORT at t %.9f: %s' % (t, exc), flush=True)
        raise SystemExit(1)
    nrej += len(rejected)
    for nsub, reason in rejected:
        print('reject t %.9f substeps %d: %s' % (t, nsub, reason), flush=True)
    for nf, ds, de, jumps, closure, Ea in records:
        nremap_tot += nf; S_tot += ds; Eremap += de
        if len(jumps): JUMP.append(jumps)
        closure_max = np.maximum(closure_max, closure)
    _, ids, areas, dAs, polys = all_geom(dict.fromkeys(keys), X)
    mTz = np.array([zstate[k_][0] for k_ in keys]); KTz = np.array([zstate[k_][1] for k_ in keys])
    t += dtn
    if s % every == 0:
        projected = np.zeros(N); np.add.at(projected, ids.ravel(), np.repeat(mTz/3, 3))
        lz = np.sum(mA*((X[:,0]-C)*Vv[:,1]-(X[:,1]-C)*Vv[:,0]))
        print('audit: mass_rel %+.3e Lz_rel %+.3e nodal_mass_l1 %.3e closure %s' % (
            (mTz.sum()-mass0)/mass0, (lz-Lz0)/max(abs(Lz0),1e-300),
            np.abs(projected-mA).sum()/mA.sum(), closure_max), flush=True)
        th = np.arctan2(X[:, 1] - C, X[:, 0] - C); out = []
        for Rr in (sR, 0.41, 0.84, 1.25):
            kk = np.where(np.abs(R0 - Rr) < 0.5*d)[0]
            if len(kk) < 10: continue
            o = kk[np.argsort(th0[kk])]; nn = len(o); dth_ = np.angle(np.exp(1j*(th[o] - th0[o]))); dth_ -= dth_.mean()
            Fm = np.abs(np.fft.rfft(dth_))*2/nn; Fm[:2] = 0; q = np.argmax(Fm); out.append('R%.2f m%d %.2e' % (Rr, q, R0[o].mean()*Fm[q]/d))
        rho = mTz/areas; Pz = KTz*rho**G_; disc = (np.hypot(*(np.array([p.mean(0) for p in polys]) - C).T) > 0.6) & (np.hypot(*(np.array([p.mean(0) for p in polys]) - C).T) < 1.9)
        cz = np.array([p.mean(0) for p in polys]); Rz = np.hypot(*(cz - C).T); qP = np.argmax(Pz/P0); qA = np.argmin(areas/d**2)
        print('   diag: max P/P0 %.2e at R %.3f (area/d2 %.3f) | min area/d2 %.3e at R %.3f (P/P0 %.2e) | zones with area<0.05d2: %d (R range %.2f-%.2f)' % (Pz[qP]/P0, Rz[qP], areas[qP]/d**2, areas[qA]/d**2, Rz[qA], Pz[qA]/P0, int((areas < 0.05*d*d).sum()), Rz[areas < 0.05*d*d].min() if (areas < 0.05*d*d).any() else 0, Rz[areas < 0.05*d*d].max() if (areas < 0.05*d*d).any() else 0), flush=True)
        if JUMP:
            jj = np.concatenate([j for j in JUMP if len(j)]); JUMP.clear(); bgi = jj[:, 1] < 0.6; dsi = (jj[:, 1] > 0.6) & (jj[:, 1] < 1.9)
            print('   remap EOS residual rms: inner bg %.2e (n %d) | disc %.2e (n %d) | max %.2e' % (np.sqrt(np.mean(jj[bgi, 0]**2)) if bgi.any() else 0, bgi.sum(), np.sqrt(np.mean(jj[dsi, 0]**2)) if dsi.any() else 0, dsi.sum(), np.abs(jj[:, 0]).max()), flush=True)
        Kd = (KTz[disc]*mTz[disc]).sum()/mTz[disc].sum()/K0; Ud = (KTz*mTz**G_*areas**(1 - G_)/(G_ - 1))[disc].sum()
        print('t %5.1f | %s | zones %d remapped(cum) %d | dE/E0 %+.2e (remap %+.2e) | P/P0 disc rms %.2e | <K>/K0 disc %.6f U_disc %.5e | S_remap %+.3e rej %d | mass %.6f  [%.0fs]' % (
            t, ' | '.join(out), len(keys), nremap_tot, (Ea - E0)/abs(E0), Eremap/abs(E0), np.sqrt(np.mean((Pz[disc]/P0 - 1)**2)), Kd, Ud, S_tot, nrej, mTz.sum(), time.time() - t0), flush=True)

print('COMPLETE t %.9f target %.9f' % (t, T), flush=True)
