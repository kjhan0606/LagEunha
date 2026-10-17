"""Static dispersion test of the alternating (odd-even) longitudinal mode on ring 1 for three force laws
(approved 2026-10-08). All geometry from lag.py (power diagram + disc closure, particle = generator).
 (1) Voronoi/EOS   : P_b = K(m/V_vor)^g,     dF = sum_b dP_b dV_b/dx_a            -> w^2_vor
 (2) Laguerre-DAE  : P_b = K(m/V^T)^g, V^T follows face kinematics dV^T_b = sum_faces A (dx_a+dx_j)/2 . e,
                     weights keep V(x,W)=V^T (zero at linear order about the IC)  -> w^2_dae
 (3) OT/Brenier    : weights re-solved so that V(x,W) = V_IC after the displacement; force (m/tau^2)(b_a - x_a)
                     with b_a the Laguerre-cell centroid                           -> k_c = restored fraction (x tau^-2)
Output in units cs^2/d^2 for (1),(2); dimensionless for (3)."""
import sys, numpy as np
from lag import ring_ic, candidates, cells, solve_weights
from shapely.geometry import Polygon
import shapely
G_ = 5/3; P0 = 1e-4
def centroids(X, w, Rc, C, nseg=256):
    N = len(X); B = np.zeros_like(X)
    circ0 = np.array([[np.cos(t), np.sin(t)] for t in np.linspace(0, 2*np.pi, nseg, endpoint=False)])*Rc
    for a in range(N):
        xa = X[a]; Bx = 8*Rc + 1.0
        cell = Polygon([(xa[0]-Bx, xa[1]-Bx), (xa[0]+Bx, xa[1]-Bx), (xa[0]+Bx, xa[1]+Bx), (xa[0]-Bx, xa[1]+Bx)])
        for j in C[a]:
            dv = X[j] - xa; dn = np.linalg.norm(dv); e = dv/dn; s = 0.5*dn + (w[a] - w[j])/(2*dn)
            pt = xa + s*e; t = np.array([-e[1], e[0]])*50; nrm = e*50
            cell = cell.intersection(Polygon([pt + t, pt - t, pt - t - nrm, pt + t - nrm]))
        cell = cell.intersection(Polygon(xa + circ0)); c = cell.centroid; B[a] = (c.x, c.y)
    return B
n = int(sys.argv[1]); modes = [int(v) for v in sys.argv[2].split(',')]
X, lab, d = ring_ic(n); Rc = d; N = len(X); C = candidates(X); w0 = np.zeros(N)
V0, G0, faces0 = cells(X, w0, Rc, C); m = V0.copy(); K = P0; cs2 = G_*P0
sel = np.where(lab == 1)[0]; Nr = len(sel); th = np.arctan2(X[sel, 1], X[sel, 0]); R = np.hypot(*X[sel].T)
et = np.stack([-np.sin(th), np.cos(th)], 1); eps = 1e-4
print('n %d N %d ring1 N %d  cs/d %.4f  (units: w^2/(cs^2/d^2); 1D chain sin^2(kd); continuum (kd)^2)' % (n, N, Nr, np.sqrt(cs2)/d), flush=True)
# face list per cell with (j, A, e) for the kinematic operator D
FL = [[] for _ in range(N)]
for a, j, A, dn in faces0:
    e = (X[j] - X[a])/dn; FL[a].append((j, A, e)); FL[j].append((a, A, -e))
def force_from_dP(dP):
    dF = np.zeros_like(X)
    for b in range(N):
        for k, g in G0[b].items(): dF[k] += dP[b]*g
    return dF
for m_ in modes:
    c_ = np.cos(m_*(th - th[0])); dX = np.zeros_like(X); dX[sel] = eps*d*c_[:, None]*et; kd = 2*np.pi*m_/Nr   # phase: Nyquist seed = (-1)^i
    proj = lambda F: 2*np.sum((F[sel]*et).sum(1)*c_)/Nr
    # (1) Voronoi: dV from exact geometry
    Xp = X + dX; Vp, Gp, _ = cells(Xp, w0, Rc, candidates(Xp)); dP = -G_*P0*(Vp - V0)/V0
    Fp = np.zeros_like(X)
    for b in range(N):
        for k, g in Gp[b].items(): Fp[k] += (P0 + dP[b])*g
    F0 = force_from_dP(np.full(N, P0)); w2_vor = -proj((Fp - F0)/m[:, None])/(eps*d)/(cs2/d/d)
    # (2) DAE: dV^T from midpoint face kinematics, force with IC geometry
    dVT = np.zeros(N)
    for a in range(N):
        for j, A, e in FL[a]: dVT[a] += A*0.5*((dX[a] + dX[j]) @ e)
    dP2 = -G_*P0*dVT/V0; w2_dae = -proj(force_from_dP(dP2)/m[:, None])/(eps*d)/(cs2/d/d)
    # (3) OT/Brenier: weights restore V0 after displacement; centroid force
    Cp = candidates(Xp); Vv, _, _ = cells(Xp, w0, Rc, Cp); wp, Vl, _, _ = solve_weights(Xp, V0*Vv.sum()/V0.sum(), Rc, Cp, verbose=False)
    Bc = centroids(Xp, wp, Rc, Cp); B0 = centroids(X, w0, Rc, C)
    kc = -proj((Bc - Xp) - (B0 - X))/(eps*d)
    print('m %3d kd/pi %.2f | Voronoi %+.4f | Laguerre-DAE %+.4f | chain %.4f | continuum %.3f | OT-centroid restored fraction %+.4f (max|w|/d^2 %.1e)' % (
        m_, kd/np.pi, w2_vor, w2_dae, np.sin(kd)**2, kd*kd, kc, np.abs(wp).max()/d/d), flush=True)
