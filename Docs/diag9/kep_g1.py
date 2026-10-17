"""Prescribed Kepler rotation (Omega=(R^2+eps^2)^-3/4, eps=0.01) of the res256 ICs (lattice / concentric rings).
Cell volume vs material volume (= own cell volume at t=0, incompressible flow) for Voronoi and G1 (metric G=F^-T F^-1,
exact deformation gradient). Cells closed by the disc phi_i <= rho^2 (rho = 1 spacing, as EUNHA_VORO_RHO=1).
Volumes by rasterization inside a 24-degree wedge; reported for cells whose generator lies in the central 12 degrees."""
import sys, struct, numpy as np
from scipy.spatial import cKDTree
ic, K = sys.argv[1], int(sys.argv[2]); PX = 20
C, RIN, ROUT, EPS = 3.0, 0.5, 2.0, 0.01; n = 384; d = 6.0/n; rho = 1.0*d
if ic == 'lattice':
    g = (np.arange(n)+0.5)*d; X, Y = np.meshgrid(g, g, indexing='ij'); x0 = np.stack([X.ravel(), Y.ravel()], 1)
else:
    b = open('/gpfs/kjhan/LagForce/glass_kepring_384.bin', 'rb').read(); m = n*n
    x0 = np.stack([np.frombuffer(b[20:20+8*m]), np.frombuffer(b[20+8*m:20+16*m])], 1)
X0 = x0 - C; R = np.hypot(*X0.T); keep = (R >= RIN) & (R <= ROUT); X0 = X0[keep]; R = R[keep]; th0 = np.arctan2(X0[:, 1], X0[:, 0])
Om = (R*R + EPS*EPS)**-0.75; dOm = -1.5*R*(R*R + EPS*EPS)**-1.75
def state(t):
    th = th0 + Om*t; x = np.stack([R*np.cos(th), R*np.sin(th)], 1)
    c, s = np.cos(Om*t), np.sin(Om*t); Rm = np.stack([np.stack([c, -s], -1), np.stack([s, c], -1)], -2)
    that = np.stack([-np.sin(th), np.cos(th)], 1); r0 = X0/R[:, None]
    F = Rm + (R*dOm*t)[:, None, None]*that[:, :, None]*r0[:, None, :]
    Fi = np.linalg.inv(F); return x, th, np.einsum('nji,njk->nik', Fi, Fi), np.linalg.cond(F)
def volumes(x, th, Gm, wedge0, K=K, half=np.radians(12)):
    # raster the wedge [wedge0-half, wedge0+half] (in the current frame), annulus margin
    a0, a1 = wedge0 - half, wedge0 + half
    rr = np.arange(RIN - 2*d, ROUT + 2*d, d/PX); out = []
    tree = cKDTree(x); V = np.zeros(len(x)); far = 0
    for r_ in rr:
        na = int((a1 - a0)*r_/(d/PX)) + 1; aa = a0 + (np.arange(na) + 0.5)*(a1 - a0)/na
        P = np.stack([r_*np.cos(aa), r_*np.sin(aa)], 1); w = (a1 - a0)*r_/na*(d/PX)
        dd, kk = tree.query(P, K); D = P[:, None, :] - x[kk]
        phi = np.einsum('pki,pkij,pkj->pk', D, Gm[kk], D); a = phi.argmin(1); pm = phi[np.arange(len(P)), a]
        far += int((a == K-1).sum()); ok = pm <= rho*rho
        np.add.at(V, kk[np.arange(len(P)), a][ok], w)
    return V, far
x, th, Gm, cond = state(0.0)
V0, far0 = volumes(x, th, np.tile(np.eye(2), (len(x), 1, 1)), 0.0, K=16, half=np.pi)   # material volume, full annulus
print('t=0 full-ring raster: sum V0 / annulus area = %.5f, far %d' % (V0.sum()/(np.pi*(ROUT**2 - RIN**2)), far0), flush=True)
for metric in ('voronoi', 'G1'):
    print('== %s IC=%s  K=%d' % (metric, ic, K), flush=True)
    for t in (20.0, 40.0, 60.0):
        x, th, Gm, cond = state(t)
        if metric == 'voronoi': Gm = np.tile(np.eye(2), (len(x), 1, 1))
        V, far = volumes(x, th, Gm, 0.0)
        sel = np.abs(np.angle(np.exp(1j*th))) < np.radians(6)
        for lo, hi in ((0.5, 0.6), (0.6, 0.8), (0.8, 1.2), (1.2, 2.0)):
            m = sel & (R >= lo) & (R < hi)
            if m.sum() < 5: continue
            e = V[m]/V0[m] - 1
            print('  t=%4.0f R %.1f-%.1f n=%4d rms dV/V %.4f max %.3f  cond(F) max %.0f  far %d' % (t, lo, hi, m.sum(), np.sqrt(np.mean(e*e)), np.abs(e).max(), cond[m].max(), far), flush=True)
