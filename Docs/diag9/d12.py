"""Astra design 9, D1/D2: prescribed trajectories, no dissipation. Voronoi volume vs material volume,
adiabatic pressure noise, and the reversible (cell-volume) force with adiabatic vs analytic pressure."""
import sys; sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks'); sys.path.insert(0, '/tmp/lageunha-push/Docs/design_audit_scripts')
import numpy as np
from sub2d import Mesh, cross, R as Rop
G = 5/3
def vols_force(x, P):
    mesh = Mesh(x % 1.0, 1.0); N = len(x); V = np.zeros(N); F = np.zeros_like(x)
    for i in range(N):
        cl = mesh.cells[i]; a = np.array([vv for vv, _, _ in cl]); n = len(a)
        V[i] = 0.5*cross(a, np.roll(a, -1, 0)).sum()
        dV = np.array([0.5*Rop(a[(k+1) % n] - a[(k-1) % n]) for k in range(n)])
        for k, (vv, t, _) in enumerate(cl):
            for al in range(3):
                F[mesh.idx[t[al]]] += P[i]*(mesh.vjac(t, al).T @ dV[k])
    return V, F
def lattice(n, amp=0.0, seed=1):
    h = 1/n; i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
    x = np.stack([(i.ravel()+.5)*h, (j.ravel()+.5)*h], 1)
    if amp: x = (x + np.random.default_rng(seed).uniform(-amp, amp, x.shape)*h) % 1.0
    return x, h
# ---------- D1: affine shear, interior cells
print('D1 affine shear x += s t (y-1/2), s=1; interior |y-1/2|<0.25; dissipation-free, P0=1')
for lab, amp in (('lattice', 0.0), ('glass .25', 0.25)):
    for n in (16, 32):
        x0, h = lattice(n, amp); V0, _ = vols_force(x0, np.zeros(len(x0))); m = V0.copy()
        out = []
        for st in (0.1, 0.25, 0.5):
            x = x0.copy(); x[:, 0] = (x0[:, 0] + st*(x0[:, 1] - 0.5)) % 1.0
            V, _ = vols_force(x, np.zeros(len(x)))
            msk = np.abs(x0[:, 1] - 0.5) < 0.25
            Pad = (V0/V)**G
            _, Fad = vols_force(x, Pad); _, F1 = vols_force(x, np.ones(len(x)))
            a_ad = np.hypot(*(Fad/m[:, None])[msk].T)
            out.append('st=%.2f rms dV/V %.3f  rms dP/P %.3f  rms|a| (P_ad) %.2e (exact 0; P=1 gives %.1e)' %
                       (st, np.sqrt(np.mean((V/V0 - 1)[msk]**2)), np.sqrt(np.mean((Pad - 1)[msk]**2)), np.sqrt(np.mean(a_ad**2)),
                        np.sqrt(np.mean(np.hypot(*(F1/m[:, None])[msk].T)**2))))
        print(' %-10s n=%d: ' % (lab, n) + ' | '.join(out))
# ---------- D2: prescribed Gresho rotation
def om(r): return np.where(r < 0.2, 5.0, np.where(r < 0.4, 2/r - 5, 0.0))
def Pg(r): return np.where(r < 0.2, 5 + 12.5*r*r, np.where(r < 0.4, 9 + 12.5*r*r - 20*r + 4*np.log(5*np.maximum(r, 1e-9)), 3 + 4*np.log(2)))
def vphi(r): return np.where(r < 0.2, 5*r, np.where(r < 0.4, 2 - 5*r, 0.0))
print('D2 prescribed Gresho rotation (theta += omega(r) t), rho0=1, P_ad = P0 (V0/V)^g vs P_analytic(x)')
for n in (32, 64):
    x0, h = lattice(n); X0 = x0 - 0.5; r = np.hypot(*X0.T); th0 = np.arctan2(X0[:, 1], X0[:, 0])
    V0, _ = vols_force(x0, np.zeros(len(x0))); m = V0.copy(); P0 = Pg(r)
    for t in (0.0, 0.25, 0.5, 1.0):
        th = th0 + om(r)*t; x = 0.5 + np.stack([r*np.cos(th), r*np.sin(th)], 1)
        V, _ = vols_force(x, np.zeros(len(x)))
        Pad = P0*(V0/V)**G; Pan = Pg(r)
        _, Fad = vols_force(x, Pad); _, Fan = vols_force(x, Pan)
        rh = np.stack([np.cos(th), np.sin(th)], 1); astar = -(vphi(r)**2/np.maximum(r, 1e-9))[:, None]*rh
        msk = r < 0.45; sm = msk & (np.abs(r - 0.2) > 2*h) & (np.abs(r - 0.4) > 2*h)
        e_ad = np.hypot(*((Fad/m[:, None]) - astar).T); e_an = np.hypot(*((Fan/m[:, None]) - astar).T); sc = np.sqrt(np.mean(np.hypot(*astar[msk].T)**2))
        print(' n=%d t=%.2f rms dV/V %.4f | rms dP_ad/P %.2e | accel err rms/|a*|: P_ad %.3f (smooth %.3f)  P_analytic %.3f (smooth %.3f)' %
              (n, t, np.sqrt(np.mean((V/V0 - 1)[msk]**2)), np.sqrt(np.mean(((Pad - Pan)/Pan)[msk]**2)),
               np.sqrt(np.mean(e_ad[msk]**2))/sc, np.sqrt(np.mean(e_ad[sm]**2))/sc, np.sqrt(np.mean(e_an[msk]**2))/sc, np.sqrt(np.mean(e_an[sm]**2))/sc), flush=True)
