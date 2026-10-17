"""Reversible Kepler ring dynamics (no Q, no HG) with disc-clipped or common-polygon-closed Voronoi cells.
Fixed masses/entropy from the initial cell areas (Sigma0 = 1, P0 = 1e-6), GM = 1, eps = 0.01, RK4 fixed dt."""
import sys, time, numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce/PolyBnd'); import polygeom as pg
mode = sys.argv[1]; d = float(sys.argv[2]); dt = float(sys.argv[3]); T = float(sys.argv[4])
G_ = 5/3; P0 = 1e-6; EPS = 0.01; RIN, ROUT = 0.5, 2.0
P = []; lab = []
K = int(round((ROUT - RIN)/d))
for k in range(K + 1 if mode == 'poly' else K):
    r = RIN + (k if mode == 'poly' else k + 0.5)*d; n = int(round(2*np.pi*r/d)); th = -np.pi + (np.arange(n) + 0.5)*2*np.pi/n
    P.append(np.stack([r*np.cos(th), r*np.sin(th)], 1)); lab += [k]*n
X = np.concatenate(P); lab = np.array(lab); N = len(X)
kw = dict(inner=np.where(lab == 0)[0][::-1], outer=np.where(lab == lab.max())[0]) if mode == 'poly' else dict(Rc=d)
g0 = pg.Geom(X, mode, **kw); m = g0.V.copy(); Kent = P0/1.0**G_
R = np.hypot(*X.T); Om = (R*R + EPS*EPS)**-0.75; Vv = np.stack([-Om*X[:, 1], Om*X[:, 0]], 1)
def acc(X):
    g = pg.Geom(X, mode, **kw); G = g.grad(); Pc = Kent*(m/g.V)**G_
    F = np.zeros_like(X)
    for a in range(N):
        for k, v in G[a].items(): F[k] += Pc[a]*v
    r2 = (X*X).sum(1) + EPS*EPS
    return F/m[:, None] - X/r2[:, None]**1.5, g.V
def energy(X, Vv, V):
    return 0.5*(m[:, None]*Vv*Vv).sum() - (m/np.sqrt((X*X).sum(1) + EPS*EPS)).sum() + (Kent*m**G_*V**(1 - G_)/(G_ - 1)).sum()
th0 = np.arctan2(X[:, 1], X[:, 0]); a0, V0 = acc(X); E0 = energy(X, Vv, V0); U0 = (Kent*m**G_*V0**(1 - G_)/(G_ - 1)).sum()
print('mode %s d %.3f N %d edge ring sizes %s  E0 %.6e U0 %.3e' % (mode, d, N, [int((lab == j).sum()) for j in (0, 1, lab.max())], E0, U0), flush=True)
t = 0.0; nstep = int(round(T/dt)); t0 = time.time()
for s in range(1, nstep + 1):
    k1x = Vv; k1v, _ = acc(X)
    k2x = Vv + 0.5*dt*k1v; k2v, _ = acc(X + 0.5*dt*k1x)
    k3x = Vv + 0.5*dt*k2v; k3v, _ = acc(X + 0.5*dt*k2x)
    k4x = Vv + dt*k3v; k4v, Vc = acc(X + dt*k3x)
    X = X + dt/6*(k1x + 2*k2x + 2*k3x + k4x); Vv = Vv + dt/6*(k1v + 2*k2v + 2*k3v + k4v); t += dt
    if s % int(round(2.0/dt)) == 0:
        _, V = acc(X); E = energy(X, Vv, V)
        R = np.hypot(*X.T); vr = (X*Vv).sum(1)/R; vk = R*(R*R + EPS*EPS)**-0.75; th = np.arctan2(X[:, 1], X[:, 0])
        out = []
        for j in (0, 1, lab.max()):
            sel = np.where(lab == j)[0]; o = sel[np.argsort(th0[sel])]; n = len(o)
            dth = np.angle(np.exp(1j*(th[o] - th0[o]))); dth -= np.angle(np.mean(np.exp(1j*dth)))
            F = np.abs(np.fft.rfft(dth))*2/n/(2*np.pi/n); F[:2] = 0; q = np.argmax(F)
            out.append('k%d vr/vK %.1e m%d %.3f' % (j, np.sqrt(np.mean((vr[sel]/vk[sel])**2)), q, F[q]))
        print('t %5.1f  (E-E0)/U0 %+.2e  | %s  [%.0fs]' % (t, (E - E0)/U0, ' | '.join(out), time.time() - t0), flush=True)
