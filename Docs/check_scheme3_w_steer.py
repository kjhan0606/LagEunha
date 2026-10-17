"""Newton (volume-constrained W) with generator steering x' = v + (c - x)/tau. V^G evolves with the fluid v only."""
import numpy as np
src = open('check_scheme3_w_drift.py').read().split('dt, nstep = 0.01, 60')[0]
__file__ = 'check_scheme3_w_drift.py'
exec(src)

def centroids(X, W):
    Xm = X % L; C = np.zeros((N, 2))
    for i in range(N):
        xi = Xm[i]; h = 1.5 * L
        poly = [(xi + np.array(p) * h, -1) for p in [(-1, -1), (1, -1), (1, 1), (-1, 1)]]
        for s, sh in enumerate(shifts):
            for j in range(N):
                if j == i and s == 4: continue
                xj = Xm[j] + sh
                if np.sum((xj - xi) ** 2) > (3.5 * dmean) ** 2: continue
                poly = clip(poly, 2 * (xj - xi), xj @ xj - xi @ xi + W[i] - W[j], (j, s))
        P = np.array([p for p, _ in poly]); x, y = P[:, 0], P[:, 1]
        cr = x * np.roll(y, -1) - np.roll(x, -1) * y; A = 0.5 * cr.sum()
        C[i] = [np.sum((x + np.roll(x, -1)) * cr) / (6 * A), np.sum((y + np.roll(y, -1)) * cr) / (6 * A)]
    D = C - Xm
    return Xm + D  # centroid in the generator's image frame

dt, nstep = 0.01, 60
for tau in (0.05,):
    X = X0.copy(); W = np.zeros(N); VG = cells(X0, W)[0].copy(); it = 0
    print(f"tau={tau}:  t   max|V-VG|/Vbar  max|dW|/d2  max|x-c|/dbar  newton_it")
    for k in range(nstep + 1):
        t = k * dt
        if k % 10 == 0:
            e = np.max(np.abs(cells(X % L, W)[0] - VG)) / VG.mean()
            c = centroids(X, W); off = np.max(np.linalg.norm(c - X % L, axis=1)) / dmean
            print(f"        {t:4.2f}  {e:12.3e}  {maxratio(X, W):9.3f}  {off:11.3f}  {it}")
        if k == nstep: break
        def st(X, W, VGs):
            W, it = newton(X, W, VGs)
            v = vel(X % L, t)
            xd = v + (centroids(X, W) - X % L) / tau
            return xd, vg_rate(X, W, v), W, it
        a = st(X, W, VG); b = st(X + dt/2*a[0], a[2], VG + dt/2*a[1])
        c_ = st(X + dt/2*b[0], b[2], VG + dt/2*b[1]); d = st(X + dt*c_[0], c_[2], VG + dt*c_[1])
        X = X + dt/6*(a[0]+2*b[0]+2*c_[0]+d[0]); VG = VG + dt/6*(a[1]+2*b[1]+2*c_[1]+d[1])
        W, it = newton(X, d[2], VG)
