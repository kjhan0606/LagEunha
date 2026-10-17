"""Integrated Wdot (differentiated volume constraint) vs algebraic Newton (semi-discrete OT).
Kinematic toy: generators advected by a smooth periodic shear+compression field.
Target volumes V^G: dV^G/dt = sum_f A_f * (u_i + u_j)/2 (face-mean normal speed, a stand-in for u*)."""
import numpy as np, sys
sys.argv = ['x']
exec(open(__file__.replace('check_scheme3_w_drift.py', 'check_scheme3_formulation.py')).read().split('V0, F0 = cells(X0, W0)')[0].replace('(2.5 * dmean)', '(3.5 * dmean)'))


def vel(X, t):
    x, y = X[:, 0], X[:, 1]
    tp = 2 * np.pi
    return np.stack([0.6 * np.sin(tp * y) + 0.25 * np.sin(tp * x),
                     0.3 * np.sin(tp * x + 1.0) + 0.15 * np.sin(tp * y + 0.5)], 1)


def geom(X, W):
    V, F = cells(X % L, W)
    Lm = np.zeros((N, N)); Jx = np.zeros((N, N, 2)); fl = []
    Xm = X % L
    for (i, j, s), (A, c) in F.items():
        xj = Xm[j] + shifts[s]; d = np.linalg.norm(xj - Xm[i]); n = (xj - Xm[i]) / d
        Lm[i, i] += A / (2 * d); Lm[i, j] -= A / (2 * d)
        Jx[i, i] += A * (c - Xm[i]) / d; Jx[i, j] -= A * (c - xj) / d
        fl.append((i, j, A, n, d))
    return V, Lm, Jx, fl


def vg_rate(X, W, v):
    V, Lm, Jx, fl = geom(X, W)
    r = np.zeros(N)
    for i, j, A, n, d in fl:
        r[i] += A * 0.5 * (v[i] + v[j]) @ n
    return r


def wdot(X, W, v, VGdot):
    V, Lm, Jx, fl = geom(X, W)
    rhs = VGdot - np.einsum('ika,ka->i', Jx, v)
    wd = np.linalg.lstsq(Lm, rhs, rcond=None)[0]
    return wd - wd.mean()


def newton(X, W, VG, tol=1e-12, itmax=30):
    for it in range(itmax):
        V, Lm, _, _ = geom(X, W)
        e = V - VG
        if np.max(np.abs(e)) < tol * VG.mean():
            return W, it
        dW = np.linalg.lstsq(Lm, -e, rcond=None)[0]
        lam = 1.0
        while True:  # damped: keep all cells nonempty
            try:
                Vn = geom(X, W + lam * dW)[0]
                if np.all(Vn > 0) and np.max(np.abs(Vn - VG)) < np.max(np.abs(e)):
                    break
            except RuntimeError:
                pass
            lam *= 0.5
            if lam < 1e-4:
                return W, -1
        W = W + lam * dW
    return W, itmax


def maxratio(X, W):
    Xm = X % L; r = 0
    for (i, j, s), _ in cells(Xm, W)[1].items():
        d2 = np.sum((Xm[j] + shifts[s] - Xm[i]) ** 2); r = max(r, abs(W[i] - W[j]) / d2)
    return r


dt, nstep = 0.01, 60
XA = X0.copy(); WA = np.zeros(N); VGA = cells(X0, WA)[0].copy()
XB = X0.copy(); WB = np.zeros(N); VGB = VGA.copy()
print(" t     | integrated Wdot: max|V-VG|/Vbar  max|dW|/d2 | Newton: max|V-VG|/Vbar  max|dW|/d2  iters")
for k in range(nstep + 1):
    t = k * dt
    if k % 10 == 0:
        try:
            eA = np.max(np.abs(cells(XA % L, WA)[0] - VGA)) / VGA.mean(); rA = maxratio(XA, WA)
        except RuntimeError:
            eA = rA = float('nan')
        eB = np.max(np.abs(cells(XB % L, WB)[0] - VGB)) / VGB.mean(); rB = maxratio(XB, WB)
        print(f" {t:4.2f}  | {eA:12.3e} {rA:14.3f}           | {eB:12.3e} {rB:12.3f}  {itB if k else 0}")
    if k == nstep:
        break
    # A: RK4 on (X, W, VG) with Wdot from the differentiated constraint
    def rhs(X, W):
        v = vel(X % L, t); vg = vg_rate(X, W, v)
        return v, wdot(X, W, v, vg), vg
    try:
        k1 = rhs(XA, WA); k2 = rhs(XA + dt/2*k1[0], WA + dt/2*k1[1])
        k3 = rhs(XA + dt/2*k2[0], WA + dt/2*k2[1]); k4 = rhs(XA + dt*k3[0], WA + dt*k3[1])
        XA = XA + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); WA = WA + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        VGA = VGA + dt/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    except RuntimeError:
        pass
    # B: RK4 on (X, VG); W from Newton at every stage
    def rhsB(X, W):
        W, it = newton(X, W, None if False else VGcur)
        v = vel(X % L, t); return v, vg_rate(X, W, v), W, it
    VGcur = VGB; a = rhsB(XB, WB)
    VGcur = VGB + dt/2*a[1]; b = rhsB(XB + dt/2*a[0], a[2])
    VGcur = VGB + dt/2*b[1]; c = rhsB(XB + dt/2*b[0], b[2])
    VGcur = VGB + dt*c[1]; d = rhsB(XB + dt*c[0], c[2])
    XB = XB + dt/6*(a[0]+2*b[0]+2*c[0]+d[0]); VGB = VGB + dt/6*(a[1]+2*b[1]+2*c[1]+d[1])
    VGcur = VGB; WB, itB = newton(XB, d[2], VGB)
