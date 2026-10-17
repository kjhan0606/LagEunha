import numpy as np, polygeom as pg
def rings(d, RIN=0.5, ROUT=2.0, jitter=0.0, seed=0):
    rng = np.random.default_rng(seed); P = []; lab = []
    K = int(round((ROUT - RIN)/d))
    for k in range(K + 1):
        r = RIN + k*d; n = int(round(2*np.pi*r/d)); th = (np.arange(n) + 0.5*(k % 2))*2*np.pi/n
        rr = r + (jitter*d*rng.uniform(-1, 1, n) if 0 < k < K else 0)
        P.append(np.stack([rr*np.cos(th), rr*np.sin(th)], 1)); lab += [k]*n
    return np.concatenate(P), np.array(lab), K
X, lab, K = rings(0.25, jitter=0.2)
inner = np.where(lab == 0)[0][::-1]; outer = np.where(lab == K)[0]
for mode in ('poly', 'disc'):
    kw = dict(inner=inner, outer=outer) if mode == 'poly' else dict(Rc=0.25)
    g0 = pg.Geom(X, mode, **kw); G = g0.grad()
    print(mode, 'N', len(X), 'sum V', g0.V.sum(), '(annulus area %.4f)' % (np.pi*(4 - 0.25)))
    rng = np.random.default_rng(1); worst = 0
    for trial in range(10):
        k = rng.integers(len(X)); comp = rng.integers(2); h = 1e-6
        Xp = X.copy(); Xp[k, comp] += h; Xm = X.copy(); Xm[k, comp] -= h
        Vp = pg.Geom(Xp, mode, **kw).V; Vm = pg.Geom(Xm, mode, **kw).V
        fd = (Vp - Vm)/(2*h)
        an = np.array([G[a].get(k, np.zeros(2))[comp] if k in G[a] else 0.0 for a in range(len(X))])
        worst = max(worst, np.abs(fd - an).max()/max(np.abs(fd).max(), 1e-12))
    print('   max |dV/dx analytic - FD| / max|FD| over 10 random generator components: %.2e' % worst)
