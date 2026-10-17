import numpy as np, arcgeom as g
rng = np.random.default_rng(1); R = 1.0; worst = 0
for trial in range(200):
    # random edge-like configuration: neighbours mostly on one side (vacuum on the other)
    k = rng.integers(3, 7); ang = rng.uniform(-0.6*np.pi, 0.6*np.pi, k); r = rng.uniform(0.9, 1.9, k)
    D = np.stack([r*np.cos(ang), r*np.sin(ang)], 1)
    pcs = g.clipped_cell(D, R)
    if pcs is None or not any(p[0] == 'arc' for p in pcs): continue
    for j in range(k):
        for comp in (0, 1):
            dD = np.zeros_like(D); dD[j, comp] = 1
            V, dV, c, dc, xs = g.tangent(pcs, D, dD, R)
            if not xs: continue
            h = 1e-6; Dp = D + h*dD; Dm = D - h*dD
            pp = g.clipped_cell(Dp, R); pm = g.clipped_cell(Dm, R)
            if [p[0] for p in pp] != [p[0] for p in pcs] or [p[0] for p in pm] != [p[0] for p in pcs]: continue
            Vp, cp, ap = g.arc_xi(pp, R); Vm, cm, am = g.arc_xi(pm, R)
            fdV = (Vp - Vm)/(2*h); fdc = (cp - cm)/(2*h); fdx = (ap[0][2] - am[0][2])/(2*h)
            e = max(abs(dV - fdV)/(abs(fdV) + 1e-3), np.abs(dc - fdc).max()/(np.abs(fdc).max() + 1e-3), abs(xs[0][1] - fdx)/(abs(fdx) + 1e-3))
            worst = max(worst, e)
print('max relative error analytic vs FD (V, c, arc xi) over random edge cells: %.2e' % worst)
