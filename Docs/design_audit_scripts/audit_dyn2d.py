"""Small 2D dynamic check of the design: periodic glass, cell-volume Lagrangian + face Q (QLAW=2 with d_r=d_u, face-subzone Z,
analytic Davis Z_eff + Z_H) + affine-filtered Hugoniot hourglass; RK4; K evolves with H_Q + H_hg.
Tests: pair approach at Mach M between two neighbours; checkerboard kick.  Reports min generator distance / h, dt, energy."""
import sys, time, numpy as np
sys.path.insert(0, '/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks')
from audit_hg2d import *
from sub2d import Mesh, cross, R as Rop


def cell_force(mesh, x, m, K):
    N = len(x); F = np.zeros_like(x); V = cell_volumes(mesh)
    rho = m/V; P = K*rho**GAM
    for i in range(N):
        cl = mesh.cells[i]; a = np.array([vv for vv, _, _ in cl]); n = len(a)
        dV = np.array([0.5*Rop(a[(k+1) % n] - a[(k-1) % n]) for k in range(n)])
        for k, (vv, t, _) in enumerate(cl):
            for al in range(3):
                F[mesh.idx[t[al]]] += P[i]*(mesh.vjac(t, al).T @ dV[k])
    U = np.sum(m*K*rho**(GAM-1)/(GAM-1))
    return F, V, P, U


def slope_ops(mesh, xm, allfaces, cen):
    """per cell: (ids tiled, positions, Bx, By) of the 2-hop weighted LSQ slope (same stencil as hg)"""
    import audit_hg2d as _ah; _ah._fa_cache.clear()
    out = []
    for i in range(len(xm)):
        gm_i = cell_geom(mesh, i, xm)
        ids, Xs, Ws = stencil(mesh, i, allfaces[i], gm_i, cen)
        Xm_ = np.column_stack([np.ones(len(Xs)), Xs - xm[i]])
        Bm = np.linalg.solve(Xm_.T @ (Ws[:, None]*Xm_), Xm_.T*Ws)
        out.append((ids, Xs, Bm[1], Bm[2]))
    return out


KAPPA_MAX = 1e4
FIT_STATS = dict(kmax=0.0, fail=0)


def fit_rows(Xs, Ws, x0, quad):
    """gradient rows (Dx, Dy) of the weighted LSQ fit (free intercept; linear or full quadratic, scaled by h)."""
    h = np.sqrt(np.sum(Ws*np.sum((Xs - x0)**2, 1))/np.sum(Ws)) if quad else 1.0
    sx = (Xs - x0)/h; wb = Ws/Ws.sum()
    cols = [np.ones(len(Xs)), sx[:, 0], sx[:, 1]]
    if quad: cols += [0.5*sx[:, 0]**2, sx[:, 0]*sx[:, 1], 0.5*sx[:, 1]**2]
    P = np.column_stack(cols); A = np.sqrt(wb)[:, None]*P
    sv = np.linalg.svd(A, compute_uv=False)
    kap = sv[0]/sv[-1] if sv[-1] > 0 else np.inf
    FIT_STATS['kmax'] = max(FIT_STATS['kmax'], kap)
    if not (sv[-1] > 64*2.2e-16*max(len(Xs), P.shape[1])*sv[0]) or kap > KAPPA_MAX:
        FIT_STATS['fail'] += 1
    B = np.linalg.lstsq(A, np.diag(np.sqrt(wb)), rcond=None)[0]     # coefficient extraction (ncol, S)
    return B[1]/h, B[2]/h


def bprime_ops_q(mesh, xm, allfaces, cen, mode):
    """(b') with gradient rows from: 'lin' (linear fit), 'quad' (quadratic fit), 'quadx' (quadratic fit excluding
    all images of both endpoints, per face). Same return format as bprime_ops."""
    import audit_hg2d as _ah; _ah._fa_cache.clear()
    N = len(xm); st = []
    for i in range(N):
        ids, Xs, Ws = stencil(mesh, i, allfaces[i], cell_geom(mesh, i, xm), cen)
        st.append((ids, Xs, Ws))
    rows = {}
    if mode in ('lin', 'quad'):
        for i in range(N):
            ids, Xs, Ws = st[i]; rows[i] = (ids, Xs, *fit_rows(Xs, Ws, xm[i], mode == 'quad'))
    res = []
    for i in range(N):
        for f in allfaces[i]:
            tj = f['j']; j = mesh.idx[tj]
            if j <= i: continue
            d = f['d']; e = (mesh.X[tj] - xm[i])/d; shift = mesh.X[tj] - xm[j]
            terms = [(i, xm[i], e.copy()), (j, mesh.X[tj], -e.copy())]
            for c, sh in ((i, 0*shift), (j, shift)):
                if mode == 'quadx':
                    ids, Xs, Ws = st[c]; keep = np.array([mesh.idx[t] not in (i, j) for t in ids])
                    ids = [t for t, k in zip(ids, keep) if k]; Xs = Xs[keep]; Ws = Ws[keep]
                    bx, by = fit_rows(Xs, Ws, xm[c], True)
                else:
                    ids, Xs, bx, by = rows[c]
                w = bx*e[0] + by*e[1]
                for k, t in enumerate(ids):
                    terms.append((mesh.idx[t], Xs[k] + sh, 0.5*d*w[k]*e))
            res.append((i, j, tj, f, terms))
    return res


def bprime_ops(mesh, xm, allfaces, sl, cen):
    """(b') residual operator per shared face, processed once: delta_f = (v_i - v_j).e + d e^T (G_i+G_j)/2 e (sign corrected vs audit 6).
    Returns list of (i, j, tj, terms) with terms = list of (particle, position(unwrapped), coeff 2-vector)."""
    N = len(xm); res = []
    for i in range(N):
        for f in allfaces[i]:
            tj = f['j']; j = mesh.idx[tj]
            if j < i or (j == i): continue
            d = f['d']; e = (mesh.X[tj] - xm[i])/d; shift = mesh.X[tj] - xm[j]
            terms = [(i, xm[i], e.copy()), (j, mesh.X[tj], -e.copy())]
            for c, sh in ((i, 0*shift), (j, shift)):
                ids, Xs, bx, by = sl[c]; w = bx*e[0] + by*e[1]
                for k, t in enumerate(ids):
                    terms.append((mesh.idx[t], Xs[k] + sh, 0.5*d*w[k]*e))   # + d e^T G e: affine (v_i - v_j).e = -d e^T A e
            res.append((i, j, tj, f, terms))
    return res


def rhs(x, v, K, m, Lb, opts):
    N = len(x); xm = x % Lb; mesh = Mesh(xm, Lb); cen = 4*N
    F, V, P, U = cell_force(mesh, xm, m, K)
    rho = m/V
    # face states per cell
    fst = {}
    allfaces = {}
    for i in range(N):
        faces, gm = face_xi_grad(mesh, i, xm, cen)
        for f in faces:
            rf = rho[i]*f['B']/f['W']; Pf = K[i]*rf**GAM
            fst[(i, f['j'])] = (rf, Pf, np.sqrt(GAM*Pf/rf), f['A'], f['d'], f['hx'])
        allfaces[i] = faces
    HQ = np.zeros(N); dtg = np.inf
    for (i, tj), (rL, pL, cl, A, d, hx) in fst.items():
        j = mesh.idx[tj]
        # j's state for the same face: tiled index of i as seen from j
        off = np.round((mesh.X[tj] - xm[j])/Lb).astype(int)
        ti = (( -off[0] + 1)*3 + (-off[1] + 1))*N + i
        rR, pR, cr, _, _, _ = fst.get((j, ti), (rL, pL, cl, A, d, hx))
        e = (mesh.X[tj] - xm[i])/d
        du = (v[i] - v[j]) @ e
        mu = rL*rR/(rL+rR); Pc = (rR*pL + rL*pR)/(rL+rR); cM = max(cl, cr)
        if du >= 0:
            Zeff = mu*(cM + du)
        else:
            Zeff = mu*cM if Pc + mu*cM*du > 0 else -Pc/du
        aH = (GAM+1)/8*abs(du); ZH = mu*(aH + np.sqrt(aH*aH + cM*cM))
        Z = max(Zeff, ZH)
        Q = Z*du if opts.get('qon', True) else 0.0
        F[i] -= A*Q*e            # j side handled when visiting (j, ti)
        HQ[i] += 0.5*A*Q*du
        dtg = min(dtg, 0.25*hx/(cl + abs(du)))
    HB = np.zeros(N)
    if opts.get('bprime', False):
        lamB = np.zeros(N)
        bmode = opts.get('bmode', 'lin')
        for (i, j, tj, f, terms) in bprime_ops_q(mesh, xm, allfaces, cen, bmode):
            dl = sum(c @ v[a] for a, _, c in terms)
            off = np.round((mesh.X[tj] - xm[j])/Lb).astype(int); ti = ((-off[0] + 1)*3 + (-off[1] + 1))*N + i
            rL, pL, cl, A, d, hx = fst[(i, tj)]; rR, pR, cr, _, _, _ = fst.get((j, ti), fst[(i, tj)])
            mu = rL*rR/(rL + rR); cM = max(cl, cr); aH = (GAM + 1)/8*abs(dl)
            Z = max(mu*(cM + max(dl, 0.0)), mu*(aH + np.sqrt(aH*aH + cM*cM)))
            if opts.get('bprime') == 2: Z = mu*aH          # Astra design 7: no acoustic term
            C = A*Z
            for a, _, c in terms: F[a] -= C*dl*c
            HB[i] += 0.5*C*dl*dl; HB[j] += 0.5*C*dl*dl
            al = sum((c @ c)/m[a] for a, _, c in terms)
            for a in set(a for a, _, _ in terms): lamB[a] += C*al   # row-sum bound over overlapping modes
        dtg = min(dtg, 1.0/max(lamB.max(), 1e-300))
    if opts.get('bulk', False):
        # Astra design 5: eta_i a_i^2 on a_i = min(0, tr G_i), G_i from the 2-hop weighted LSQ affine fit (same stencil as hg)
        import audit_hg2d as _ah; _ah._fa_cache.clear()
        for i in range(N):
            faces_i = allfaces[i]; gm_i = cell_geom(mesh, i, xm)
            ids, Xs, Ws = stencil(mesh, i, faces_i, gm_i, cen)
            Xm_ = np.column_stack([np.ones(len(Xs)), Xs - xm[i]])
            Bm = np.linalg.solve(Xm_.T @ (Ws[:, None]*Xm_), Xm_.T*Ws)     # (3,S)
            vs = np.array([v[mesh.idx[t]] for t in ids])
            dv = Bm[1] @ vs[:, 0] + Bm[2] @ vs[:, 1]
            if dv >= 0: continue
            hh = np.sqrt(V[i]); ci = np.sqrt(GAM*P[i]/rho[i])
            eta = m[i]*hh*(ci + (GAM+1)/4*hh*abs(dv))
            for k, t in enumerate(ids):
                F[mesh.idx[t], 0] -= eta*dv*Bm[1, k]; F[mesh.idx[t], 1] -= eta*dv*Bm[2, k]
            HB[i] += eta*dv*dv
            al = sum((Bm[1, k]**2 + Bm[2, k]**2)/m[mesh.idx[t]] for k, t in enumerate(ids))
            dtg = min(dtg, 1.0/(eta*al))
    HQ = HQ + HB
    H = np.zeros(N)
    if opts.get('hg', True):
        Fh, H, _, _ = hg_eval(xm, Lb, m, K, lambda t: v[mesh.idx[t]], mesh=mesh, filt=opts.get('filt', True), hug=opts.get('hug', True))
        F += Fh
    G = m*rho**(GAM-1)/(GAM-1)
    return v, F/m[:, None], (HQ + H)/G, dict(dtg=dtg, HQ=HQ.sum(), H=H.sum(), U=U, V=V)


def mind(x, Lb):
    d = x[:, None, :] - x[None, :, :]; d -= Lb*np.round(d/Lb)
    r = np.sqrt((d**2).sum(-1)); np.fill_diagonal(r, np.inf); return r.min()


def run(name, x0, v0, m, K0, Lb, h, T, opts, c0):
    x = x0.copy(); v = v0.copy(); K = K0.copy(); t = 0.0
    _, _, _, a = rhs(x, v, K, m, Lb, opts)
    E0 = a['U'] + 0.5*np.sum(m[:, None]*v*v); KE0 = 0.5*np.sum(m[:, None]*v*v)
    md = [mind(x % Lb, Lb)/h]; dts = []; cum = dict(HQ=0.0, H=0.0); n = 0; t0 = time.time(); ok = 'ok'
    while t < T:
        _, _, _, a = rhs(x, v, K, m, Lb, opts)
        dt = min(a['dtg'], T - t); dts.append(dt)
        if dt < 1e-9*h/c0:
            ok = f'dt collapse {dt:.1e}'; break
        ks = []
        xx, vv, KK = x, v, K
        for cfac in (0.0, 0.5, 0.5, 1.0):
            if ks:
                xx = x + cfac*dt*ks[-1][0]; vv = v + cfac*dt*ks[-1][1]; KK = K + cfac*dt*ks[-1][2]
            ks.append(rhs(xx, vv, KK, m, Lb, opts))
        w = (1, 2, 2, 1)
        x = x + dt/6*sum(wi*k[0] for wi, k in zip(w, ks)); v = v + dt/6*sum(wi*k[1] for wi, k in zip(w, ks)); K = K + dt/6*sum(wi*k[2] for wi, k in zip(w, ks))
        for wi, k in zip(w, ks):
            cum['HQ'] += dt/6*wi*k[3]['HQ']; cum['H'] += dt/6*wi*k[3]['H']
        t += dt; n += 1
        md.append(mind(x % Lb, Lb)/h)
        if np.any(K <= 0) or md[-1] < 1e-3:
            ok = 'COLLISION'; break
        if time.time() - t0 > 5400:
            ok = f'time limit at t={t/(h/c0):.2f} h/c'; break
    _, _, _, a = rhs(x, v, K, m, Lb, opts)
    E = a['U'] + 0.5*np.sum(m[:, None]*v*v); KE = 0.5*np.sum(m[:, None]*v*v)
    dts = np.array(dts)
    print(f"{name} {opts}: {ok} steps={n} t={t/(h/c0):.2f}h/c min dist/h: t0 {md[0]:.3f} -> min {min(md):.3f} -> end {md[-1]:.3f}; "
          f"dt/(h/c) first/min/last = {dts[0]/(h/c0):.2e}/{dts.min()/(h/c0):.2e}/{dts[-1]/(h/c0):.2e}; KE/KE0={KE/KE0:.3f} dE/E0={(E-E0)/abs(E0):+.1e} "
          f"HQ/KE0={cum['HQ']/KE0:.3f} Hhg/KE0={cum['H']/KE0:.3f} [{time.time()-t0:.0f}s]", flush=True)


OPTS = dict(bulk=(dict(qon=False, bulk=True),),
            bp=(dict(qon=False, bulk=True, bprime=True), dict(qon=False, bulk=True, bprime=True, hg=False)),
            bp2=(dict(qon=False, bulk=True, bprime=2), dict(qon=False, bulk=True, bprime=2, hg=False)),
            bpq=(dict(qon=False, bulk=True, bprime=1, bmode='quad'),),
            bpqx=(dict(qon=False, bulk=True, bprime=1, bmode='quadx'),),
            bpl=(dict(qon=False, bulk=True, bprime=1, bmode='lin'),))

if __name__ == '__main__':
    which = sys.argv[1]
    n = 8; Lb = 1.0; amp = 0.2
    x0, h = glass(n, amp, 5); N = len(x0); m = np.full(N, h*h)
    mesh = Mesh(x0, Lb); V = cell_volumes(mesh); K0 = 1.0/(m/V)**GAM; c0 = np.sqrt(GAM*1.0/(m/V).mean())
    cen = 4*N
    i0 = N//2 + n//2; faces0, _ = face_xi_grad(mesh, i0, x0, cen); tj = faces0[0]['j']; j0 = mesh.idx[tj]
    e = mesh.X[tj] - x0[i0]; e /= np.linalg.norm(e)
    if which.startswith('pair'):
        M = float(which[4:])
        v0 = np.zeros_like(x0); v0[i0] = 0.5*M*c0*e; v0[j0] = -0.5*M*c0*e
        T = 4.0*h/c0
        for opts in (OPTS[sys.argv[2]] if len(sys.argv) > 2 else (dict(), dict(hg=False), dict(qon=False))):
            run(f"PAIR2D M={M}", x0, v0, m, K0, Lb, h, T, opts, c0)
    elif which == 'chk':
        ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing='ij'); chk = (-1.0)**(ii.ravel() + jj.ravel())
        v0 = 0.3*c0*np.stack([0*chk, chk], 1)
        for opts in (OPTS[sys.argv[2]] if len(sys.argv) > 2 else (dict(), dict(hg=False), dict(qon=False))):
            run("CHK2D Mach0.3", x0, v0, m, K0, Lb, h, 4.0*h/c0, opts, c0)
