"""2D audit of design_astra.md sections 3-5 on the periodic Voronoi reference geometry (sub2d.Mesh).
xi_if = ln(h_if(c_i)/h_if(x_i)) with analytic gradient (centroid chain + normal rotation),
2-hop weighted-LSQ affine filter R_i, L_if = g_if R_i, Hugoniot C(q), F = -sum C r L^T, H_i = sum C r^2.
"""
import sys; sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks')
import numpy as np
from sub2d import Mesh, G, cross, R as Rop, L as Lop, total, U0

GAM = G


def cell_geom(mesh, i, x):
    cl = mesh.cells[i]; n = len(cl)
    a = np.array([vv for vv, _, _ in cl]); b = np.roll(a, -1, axis=0); P = x[i]
    crs = cross(a, b); V = 0.5*crs.sum(); c = (crs[:, None]*(a + b)).sum(0)/(6*V)
    W = 0.5*cross(a - P, b - P); B = 0.5*cross(a - c, b - c)
    # centroid chain: dc/dvertex_k (2x2) and vertex jacobians
    dV = np.array([0.5*Rop(a[(k+1) % n] - a[(k-1) % n]) for k in range(n)])
    Jc = []
    for k in range(n):
        ak, an, ap = a[k], a[(k+1) % n], a[(k-1) % n]
        ck, cp = cross(ak, an), cross(ap, ak)
        dJ = (np.outer(ak + an, Rop(an)) + ck*np.eye(2) + np.outer(ap + ak, Lop(ap)) + cp*np.eye(2))/6.0
        Jc.append((dJ - np.outer(c, dV[k]))/V)
    J = [[mesh.vjac(t, al) for al in range(3)] for (_, t, _) in cl]
    # dc/dX_a for tiled generator a (dict tiled index -> 2x2)
    dc = {}
    for k in range(n):
        t = cl[k][1]
        for al in range(3):
            dc[t[al]] = dc.get(t[al], 0) + Jc[k] @ J[k][al]
    return dict(n=n, a=a, b=b, P=P, V=V, c=c, W=W, B=B, dc=dc, nbr=mesh.nbr[i])


def face_xi_grad(mesh, i, x, cen):
    """returns per face: xi, h_x, W, rho-ratio B/W, and g (dict tiled index -> 2-vector), plus geometry"""
    gm = cell_geom(mesh, i, x); P = gm['P']; c = gm['c']
    faces = []
    for k, tj in enumerate(gm['nbr']):
        Xj = mesh.X[tj]; dvec = Xj - P; d = np.linalg.norm(dvec); e = dvec/d
        hx = 0.5*d; hc = 0.5*d - (c - P) @ e
        xi = np.log(hc/hx)
        Pe = np.eye(2) - np.outer(e, e)
        g = {}
        # direct terms
        gj = (0.5*e - Pe @ (c - P)/d)/hc - (0.5*e)/hx
        gi = (0.5*e + Pe @ (c - P)/d)/hc - (-0.5*e)/hx
        g[tj] = g.get(tj, 0) + gj
        g[cen + i] = g.get(cen + i, 0) + gi
        # centroid chain: -(dc/dX_a)^T e / hc
        for ta, M in gm['dc'].items():
            g[ta] = g.get(ta, 0) - (M.T @ e)/hc
        faces.append(dict(j=tj, xi=xi, hx=hx, W=gm['W'][k], B=gm['B'][k], d=d, A=np.linalg.norm(gm['b'][k] - gm['a'][k]), g=g))
    return faces, gm


def stencil(mesh, i, faces, gm, cen):
    """2-hop stencil with design weights; returns list of tiled indices, positions, weights"""
    N = mesh.N
    w = {}
    wdir = {}
    for f in faces:
        wdir[f['j']] = f['A']/f['d']**2
    hi = np.sqrt(gm['V'])
    for tj, wij in wdir.items():
        w[tj] = w.get(tj, 0) + wij
        j = mesh.idx[tj]; shift = mesh.X[tj] - mesh.x[j]
        # neighbours of j: mesh.nbr[j] are tiled indices relative to j's own central copy; shift them
        for tk in mesh.nbr[j]:
            pos = mesh.X[tk] + shift
            # find tiled index with this position (same particle k = idx[tk])
            k = mesh.idx[tk]
            # tile offset
            off = np.round((pos - mesh.x[k])/mesh.L).astype(int)
            if np.any(np.abs(off) > 1):
                continue
            tile = (off[0] + 1)*3 + (off[1] + 1)
            tkk = tile*N + k
            if tkk == cen + i:
                continue
            # A_jk, d_jk from j's face list: need face length between j and k -> recompute from j's cell
            Aj, dj = face_A_d(mesh, j, tk)
            w[tkk] = w.get(tkk, 0) + hi*wij*Aj/dj**2
    w[cen + i] = sum(wdir.values())
    ids = sorted(w); X = np.array([mesh.X[t] for t in ids]); W = np.array([w[t] for t in ids])
    return ids, X, W


_fa_cache = {}


def face_A_d(mesh, j, tk):
    key = (j, tk)
    if key in _fa_cache:
        return _fa_cache[key]
    cl = mesh.cells[j]; a = np.array([vv for vv, _, _ in cl]); b = np.roll(a, -1, axis=0)
    for k, tn in enumerate(mesh.nbr[j]):
        _fa_cache[(j, tn)] = (np.linalg.norm(b[k] - a[k]), np.linalg.norm(mesh.X[tn] - mesh.x[j]))
    return _fa_cache[key]


def fit_operator(X, Wt, x0):
    Xm = np.column_stack([np.ones(len(X)), X - x0])
    A = Xm.T @ (Wt[:, None]*Xm)
    # QR/SVD-robust solve
    Pm = Xm @ np.linalg.lstsq(A, Xm.T*Wt, rcond=None)[0]
    return Pm


def hg_eval(x, Lb, m, K, vfun, mesh=None, filt=True, hug=True, hgscale=1.0, coef=1.0, return_ops=False):
    """vfun(tiled_index) -> velocity 2-vector (allows unwrapped affine fields). Returns F (N,2), H (N,), list of (i, r, C, ids, LT)."""
    if mesh is None:
        mesh = Mesh(x, Lb)
    N = len(x); cen = 4*N
    _fa_cache.clear()
    F = np.zeros((N, 2)); H = np.zeros(N); ops = []
    for i in range(N):
        faces, gm = face_xi_grad(mesh, i, x, cen)
        ids, Xs, Ws = stencil(mesh, i, faces, gm, cen)
        pos = {t: k for k, t in enumerate(ids)}
        Pm = fit_operator(Xs, Ws, x[i]) if filt else np.zeros((len(ids), len(ids)))
        Rm = np.eye(len(ids)) - Pm
        vs = np.array([vfun(t) for t in ids])          # (S,2)
        Rv = Rm @ vs
        rho = m[i]/gm['V']
        for f in faces:
            gvec = np.zeros((len(ids), 2))
            for t, gg in f['g'].items():
                gvec[pos[t]] += gg
            r = (gvec*Rv).sum()
            LT = Rm.T @ gvec                            # (S,2)
            rhof = rho*f['B']/f['W']; cf = np.sqrt(GAM*K[i]*rhof**(GAM-1))
            h = f['hx']; W = f['W']; q = abs(r)
            if hug:
                b = hgscale*(GAM+1)/4*h
                C = coef*rhof*h*W*(b*q + np.sqrt(b*b*q*q + cf*cf))
            else:
                C = coef*rhof*cf*h*W
            H[i] += C*r*r
            for k, t in enumerate(ids):
                F[mesh.idx[t]] -= C*r*LT[k]
            if return_ops:
                ops.append((i, f, r, C, ids, LT, gvec, cf, rhof, h, W))
    return F, H, mesh, ops


def rayleigh(ops, hugscale=1.0, hug=True):
    Rt = 0.0
    for (i, f, r, C, ids, LT, gvec, cf, rhof, h, W) in ops:
        q = abs(r)
        if hug:
            b = hugscale*(GAM+1)/4*h
            Rt += rhof*h*W*(b*q**3/3 + ((b*b*q*q + cf*cf)**1.5 - cf**3)/(3*b*b))
        else:
            Rt += 0.5*rhof*cf*h*W*q*q
    return Rt


def glass(n, amp, seed=5, Lb=1.0):
    rng = np.random.default_rng(seed); h = Lb/n
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
    lat = np.stack([(i.ravel()+0.5)*h, (j.ravel()+0.5)*h], 1)
    x = (lat + rng.uniform(-amp, amp, lat.shape)*h) % Lb
    return x, h


def cell_volumes(mesh):
    return np.array([0.5*cross(np.array([vv for vv, _, _ in c]), np.roll(np.array([vv for vv, _, _ in c]), -1, 0)).sum() for c in mesh.cells])
