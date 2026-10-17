"""2D periodic prototype: formulation (G) (Godunov volume + OT weights, particle = generator, no steering)
vs scheme (1) GFS-like (Voronoi, W=0, work with geometric face velocity incl. rotation).
Same first-order HLLC-type Riemann (P*, u*), same RK4. Usage: python proto_G.py gresho G|1 N tend"""
import sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import ConvexHull

GAM = 5.0 / 3.0
LB = 1.0
GRAV = None       # function X -> (acc, phi)
PFLOOR = None


def diagram(X, W, margin):
    N = len(X)
    offs = [(0, 0)] + [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
    P, idx = [], []
    for a, b in offs:
        Y = X + np.array([a, b]) * LB
        m = np.ones(N, bool) if (a, b) == (0, 0) else \
            (Y[:, 0] > -margin) & (Y[:, 0] < LB + margin) & (Y[:, 1] > -margin) & (Y[:, 1] < LB + margin)
        P.append(Y[m]); idx.append(np.nonzero(m)[0])
    P = np.vstack(P); idx = np.concatenate(idx); Wa = W[idx]
    z = np.sum(P * P, 1) - Wa
    hull = ConvexHull(np.c_[P, z])
    T = hull.simplices[hull.equations[:, 2] < 0]
    pa, pb, pc = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    M = np.stack([2 * (pb - pa), 2 * (pc - pa)], 1)             # (T,2,2)
    r = np.stack([z[T[:, 1]] - z[T[:, 0]], z[T[:, 2]] - z[T[:, 0]]], 1)
    C = np.linalg.solve(M, r[..., None])[..., 0]                  # power centres
    E = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    tid = np.tile(np.arange(len(T)), 3)
    E.sort(1)
    key = E[:, 0].astype(np.int64) * len(P) + E[:, 1]
    o = np.argsort(key, kind='stable'); key, E, tid = key[o], E[o], tid[o]
    dup = np.nonzero(key[1:] == key[:-1])[0]
    a, b = E[dup, 0], E[dup, 1]
    keep = a < N                                                   # owner real (reals have smallest ids)
    a, b, t1, t2 = a[keep], b[keep], tid[dup][keep], tid[dup + 1][keep]
    c1, c2 = C[t1], C[t2]
    A = np.linalg.norm(c1 - c2, axis=1)
    cf = 0.5 * (c1 + c2)
    dv = P[b] - P[a]; d = np.linalg.norm(dv, axis=1); n = dv / d[:, None]
    return dict(a=a, b=b, jb=idx[b], breal=b < N, A=A, cf=cf, d=d, n=n, Pa=P[a], N=N, W=W)


def volumes(F):
    N, a, jb, A, d = F['N'], F['a'], F['jb'], F['A'], F['d']
    W = F['W']
    s = (d * d + W[a] - W[jb]) / (2 * d)
    V = np.zeros(N)
    np.add.at(V, a, 0.5 * A * s)
    br = F['breal']
    np.add.at(V, jb[br], 0.5 * A[br] * (d[br] - s[br]))
    return V, s


def hessian(F):
    N, a, jb, A, d, br = F['N'], F['a'], F['jb'], F['A'], F['d'], F['breal']
    w = A / (2 * d)
    rows = [a, a]; cols = [a, jb]; vals = [w, -w]
    rows += [jb[br], jb[br]]; cols += [jb[br], a[br]]; vals += [w[br], -w[br]]
    H = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N)).tocsr()
    return H


H_TAUV = 0.3
H_TAUC = 0.02
H_R0 = 0.5
H_DV = 0.1
H_KB = 10.0
H_TAUB = 0.02
CAP = None          # None: plain (G). Else dict(delta=, kappa=)


def pairs_from(F):
    """Fixed cap pair set (Voronoi/current neighbours at solve start); image pairs appear twice -> weight 1/2."""
    return dict(a=F['a'].copy(), jb=F['jb'].copy(), d=F['d'].copy(), w=np.where(F['breal'], 1.0, 0.5))


def cap_terms(Pr, W):
    dW = W[Pr['a']] - W[Pr['jb']]
    c = (1 - CAP['delta']) * Pr['d']**2
    ex = np.maximum(0.0, np.abs(dW) - c)
    t = CAP['kappa'] * Pr['w'] * ex * np.sign(dW)   # transfer: a loses t, jb gains t
    return t, ex > 0, ex


def phi(F, W, VG, Pr):
    """Concave OT dual minus cap penalty (to be maximised)."""
    a, jb, br, A, d = F['a'], F['jb'], F['breal'], F['A'], F['d']
    s_ = (d * d + W[a] - W[jb]) / (2 * d)
    # second moments of triangles (generator, face) : signed area 1/2 A s ; M2 = area/6 (|p|^2+|q|^2+p.q)
    # face endpoints relative to generators: centroid cf, half-length along tangent
    n = F['n']; tg = np.stack([-n[:, 1], n[:, 0]], 1)
    pa = F['cf'] - F['Pa']
    p1 = pa + 0.5 * A[:, None] * tg; p2 = pa - 0.5 * A[:, None] * tg
    m2a = 0.5 * A * s_ / 6 * (np.sum(p1 * p1, 1) + np.sum(p2 * p2, 1) + np.sum(p1 * p2, 1))
    pb = pa - n * d[:, None]
    q1 = pb + 0.5 * A[:, None] * tg; q2 = pb - 0.5 * A[:, None] * tg
    m2b = 0.5 * A * (d - s_) / 6 * (np.sum(q1 * q1, 1) + np.sum(q2 * q2, 1) + np.sum(q1 * q2, 1))
    V, _ = volumes(F)
    val = np.sum(m2a) + np.sum(m2b[br]) - np.sum(W * V) + np.sum(W * VG)
    if CAP is not None and Pr is not None:
        _, _, ex = cap_terms(Pr, W)
        val -= 0.5 * CAP['kappa'] * np.sum(Pr['w'] * ex**2)
    return val


def resid(F, W, VG, Pr):
    V, _ = volumes(F)
    Gv = V - VG
    if CAP is not None:
        t, act, _ = cap_terms(Pr, W)
        np.add.at(Gv, Pr['a'], t); np.add.at(Gv, Pr['jb'], -t)
    return Gv


def jac(F, W, Pr):
    H = hessian(F)
    if CAP is not None:
        t, act, _ = cap_terms(Pr, W)
        a, jb = Pr['a'][act], Pr['jb'][act]
        k = CAP['kappa'] * Pr['w'][act]
        N = F['N']
        H = H + sp.coo_matrix((np.concatenate([k, -k, k, -k]),
                               (np.concatenate([a, a, jb, jb]), np.concatenate([a, jb, jb, a]))), shape=(N, N)).tocsr()
    return H


def newton(X, VG, W, margin, tol=1e-10, itmax=60):
    F = diagram(X, W, margin)
    Pr = pairs_from(diagram(X, np.zeros_like(W), margin)) if CAP is not None else None
    ph = phi(F, W, VG, Pr)
    for it in range(itmax):
        e = resid(F, W, VG, Pr)
        en = np.max(np.abs(e)) / VG.mean()
        if en < tol:
            return W, F, it, Pr
        H = jac(F, W, Pr) + sp.identity(len(W)) * 1e-12 * H_scale(F)
        dW = spla.spsolve(H.tocsc(), -e); dW -= dW.mean()
        slope = -np.dot(e, dW)            # d(phi)/d(lam) at 0 = (VG - V - cap).dW = -e.dW > 0
        lam = 1.0
        while True:
            Wn = W + lam * dW
            Fn = diagram(X, Wn, margin); Vn, _ = volumes(Fn)
            if np.all(Vn > 0):
                phn = phi(Fn, Wn, VG, Pr)
                if phn >= ph + 1e-4 * lam * slope or \
                        np.max(np.abs(resid(Fn, Wn, VG, Pr))) / VG.mean() < 0.9 * en:
                    break
            lam *= 0.5
            if lam < 1e-8:
                if en < 1e-8:
                    return W, F, it, Pr
                raise RuntimeError(f"Newton line search failed, err={en:.3e}")
        W, F, ph = Wn, Fn, phn
    if en < 1e-8:
        return W, F, itmax, Pr
    raise RuntimeError(f"Newton no convergence, err={en:.3e}")


def newton_plain(X, VG, W, margin, tol=1e-10, itmax=40):
    for it in range(itmax):
        F = diagram(X, W, margin); V, _ = volumes(F)
        e = V - VG
        en = np.max(np.abs(e)) / VG.mean()
        if en < tol:
            return W, F, it
        H = hessian(F) + sp.identity(len(W)) * 1e-12 * H_scale(F)
        dW = spla.spsolve(H.tocsc(), -e); dW -= dW.mean()
        lam = 1.0
        while True:
            Wn = W + lam * dW
            Fn = diagram(X, Wn, margin); Vn, _ = volumes(Fn)
            if np.all(Vn > 0) and np.max(np.abs(Vn - VG)) / VG.mean() < en:
                break
            lam *= 0.5
            if lam < 1e-6:
                if en < 1e-8:                      # round-off floor: accept
                    return W, F, it
                raise RuntimeError(f"Newton line search failed, err={en:.3e}")
        W = Wn
    raise RuntimeError(f"Newton no convergence, err={en:.3e}")


def H_scale(F):
    return np.mean(F['A'] / F['d'])


def riemann(rL, pL, uL, rR, pR, uR):
    cL, cR = np.sqrt(GAM * pL / rL), np.sqrt(GAM * pR / rR)
    SL = np.minimum(uL - cL, uR - cR); SR = np.maximum(uL + cL, uR + cR)
    us = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))
    ps = pL + rL * (SL - uL) * (us - uL)
    ps = np.maximum(ps, 0.0)          # HLLC-type P* can be < 0 in strong rarefaction (Astra); clamp keeps G_f >= 0
    return ps, us


MUSCL = True


def recon(F, qs):
    """LS gradients (weights A/d) + Barth-Jespersen, values at face centroids. Returns list of (qa, qb)."""
    N, a, jb, br, A, d = F['N'], F['a'], F['jb'], F['breal'], F['A'], F['d']
    r = F['n'] * d[:, None]                      # x_b - x_a (image frame)
    w = A / d
    M = np.zeros((N, 2, 2))
    rr = w[:, None, None] * r[:, :, None] * r[:, None, :]
    np.add.at(M, a, rr); np.add.at(M, jb[br], rr[br])
    Minv = np.linalg.pinv(M)            # degenerate neighbour geometry -> least-norm gradient
    xa = F['cf'] - F['Pa']
    xb = F['cf'] - (F['Pa'] + r)
    out = []
    for q in qs:
        dq = q[jb] - q[a]
        bvec = np.zeros((N, 2))
        np.add.at(bvec, a, (w * dq)[:, None] * r)
        np.add.at(bvec, jb[br], (w * dq)[br][:, None] * r[br])
        g = np.einsum('nij,nj->ni', Minv, bvec)
        if not MUSCL:
            g[:] = 0
        qmax = q.copy(); qmin = q.copy()
        np.maximum.at(qmax, a, q[jb]); np.minimum.at(qmin, a, q[jb])
        np.maximum.at(qmax, jb[br], q[a][br]); np.minimum.at(qmin, jb[br], q[a][br])
        alpha = np.ones(N)
        for side, idx, xr, msk in (('a', a, xa, np.ones(len(a), bool)), ('b', jb, xb, br)):
            D = np.sum(g[idx] * xr, 1)
            lim = np.where(D > 1e-300, (qmax[idx] - q[idx]) / np.where(D > 1e-300, D, 1),
                           np.where(D < -1e-300, (qmin[idx] - q[idx]) / np.where(D < -1e-300, D, 1), 1.0))
            lim = np.clip(lim, 0, 1)
            np.minimum.at(alpha, idx[msk], lim[msk])
        g *= alpha[:, None]
        out.append((q[a] + np.sum(g[a] * xa, 1), q[jb] + np.sum(g[jb] * xb, 1)))
    return out


def rates(state, scheme, m, margin):
    X, v, E, VG, W = state
    X = X % LB
    if scheme == 'G':
        W, F, it, Pr = newton(X, VG, W, margin)
        Vuse = VG
    elif scheme == 'H':
        F = diagram(X, W, margin); it = 0
        Vuse, _ = volumes(F)
    else:
        W = np.zeros_like(W); F = diagram(X, W, margin); it = 0
        Vuse, _ = volumes(F)
    rho = m / Vuse
    U = E - 0.5 * m * np.sum(v * v, 1)
    P = (GAM - 1) * U / Vuse
    if PFLOOR is not None:
        P = np.maximum(P, PFLOOR)       # stage pressure floor (as the GPU path), not booked
    a, jb, br, A, n, d = F['a'], F['jb'], F['breal'], F['A'], F['n'], F['d']
    ua, ub = np.sum(v[a] * n, 1), np.sum(v[jb] * n, 1)
    (rA, rB), (pA, pB), (vxA, vxB), (vyA, vyB) = recon(F, [rho, P, v[:, 0], v[:, 1]])
    uA = vxA * n[:, 0] + vyA * n[:, 1]; uB = vxB * n[:, 0] + vyB * n[:, 1]
    ps, us = riemann(rA, pA, uA, rB, pB, uB)
    fa = -(ps * A)[:, None] * n
    Dstat = None
    if scheme == 'H':
        N = len(m)
        phi = 0.5 + (W[a] - W[jb]) / (2 * d * d)
        t = F['cf'] - (F['Pa'] + (phi * d)[:, None] * n)
        dP = P[jb] - P[a]
        fa = -(ps * A)[:, None] * n - (A * dP / d)[:, None] * t
        ut = phi * ua + (1 - phi) * ub
        rot = np.sum(t * (v[a] - v[jb]), 1) / d
        kin = np.zeros(N); np.add.at(kin, a, A * (ut + rot)); np.add.at(kin, jb[br], -A[br] * (ut + rot)[br])
        god = np.zeros(N); np.add.at(god, a, A * us); np.add.at(god, jb[br], -A[br] * us[br])
        Vvor, _ = volumes(diagram(X, np.zeros_like(W), margin))
        VT = VG
        dVT = god + (Vvor - VT) / H_TAUV
        rhs = dVT - kin + (VT - Vuse) / H_TAUC
        Hm = hessian(F) + sp.identity(N) * 1e-12 * H_scale(F)
        # barrier inside the Wdot law: pairs with |dW|/d^2 > H_R0 are driven back to the threshold over H_TAUB
        dWf = W[a] - W[jb]
        exc = np.abs(dWf) - H_R0 * d * d
        act = exc > 0
        nbar = int(act.sum())
        if nbar:
            k = H_KB * A[act] / d[act]
            ddot_b = ub[act] - ua[act]
            tgt = (dWf[act] / d[act]**2) * 2 * d[act] * ddot_b - np.sign(dWf[act]) * exc[act] / H_TAUB   # keep dW/d^2, relax excess
            aa, bb, brr = a[act], jb[act], br[act]
            Hm = Hm + sp.coo_matrix((np.concatenate([k, -k, k[brr], -k[brr]]),
                                     (np.concatenate([aa, aa, bb[brr], bb[brr]]),
                                      np.concatenate([aa, bb, bb[brr], aa[brr]]))), shape=(N, N)).tocsr()
            np.add.at(rhs, aa, k * tgt); np.add.at(rhs, bb[brr], -(k * tgt)[brr])
        Wd = spla.spsolve(Hm.tocsc(), rhs); Wd -= Wd.mean()
        dl = (Wd[a] - Wd[jb]) / (2 * d)
        wf = ut + rot + dl
        Pt = phi * P[jb] + (1 - phi) * P[a]
        D = A * ((P[a] - P[jb]) * dl + (ps - Pt) * (ua - ub))
        ga = A * np.maximum(0.0, -(ps - pA) * (us - uA)); gb = A * np.maximum(0.0, (ps - pB) * (us - uB))
        gs = ga + gb
        tha = np.where(gs > 0, ga / np.where(gs > 0, gs, 1), 0.5)
        dEa = np.sum(v[a] * fa, 1) - P[a] * A * wf + tha * D
        dEb = -np.sum(v[jb] * fa, 1) + P[jb] * A * wf + (1 - tha) * D
        dv = np.zeros_like(v); dE = np.zeros(N)
        np.add.at(dv, a, fa); np.add.at(dE, a, dEa)
        np.add.at(dv, jb[br], -fa[br]); np.add.at(dE, jb[br], dEb[br])
        Dw = np.where(br, 1.0, 0.5)
        Dstat = (np.mean(D < 0), np.sum(Dw * np.minimum(D, 0)), np.sum(Dw * np.maximum(D, 0)))
        if GRAV is not None:
            acc, _ = GRAV(X)
            dE += m * np.sum(v * acc, 1)
            dv += m[:, None] * acc
        dv /= m[:, None]
        # validity time-step limit: keep |W_a - W_b| < (1 - H_DV) d^2 on every face over the step
        dWd = Wd[a] - Wd[jb]; ddot = ub - ua
        lim = (1 - H_DV) * d * d
        mg = lim - np.abs(dWf)
        gr = np.where(dWf != 0, np.sign(dWf) * dWd, np.abs(dWd)) - (1 - H_DV) * 2 * d * ddot
        ok = (mg > 0) & (gr > 0)
        dtlim = 0.5 * np.min(mg[ok] / gr[ok]) if ok.any() else 1e30
        info = dict(W=W, F=F, it=nbar, P=P, rho=rho, V=Vuse, U=U, Pr=None, Dstat=Dstat, dtlim=dtlim)
        return (v, dv, dE, dVT, Wd), info
    if scheme == 'F':
        # (F), W = 0: VPH force + Riemann correction; dU = -P dV_geo + sigma (Docs §3)
        t = F['cf'] - (F['Pa'] + (0.5 * d)[:, None] * n)
        dP = P[jb] - P[a]
        fa = fa - (A * dP / d)[:, None] * t
        wf = 0.5 * (ua + ub) + np.sum(t * (v[a] - v[jb]), 1) / d
        D = A * (ps - 0.5 * (P[a] + P[jb])) * (ua - ub)
        ga = A * np.maximum(0.0, -(ps - pA) * (us - uA)); gb = A * np.maximum(0.0, (ps - pB) * (us - uB))
        gs = ga + gb
        tha = np.where(gs > 0, ga / np.where(gs > 0, gs, 1), 0.5)
        dEa = np.sum(v[a] * fa, 1) - P[a] * A * wf + tha * D
        dEb = -np.sum(v[jb] * fa, 1) + P[jb] * A * wf + (1 - tha) * D
        dv = np.zeros_like(v); dE = np.zeros(len(m)); dVG = np.zeros(len(m))
        np.add.at(dv, a, fa); np.add.at(dE, a, dEa)
        np.add.at(dv, jb[br], -fa[br]); np.add.at(dE, jb[br], dEb[br])
        Dw = np.where(br, 1.0, 0.5)
        Dstat = (np.mean(D < 0), np.sum(Dw * np.minimum(D, 0)), np.sum(Dw * np.maximum(D, 0)))
        if GRAV is not None:
            acc, _ = GRAV(X)
            dE += m * np.sum(v * acc, 1)
            dv += m[:, None] * acc
        dv /= m[:, None]
        info = dict(W=W, F=F, it=it, P=P, rho=rho, V=Vuse, U=U, Pr=None, Dstat=Dstat)
        return (v, dv, dE, dVG), info
    if scheme == 'G':
        wf = us
    else:                                   # geometric face-centroid normal speed (W=0)
        phi = 0.5
        t = F['cf'] - (F['Pa'] + (phi * d)[:, None] * n)
        wf = phi * ua + (1 - phi) * ub + np.sum(t * (v[a] - v[jb]), 1) / d
    dE_face = -ps * A * wf
    dv = np.zeros_like(v); dE = np.zeros(len(m)); dVG = np.zeros(len(m))
    np.add.at(dv, a, fa); np.add.at(dE, a, dE_face); np.add.at(dVG, a, A * us)
    np.add.at(dv, jb[br], -fa[br]); np.add.at(dE, jb[br], -dE_face[br]); np.add.at(dVG, jb[br], -A[br] * us[br])
    if GRAV is not None:
        acc, _ = GRAV(X)
        dE += m * np.sum(v * acc, 1)
        dv += m[:, None] * acc
    dv /= m[:, None]
    info = dict(W=W, F=F, it=it, P=P, rho=rho, V=Vuse, U=U, Pr=Pr if scheme == 'G' else None, Dstat=None)
    info['fd'] = dict(a=a, jb=jb, cf=F['cf'], dun=uB - uA, pA=pA, pB=pB, rA=rA, rB=rB, ps=ps)
    return (v, dv, dE, dVG), info


def rk4_step(state, dt, scheme, m, margin):
    k1, i1 = rates(state, scheme, m, margin)
    nI = len(k1)

    def adv(k, c, info):
        s_ = [state[j] + c * k[j] for j in range(4)]
        s_.append(state[4] + c * k[4] if nI == 5 else info['W'])
        return s_
    k2, i2 = rates(adv(k1, dt / 2, i1), scheme, m, margin)
    k3, i3 = rates(adv(k2, dt / 2, i2), scheme, m, margin)
    k4, i4 = rates(adv(k3, dt, i3), scheme, m, margin)
    new = [state[j] + dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(4)]
    new.append(state[4] + dt / 6 * (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4]) if nI == 5 else i4['W'])
    return new


def gresho(N1, jitter=1e-3):
    g = (np.arange(N1) + 0.5) / N1
    X = np.array([(x, y) for y in g for x in g])
    X += np.random.default_rng(1).uniform(-jitter, jitter, X.shape) / N1
    R = X - 0.5; r = np.linalg.norm(R, axis=1)
    vphi = np.where(r < 0.2, 5 * r, np.where(r < 0.4, 2 - 5 * r, 0.0))
    P = np.where(r < 0.2, 5 + 12.5 * r**2,
                 np.where(r < 0.4, 9 + 12.5 * r**2 - 20 * r + 4 * np.log(np.maximum(r, 1e-30) / 0.2), 3 + 4 * np.log(2)))
    eph = np.stack([-R[:, 1], R[:, 0]], 1) / np.maximum(r, 1e-30)[:, None]
    v = vphi[:, None] * eph
    return X, v, P


def vphi_exact(r):
    return np.where(r < 0.2, 5 * r, np.where(r < 0.4, 2 - 5 * r, 0.0))


def diag(state, m, info, t):
    X, v, E, VG, W = state
    R = (X % LB) - 0.5; r = np.linalg.norm(R, axis=1)
    eph = np.stack([-R[:, 1], R[:, 0]], 1) / np.maximum(r, 1e-30)[:, None]
    vp = np.sum(v * eph, 1)
    V = info['V']
    L1 = np.sum(V * np.abs(vp - vphi_exact(r))) / np.sum(V)
    Lz = np.sum(m * (R[:, 0] * v[:, 1] - R[:, 1] * v[:, 0]))
    F = info['F']
    _, s = volumes(F)
    frac = np.abs(F['W'][F['a']] - F['W'][F['jb']]) / F['d']**2
    frac = np.array([frac.max(), np.percentile(frac, 99)])
    out = np.zeros(len(m), bool)
    out[F['a'][s < 0]] = True
    br = F['breal']; out[F['jb'][br & (s > F['d'])]] = True
    return L1, np.sum(E), Lz, frac, out.sum()


CAPLOG = dict(vol=0.0, sneg=0.0, nact=0)


def cap_transfer(state, m, margin, info):
    if CAP is None:
        return state
    Pr = info['Pr']; W = info['W']
    t, act, _ = cap_terms(Pr, W)
    if not act.any():
        return state
    a, jb = Pr['a'], Pr['jb']
    P = info['P']
    Pf = 0.5 * (P[a] + P[jb])
    dV = np.zeros(len(m)); dU = np.zeros(len(m))
    np.add.at(dV, a, -t); np.add.at(dV, jb, t)
    np.add.at(dU, a, Pf * t); np.add.at(dU, jb, -Pf * t)
    X, v, E, VG, Wst = state
    CAPLOG['vol'] += np.sum(np.abs(t)); CAPLOG['nact'] += int(act.sum())
    CAPLOG['sneg'] += np.sum(np.minimum(0.0, dU + P * dV))
    return [X, v, E + dU, VG + dV, Wst]


def run(scheme, N1, tend, C=0.3, report=None):
    X, v, P = gresho(N1)
    margin = 4.0 / N1
    W = np.zeros(len(X))
    F = diagram(X, W, margin); V, _ = volumes(F)
    m = 1.0 * V                              # rho = 1
    E = 0.5 * m * np.sum(v * v, 1) + P * V / (GAM - 1)
    state = [X, v, E, V.copy(), W]
    t, step = 0.0, 0
    _, info = rates(state, scheme, m, margin)
    L1, E0, Lz0, fr, no = diag(state, m, info, t)
    print(f"[{scheme} N={N1}] t=0 L1={L1:.4e} E={E0:.10f} Lz={Lz0:.6e}", flush=True)
    t0 = time.time()
    report = report or tend / 10
    tnext = report
    while t < tend - 1e-12:
        c = np.sqrt(GAM * info['P'] / info['rho'])
        dt = min(C * np.min(np.sqrt(info['V']) / (c + np.linalg.norm(state[1], axis=1))), tend - t)
        state = rk4_step(state, dt, scheme, m, margin)
        state[0] = state[0] % LB
        t += dt; step += 1
        _, info = rates(state, scheme, m, margin)
        if CAP is not None:
            state = cap_transfer(state, m, margin, info); _, info = rates(state, scheme, m, margin)
        if t >= tnext - 1e-12 or t >= tend - 1e-12:
            L1, Et, Lz, fr, no = diag(state, m, info, t)
            print(f"[{scheme} N={N1}] t={t:.3f} step={step} dt={dt:.2e} L1(vphi)={L1:.4e} dE/E0={(Et-E0)/E0:+.2e} "
                  f"dLz/Lz0={(Lz-Lz0)/Lz0:+.2e} max|dW|/d2={fr[0]:.3f} p99={fr[1]:.3f} n_out={no} newton_it={info['it']} "
                  f"minP={info['P'].min():.3f} {'' if info['Dstat'] is None else 'Dneg_frac=%.3f Dneg/Dpos=%.2e ' % (info['Dstat'][0], -info['Dstat'][1]/max(info['Dstat'][2],1e-300))}cap_vol/V={CAPLOG['vol']/np.sum(info['V']):.2e} cap_sneg/E={CAPLOG['sneg']/E0:.2e} nact={CAPLOG['nact']} wall={time.time()-t0:.0f}s", flush=True)
            tnext += report
    np.savez(f"final_{scheme}_N{N1}.npz", X=state[0], v=state[1], E=state[2], VG=state[3], W=state[4], m=m, t=t)
    return state


# ---------------- Kepler disk (production "old 1/R" IC) ----------------
KEPS = 0.047
KC = np.array([12.0, 12.0])


def kgrav(X):
    R = X - KC
    q = np.sum(R * R, 1) + KEPS**2
    return -R / q[:, None]**1.5, -1.0 / np.sqrt(q)


def kepler_ic(N1):
    g = (np.arange(N1) + 0.5) / N1 * LB
    X = np.array([(x, y) for y in g for x in g])
    X += np.random.default_rng(1).uniform(-1e-6, 1e-6, X.shape) * LB / N1   # regular grid as the code; jitter only for Qhull
    R = X - KC; r = np.linalg.norm(R, axis=1)
    inner = 0.5 * (1 + np.tanh((r - 2.0) / 0.25)); outer = 0.5 * (1 + np.tanh((8.0 - r) / 0.25))
    rho = 1e-3 + (2.0 / (r + 1e-6)) * inner * outer
    vamp = r * (r * r + KEPS**2) ** -0.75
    v = np.stack([-R[:, 1], R[:, 0]], 1) / np.maximum(r, 1e-12)[:, None] * vamp[:, None]
    return X, v, rho, np.full(len(r), 1e-6)


RBINS = np.array([1.0, 2.0, 4.0, 8.0])        # bins: <1, 1-2, 2-4, 4-8, >8
FBIN = np.zeros(5); NBIN = np.zeros(5, int)
IN0 = None


def kdiag(state, info, m, t, aE, in0):
    X, v = state[0] % LB, state[1]
    r = X - KC; R = np.linalg.norm(r, axis=1)
    lz = m * (r[:, 0] * v[:, 1] - r[:, 1] * v[:, 0])
    ib = np.digitize(R, RBINS)
    Lb = np.bincount(ib, weights=lz, minlength=5); Nb = np.bincount(ib, minlength=5)
    fd = info['fd']; a, jb = fd['a'], fd['jb']
    Rf = np.linalg.norm(fd['cf'] - KC, axis=1); fb = np.digitize(Rf, RBINS)
    cs = np.sqrt(GAM * np.maximum(np.minimum(fd['pA'], fd['pB']), 1e-300) / np.maximum(np.maximum(fd['rA'], fd['rB']), 1e-300))
    mach = np.abs(fd['dun']) / cs
    comp = fd['dun'] < 0
    pr = fd['ps'] / np.maximum(0.5 * (fd['pA'] + fd['pB']), 1e-300)
    F = info['F']; dmin = F['d'].min()
    print(f"[DIAG] t={t:.3f} minR={R.min():.4f} dmin={dmin:.4f} inner24: meanR={R[in0].mean():.3f} Lz={lz[in0].sum():.4e} "
          f"floor_by_R(/|E0|)={' '.join('%.2e' % (x / aE) for x in FBIN)} nfloor_by_R={' '.join(str(x) for x in NBIN)}", flush=True)
    print(f"[DIAG] t={t:.3f} Lz_by_R={' '.join('%.4e' % x for x in Lb)} N_by_R={' '.join(str(x) for x in Nb)}", flush=True)
    for b, name in enumerate(('R<1', '1-2', '2-4', '4-8', '>8')):
        sel = fb == b
        if not sel.any():
            continue
        mc = mach[sel & comp]; prc = pr[sel & comp]
        print(f"[FACE] t={t:.3f} {name:4s} nf={sel.sum():5d} comp={int((sel & comp).sum()):5d} "
              f"Mach_comp[med,p90,max]=[{np.median(mc) if mc.size else 0:.2f},{np.percentile(mc,90) if mc.size else 0:.2f},{mc.max() if mc.size else 0:.2f}] "
              f"frac_supersonic_comp={np.mean(mc > 1) if mc.size else 0:.3f} Pstar/P[med,p90,max]=[{np.median(prc) if prc.size else 0:.2f},{np.percentile(prc,90) if prc.size else 0:.2f},{prc.max() if prc.size else 0:.2f}]", flush=True)


def run_kepler(scheme, N1, tend, report, C=0.4, eta=0.06):
    global LB, GRAV, PFLOOR
    LB = 24.0; GRAV = kgrav; PFLOOR = 1e-6
    X, v, rho0, P0 = kepler_ic(N1)
    margin = 4.0 * LB / N1
    W = np.zeros(len(X))
    F = diagram(X, W, margin); V, _ = volumes(F)
    m = rho0 * V
    E = 0.5 * m * np.sum(v * v, 1) + P0 * V / (GAM - 1)
    state = [X, v, E, V.copy(), W]
    Etot0 = np.sum(E) + np.sum(m * kgrav(X)[1]); aE = abs(Etot0)
    global IN0
    IN0 = np.linalg.norm(X - KC, axis=1) < 1.0
    t, step, floor_cum, nfl = 0.0, 0, 0.0, 0
    import os
    if os.environ.get('KEP_RESTART'):
        z = np.load(os.environ['KEP_RESTART'])
        state = [z['X'], z['v'], z['E'], z['VG'], z['W']]; m = z['m']
        t, step, floor_cum, nfl = float(z['t']), int(z['step']), float(z['floor_cum']), int(z['nfl'])
        Etot0 = float(z['Etot0']); aE = abs(Etot0)
        print(f"[K{scheme} N={N1}] restart from {os.environ['KEP_RESTART']} t={t}", flush=True)
    tck = np.floor(t) + 1.0
    _, info = rates(state, scheme, m, margin)
    print(f"[K{scheme} N={N1}] t=0 Etot0={Etot0:.6f} Epot={np.sum(m*kgrav(X)[1]):.6f}", flush=True)
    t0 = time.time(); tnext = (np.floor(t / report) + 1) * report
    while t < tend - 1e-12:
        c = np.sqrt(GAM * info['P'] / info['rho'])
        Rq = np.sum((state[0] - KC)**2, 1) + KEPS**2
        dt = min(C * np.min(np.sqrt(info['V']) / (c + np.linalg.norm(state[1], axis=1))),
                 eta * np.min(Rq**0.75), info.get('dtlim', 1e30), tend - t)
        prev = [x.copy() for x in state]
        try:
            state = rk4_step(state, dt, scheme, m, margin)
            if not all(np.all(np.isfinite(x)) for x in state):
                raise FloatingPointError("non-finite state")
        except Exception as ex:
            np.savez(f"crash_prev_{scheme}_N{N1}.npz", X=prev[0], v=prev[1], E=prev[2], VG=prev[3], W=prev[4], m=m,
                     t=t, step=step, dt=dt, floor_cum=floor_cum, nfl=nfl, Etot0=Etot0)
            print(f"[K{scheme} N={N1}] CRASH at t={t:.5f} step={step} dt={dt:.3e}: {ex!r}; saved previous state", flush=True)
            raise
        state[0] = state[0] % LB
        t += dt; step += 1
        _, info = rates(state, scheme, m, margin)
        if t >= tck - 1e-12:
            np.savez(f"ckpt_{scheme}_N{N1}.npz", X=state[0], v=state[1], E=state[2], VG=state[3], W=state[4], m=m,
                     t=t, step=step, floor_cum=floor_cum, nfl=nfl, Etot0=Etot0)
            tck += 1.0
        if CAP is not None:
            state = cap_transfer(state, m, margin, info); _, info = rates(state, scheme, m, margin)
        # end-of-step floor as exam.c RK4 path: ie <= 0 -> P = 1e-6, booked in floor_cum
        U = state[2] - 0.5 * m * np.sum(state[1]**2, 1)
        low = U <= 0
        if low.any():
            Ufl = 1e-6 * info['V'] / (GAM - 1)
            floor_cum += np.sum(Ufl[low] - U[low]); nfl += int(low.sum())
            Rl = np.linalg.norm(state[0][low] - KC, axis=1)
            ib = np.digitize(Rl, RBINS)
            np.add.at(FBIN, ib, (Ufl[low] - U[low])); np.add.at(NBIN, ib, 1)
            state[2] = np.where(low, state[2] + (Ufl - U), state[2])
            _, info = rates(state, scheme, m, margin)
        if t >= tnext - 1e-12 or t >= tend - 1e-12:
            Et = np.sum(state[2]) + np.sum(m * kgrav(state[0])[1])
            Ud = state[2] - 0.5 * m * np.sum(state[1]**2, 1)
            Fd = info['F']; _, s_ = volumes(Fd)
            fr = np.abs(Fd['W'][Fd['a']] - Fd['W'][Fd['jb']]) / Fd['d']**2
            nout = np.zeros(len(m), bool); nout[Fd['a'][s_ < 0]] = True
            br = Fd['breal']; nout[Fd['jb'][br & (s_ > Fd['d'])]] = True
            Rr = np.linalg.norm(state[0] - KC, axis=1); dsk = (Rr > 2.5) & (Rr < 7.5)
            print(f"[K{scheme} N={N1}] t={t:.3f} step={step} dt={dt:.2e} dE/|E0|={(Et-Etot0)/aE:+.3e} "
                  f"floor_cum/|E0|={floor_cum/aE:.3e} nfloor={nfl} Eint={np.sum(Ud):.3e} "
                  f"max|dW|/d2={fr.max():.3f} p99={np.percentile(fr,99):.3f} n_out={nout.sum()} "
                  f"newton_it={info['it']} {'' if info['Dstat'] is None else 'Dneg_frac=%.3f Dneg/Dpos=%.2e ' % (info['Dstat'][0], -info['Dstat'][1]/max(info['Dstat'][2],1e-300))}cap_vol/V={CAPLOG['vol']/np.sum(info['V']):.2e} cap_sneg/|E0|={CAPLOG['sneg']/aE:.2e} nact={CAPLOG['nact']} diskP[min,med,max]=[{info['P'][dsk].min():.2e},{np.median(info['P'][dsk]):.2e},{info['P'][dsk].max():.2e}] "
                  f"wall={time.time()-t0:.0f}s", flush=True)
            kdiag(state, info, m, t, aE, IN0)
            tnext += report
    np.savez(f"kep_{scheme}_N{N1}.npz", X=state[0], v=state[1], E=state[2], VG=state[3], W=state[4], m=m, t=t)


if __name__ == '__main__':
    import os
    if os.environ.get('H_TAUV'):
        H_TAUV = float(os.environ['H_TAUV']); H_TAUC = float(os.environ.get('H_TAUC', '0.02'))
        H_R0 = float(os.environ.get('H_R0', '0.5')); H_KB = float(os.environ.get('H_KB', '10')); H_TAUB = float(os.environ.get('H_TAUB', '0.02'))
        print(f"H_R0={H_R0} H_KB={H_KB} H_TAUB={H_TAUB}", flush=True)
        print(f"H_TAUV={H_TAUV} H_TAUC={H_TAUC}", flush=True)
    if os.environ.get('CAP_DELTA'):
        CAP = dict(delta=float(os.environ['CAP_DELTA']), kappa=float(os.environ.get('CAP_KAPPA', '1000')))
        print(f"CAP={CAP}", flush=True)
    if sys.argv[1] == 'kepler':
        MUSCL = sys.argv[5] != '0' if len(sys.argv) > 5 else True
        run_kepler(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[6]) if len(sys.argv) > 6 else 0.5)
        sys.exit(0)
    scheme = sys.argv[1]; N1 = int(sys.argv[2]); tend = float(sys.argv[3])
    if len(sys.argv) > 5: MUSCL = sys.argv[5] != '0'
    print(f"MUSCL={MUSCL}", flush=True)
    run(scheme, N1, tend, report=float(sys.argv[4]) if len(sys.argv) > 4 else None)
