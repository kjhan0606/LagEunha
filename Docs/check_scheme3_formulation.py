"""Finite-difference checks of Docs/gfs_scheme3_formulation.md (2D periodic Laguerre cells)."""
import numpy as np
from collections import defaultdict

rng = np.random.default_rng(20260929)
L = 1.0
n1 = 6
N = n1 * n1
g = (np.arange(n1) + 0.5) / n1
X0 = np.array([(a, b) for a in g for b in g]) + rng.uniform(-0.25, 0.25, (N, 2)) / n1
X0 %= L
dmean = L / n1
W0 = rng.uniform(0, 0.25, N) * dmean**2
shifts = np.array([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)], float) * L


def clip(poly, a, b, lab):
    out = []
    m = len(poly)
    for k in range(m):
        P, lp = poly[k]
        Q, _ = poly[(k + 1) % m]
        fP, fQ = a @ P - b, a @ Q - b
        inP, inQ = fP <= 0, fQ <= 0
        if inP:
            out.append((P, lp))
        if inP != inQ:
            I = P + (Q - P) * (fP / (fP - fQ))
            out.append((I, lab if inP else lp))
    return out


def cells(X, W):
    """Return volumes V[i] and faces dict (i, j, s) -> (A, c) with s the image shift index."""
    V = np.zeros(N)
    faces = {}
    for i in range(N):
        xi = X[i]
        h = 1.5 * L
        poly = [(xi + np.array(p) * h, -1) for p in [(-1, -1), (1, -1), (1, 1), (-1, 1)]]
        for s, sh in enumerate(shifts):
            for j in range(N):
                if j == i and s == 4:
                    continue
                xj = X[j] + sh
                if np.sum((xj - xi) ** 2) > (2.5 * dmean) ** 2:
                    continue
                a = 2 * (xj - xi)
                b = xj @ xj - xi @ xi + W[i] - W[j]
                poly = clip(poly, a, b, (j, s))
        P = np.array([p for p, _ in poly])
        x, y = P[:, 0], P[:, 1]
        V[i] = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
        m = len(poly)
        for k in range(m):
            lab = poly[k][1]
            if lab == -1:
                raise RuntimeError("cell touches initial box")
            A = np.linalg.norm(poly[(k + 1) % m][0] - poly[k][0])
            if A > 1e-14:
                faces[(i, lab[0], lab[1])] = (A, 0.5 * (poly[k][0] + poly[(k + 1) % m][0]))
    return V, faces


V0, F0 = cells(X0, W0)
print(f"N={N}, sum V = {V0.sum():.15f} (box {L*L}), faces = {len(F0)}")
h = 1e-6


def rel(a, b):
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)


# 1. dV/dW = 1/2 L_{A/d}
J_fd = np.zeros((N, N))
for k in range(N):
    e = np.zeros(N); e[k] = h
    J_fd[:, k] = (cells(X0, W0 + e)[0] - cells(X0, W0 - e)[0]) / (2 * h)
J_th = np.zeros((N, N))
for (i, j, s), (A, c) in F0.items():
    d = np.linalg.norm(X0[j] + shifts[s] - X0[i])
    J_th[i, i] += A / (2 * d)
    J_th[i, j] -= A / (2 * d)
print(f"[1] dV/dW = L_(A/d)/2        max rel err = {rel(J_fd, J_th):.2e}")

# 2. dV/dx
Jx_fd = np.zeros((N, N, 2))
for k in range(N):
    for a in range(2):
        e = np.zeros((N, 2)); e[k, a] = h
        Jx_fd[:, k, a] = (cells(X0 + e, W0)[0] - cells(X0 - e, W0)[0]) / (2 * h)
Jx_th = np.zeros((N, N, 2))
for (i, j, s), (A, c) in F0.items():
    xj = X0[j] + shifts[s]
    d = np.linalg.norm(xj - X0[i])
    Jx_th[i, i] += A * (c - X0[i]) / d
    Jx_th[i, j] -= A * (c - xj) / d
print(f"[2] dV/dx                     max rel err = {rel(Jx_fd, Jx_th):.2e}")

# 3. face velocity decomposition w = u~ + rho + delta, dV/dt = sum A w
v = rng.normal(0, 1, (N, 2))
Wd = rng.normal(0, 1, N) * dmean**2
dV_fd = (cells(X0 + h * v, W0 + h * Wd)[0] - cells(X0 - h * v, W0 - h * Wd)[0]) / (2 * h)
dV_th = np.zeros(N)
geo = {}
for (i, j, s), (A, c) in F0.items():
    xj = X0[j] + shifts[s]
    d = np.linalg.norm(xj - X0[i]); n = (xj - X0[i]) / d
    phi = 0.5 + (W0[i] - W0[j]) / (2 * d * d)
    t = c - (X0[i] + phi * d * n)
    ui, uj = v[i] @ n, v[j] @ n
    ut = phi * ui + (1 - phi) * uj
    rho = t @ (v[i] - v[j]) / d
    dl = (Wd[i] - Wd[j]) / (2 * d)
    dV_th[i] += A * (ut + rho + dl)
    geo[(i, j, s)] = (A, d, n, phi, t, ui, uj, ut, dl)
print(f"[3] dV/dt = sum A (u~+rho+dl) max rel err = {rel(dV_fd, dV_th):.2e}")

# 4. force (F) pairwise antisymmetric; energy identity sum(v.F + dU) = 0 with dU = -P dV + sigma
P = np.exp(rng.normal(0, 1.5, N))
Pstar = {}
for key in geo:
    i, j, s = key
    kk = (j, i, (8 - s))
    Pstar[key] = Pstar.get(kk, np.exp(rng.normal(0, 1.5)))
Fvec = np.zeros((N, 2)); asym = 0.0; Dsum = 0.0
for key, (A, d, n, phi, t, ui, uj, ut, dl) in geo.items():
    i, j, s = key
    f = -A * (Pstar[key] * n + (P[j] - P[i]) * t / d)
    Fvec[i] += f
    kk = (j, i, 8 - s)
    A2, d2, n2, phi2, t2, *_ = geo[kk]
    f2 = -A2 * (Pstar[kk] * n2 + (P[i] - P[j]) * t2 / d2)
    asym = max(asym, np.linalg.norm(f + f2) / np.linalg.norm(f))
    Pt = phi * P[j] + (1 - phi) * P[i]
    Dsum += 0.5 * A * ((P[i] - P[j]) * dl + (Pstar[key] - Pt) * (ui - uj))  # each face seen twice
print(f"[4a] f_ij + f_ji                max rel     = {asym:.2e}")
print(f"     sum F                      = {np.abs(Fvec.sum(0)).max():.2e}")
lhs = np.sum(v * Fvec) - np.sum(P * dV_fd) + Dsum
print(f"[4b] sum(v.F) - sum(P dV_fd) + sum D_f = {lhs:.2e}  (scale {np.sum(np.abs(P*dV_fd)):.2e})")
# VPH part equals sum_k P_k dV_k/dx_i
F_vph_fd = np.einsum('k,kia->ia', P, Jx_fd)
F_vph = np.zeros((N, 2))
for key, (A, d, n, phi, t, *_ ) in geo.items():
    i, j, s = key
    Pt = phi * P[j] + (1 - phi) * P[i]
    F_vph[i] += -A * (Pt * n + (P[j] - P[i]) * t / d)
print(f"[4c] VPH pair force vs sum P dV/dx (FD) max rel err = {rel(F_vph, F_vph_fd):.2e}")

# 5. D_f(eps=0) = G_f, any phi
m = 100000
PL, PR = np.exp(rng.normal(0, 2, m)), np.exp(rng.normal(0, 2, m))
uL, uR, us, ps = rng.normal(0, 1, (4, m))
phi = rng.uniform(0.1, 0.9, m)
ut = phi * uL + (1 - phi) * uR; Pt = phi * PR + (1 - phi) * PL
D = (PL - PR) * (us - ut) + (ps - Pt) * (uL - uR)
G = -(ps - PL) * (us - uL) + (ps - PR) * (us - uR)
print(f"[5] D_f(eps=0) - G_f           max rel     = {np.max(np.abs(D-G)/(np.abs(G)+np.abs(D)+1e-300)):.2e}")


# 6. Lax monotonicity: code's hll_star_state vs HLLC pressure
def hll_code(rL, pL, uL, cL, rR, pR, uR, cR, gam):
    SL = np.minimum(uL - cL, uR - cR); SR = np.maximum(uL + cL, uR + cR)
    gm1 = gam - 1
    EL = rL * (pL / (gm1 * rL) + 0.5 * uL**2); ER = rR * (pR / (gm1 * rR) + 0.5 * uR**2)
    den = SR - SL
    rS = (SR * rR - SL * rL - (rR * uR - rL * uL)) / den
    mS = (SR * rR * uR - SL * rL * uL - ((rR * uR**2 + pR) - (rL * uL**2 + pL))) / den
    ES = (SR * ER - SL * EL - (uR * (ER + pR) - uL * (EL + pL))) / den
    vS = mS / rS
    PS = gm1 * (ES - 0.5 * rS * vS**2)
    bad = ~(PS > 0)
    PS = np.where(bad, 0.5 * (pL + pR), PS); vS = np.where(bad, 0.5 * (uL + uR), vS)
    return PS, vS


def hllc(rL, pL, uL, cL, rR, pR, uR, cR):
    SL = np.minimum(uL - cL, uR - cR); SR = np.maximum(uL + cL, uR + cR)
    S = (pR - pL + rL * uL * (SL - uL) - rR * uR * (SR - uR)) / (rL * (SL - uL) - rR * (SR - uR))
    return pL + rL * (SL - uL) * (S - uL), S


gam = 5 / 3
for name, sig_u in (("random", 1.0), ("near-rest |du|<<c", 1e-3)):
    rL, rR = np.exp(rng.normal(0, 2, m)), np.exp(rng.normal(0, 2, m))
    pL, pR = np.exp(rng.normal(0, 3, m)), np.exp(rng.normal(0, 3, m))
    cL, cR = np.sqrt(gam * pL / rL), np.sqrt(gam * pR / rR)
    uL = rng.normal(0, sig_u, m) * np.maximum(cL, cR); uR = rng.normal(0, sig_u, m) * np.maximum(cL, cR)
    for sol in ("hll_code", "hllc"):
        if sol == "hll_code":
            ps, us = hll_code(rL, pL, uL, cL, rR, pR, uR, cR, gam)
        else:
            ps, us = hllc(rL, pL, uL, cL, rR, pR, uR, cR)
        vl = (ps - pL) * (us - uL) > 0
        vr = (ps - pR) * (us - uR) < 0
        G = -(ps - pL) * (us - uL) + (ps - pR) * (us - uR)
        print(f"[6] {name:18s} {sol:8s}: (L) left viol {vl.mean():.3f}, right viol {vr.mean():.3f}, G_f<0 frac {(G<0).mean():.3f}")
