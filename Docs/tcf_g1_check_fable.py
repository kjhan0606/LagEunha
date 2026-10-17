import numpy as np
from scipy.spatial import Voronoi

rng = np.random.default_rng(3)
N = 40
L = 1.0
X0 = rng.random((N, 2)) * L

def tile(X):
    pts = []
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            pts.append(X + np.array([a, b]) * L)
    return np.vstack(pts)  # central copy = index 4*N .. 5*N-1

def poly_area(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)

def volumes(X):
    vor = Voronoi(tile(X))
    V = np.zeros(N)
    for k in range(N):
        reg = vor.regions[vor.point_region[4 * N + k]]
        P = vor.vertices[reg]
        V[k] = abs(poly_area(P))
    return V, vor

def faces(X):
    """per face of central cells: (i, j(mod N), A, e, t, d) with t = c - m"""
    vor = Voronoi(tile(X))
    out = []
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        if -1 in rv:
            continue
        for a, b in ((p, q), (q, p)):
            if 4 * N <= a < 5 * N:
                i = a - 4 * N
                xi, xj = vor.points[a], vor.points[b]
                v1, v2 = vor.vertices[rv[0]], vor.vertices[rv[1]]
                A = np.linalg.norm(v2 - v1)
                c = 0.5 * (v1 + v2)
                d = np.linalg.norm(xj - xi)
                e = (xj - xi) / d
                m = 0.5 * (xi + xj)
                t = c - m
                out.append((i, b % N, A, e, t, d))
    return out

# 1. finite-difference dV_k/dx_i vs formula
V0, _ = volumes(X0)
F = faces(X0)
eps = 1e-6
# analytic: dV_i/dx_j = A(e/2 - t/d); dV_i/dx_i = sum_f A(e/2 + t/d)
dVdx_formula = {}
for (i, j, A, e, t, d) in F:
    dVdx_formula[(i, j)] = dVdx_formula.get((i, j), 0) + A * (e / 2 - t / d)
    dVdx_formula[(i, i)] = dVdx_formula.get((i, i), 0) + A * (e / 2 + t / d)
err_formula = 0; err_docsign = 0; nchk = 0
for (i, j), g in list(dVdx_formula.items())[:200]:
    fd = np.zeros(2)
    for a in range(2):
        Xp = X0.copy(); Xp[j, a] += eps
        Xm = X0.copy(); Xm[j, a] -= eps
        fd[a] = (volumes(Xp)[0][i] - volumes(Xm)[0][i]) / (2 * eps)
    err_formula = max(err_formula, np.abs(fd - g).max())
    # doc sign: dV_i/dx_j = A(e/2 + t/d)
    if i != j:
        gd = sum(A * (e / 2 + t / d) for (ii, jj, A, e, t, d) in F if ii == i and jj == j)
        err_docsign = max(err_docsign, np.abs(fd - gd).max())
    nchk += 1
print("checked", nchk, "derivatives; max|FD - A(e/2 - t/d)| =", err_formula,
      "; max|FD - doc A(e/2 + t/d)| =", err_docsign)

# 2. energy identity: sum v.F - sum P dV/dt = 0, with F_i(f) = -A[Pbar e + (Pj - Pi) t/d]
P = rng.random(N) + 0.5
v = rng.normal(size=(N, 2))
Fi = np.zeros((N, 2)); Fi_doc = np.zeros((N, 2))
per_face = []
for (i, j, A, e, t, d) in F:
    Fi[i] += -A * (0.5 * (P[i] + P[j]) * e + (P[j] - P[i]) * t / d)
    Fi_doc[i] += -A * (0.5 * (P[i] + P[j]) * e + (P[i] - P[j]) * t / d)
Vp = volumes(X0 + eps * v)[0]; Vm = volumes(X0 - eps * v)[0]
dVdt = (Vp - Vm) / (2 * eps)
print("dE/dt (code sign)  =", np.sum(v * Fi) - np.sum(P * dVdt), " scale", np.sum(np.abs(P * dVdt)))
print("dE/dt (doc sign)   =", np.sum(v * Fi_doc) - np.sum(P * dVdt))
print("sum F (momentum)   =", Fi.sum(axis=0))
# per-face residual (i-side + j-side) using w_f = ubar + (vi - vj).t/d
res = []
for (i, j, A, e, t, d) in F:
    if i < j:  # count each face once (periodic images may duplicate; fine for a test)
        w = 0.5 * (v[i] + v[j]) @ e + (v[i] - v[j]) @ t / d
        fi = -A * (0.5 * (P[i] + P[j]) * e + (P[j] - P[i]) * t / d)
        r = (v[i] - v[j]) @ fi - P[i] * A * w + P[j] * A * w
        res.append(r)
print("per-face residual max|.| =", np.max(np.abs(res)), " (nonzero: cancels only cell-wise)")

# 3. angular momentum (free cloud, bounded cells only)
Xf = rng.random((N, 2))
vor = Voronoi(Xf)
bounded = [k for k in range(N) if -1 not in vor.regions[vor.point_region[k]]]
Pf = np.zeros(N); Pf[bounded] = rng.random(len(bounded)) + 0.5
def vol_free(X):
    vor = Voronoi(X); V = np.zeros(N)
    for k in bounded:
        V[k] = abs(poly_area(vor.vertices[vor.regions[vor.point_region[k]]]))
    return V
Ff = np.zeros((N, 2))
for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
    if -1 in rv: continue
    for a, b in ((p, q), (q, p)):
        xi, xj = Xf[a], Xf[b]
        v1, v2 = vor.vertices[rv[0]], vor.vertices[rv[1]]
        A = np.linalg.norm(v2 - v1); c = 0.5 * (v1 + v2)
        d = np.linalg.norm(xj - xi); e = (xj - xi) / d; t = c - 0.5 * (xi + xj)
        # F_a = sum_k P_k dV_k/dx_a over this face: P_a dV_a/dx_a + P_b dV_b/dx_a
        Ff[a] += Pf[a] * A * (e / 2 + t / d) + Pf[b] * A * (-e / 2 - t / d)
# verify against FD of sum_k P_k V_k
rot = np.stack([-Xf[:, 1], Xf[:, 0]], axis=1)
Vp = vol_free(Xf + eps * rot); Vm = vol_free(Xf - eps * rot)
print("rotation: d(sum P V)/dtheta FD =", np.sum(Pf * (Vp - Vm)) / (2 * eps),
      " torque sum x x F =", np.sum(Xf[:, 0] * Ff[:, 1] - Xf[:, 1] * Ff[:, 0]))
# pair torque per face nonzero?
# 4. dissipation split: entropy per cell
Q = rng.normal(size=len(F))
TdS = np.zeros(N); dE = 0
for n, (i, j, A, e, t, d) in enumerate(F):
    ui, uj = v[i] @ e, v[j] @ e
    q = Q[n] if Q[n] * (ui - uj) >= 0 else 0.0
    TdS[i] += 0.5 * A * q * (ui - uj); TdS[j] += 0.5 * A * q * (ui - uj)
    dE += v[i] @ (-A * q * e) + v[j] @ (A * q * e) + A * q * (ui - uj)
print("dissipation: min TdS =", TdS.min(), " dE =", dE)
