"""2D reference implementation of the subzonal potential (Docs/subzonal_derivation.md, 3.2).

Periodic box [0,L)^2, W=0 (Voronoi), ideal gas e = K rho^(g-1)/(g-1).
U_i = sum_f M_f e(rho_f, K_i), M_f = m_i B_f(c_i)/V_i, rho_f = M_f/W_f, W_f = B_f(x_i),
B_f(p) = 1/2 (a-p) x (b-p) for edge f = (a, b) (CCW).
Analytic gradient: dU_i = sum_f [-P_f dW_f + H_f dM_f], chained through the
centroid c_i = J_i/V_i and the Voronoi vertices (circumcentres of Delaunay triangles).
"""
import numpy as np
from scipy.spatial import Delaunay

G = 5.0/3.0

def cross(u, w):
    return u[..., 0]*w[..., 1] - u[..., 1]*w[..., 0]

def R(w):   # d cross(u, w)/du
    return np.array([w[1], -w[0]])

def L(u):   # d cross(u, w)/dw
    return np.array([-u[1], u[0]])


class Mesh:
    """Periodic Voronoi cells from a 3x3 tiling. cells[i] = list of (vertex xy, generator triple)."""
    def __init__(self, x, Lbox):
        self.x, self.L, self.N = x, Lbox, len(x)
        sh = np.array([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)], float)*Lbox
        X = np.concatenate([x + s for s in sh]); idx = np.tile(np.arange(self.N), 9)
        self.X, self.idx = X, idx
        tri = Delaunay(X)
        self.cells = [[] for _ in range(self.N)]
        cen = 4*self.N   # tile (0,0) is index 4
        for t in tri.simplices:
            p = X[t]
            A = 2*np.array([p[1]-p[0], p[2]-p[0]])
            b = np.array([p[1]@p[1]-p[0]@p[0], p[2]@p[2]-p[0]@p[0]])
            v = np.linalg.solve(A, b)
            for k in range(3):
                if cen <= t[k] < cen + self.N:
                    self.cells[t[k]-cen].append((v, t, k))
        # order the triangles around each generator by adjacency (robust to coincident
        # circumcentres): in CCW triangle (i, p, q) the next triangle is (i, q, r);
        # the neighbour across the Voronoi edge between them is q.
        self.nbr = []
        for i in range(self.N):
            me = cen + i; nxt = {}
            for v, t, k in self.cells[i]:
                o = [g for g in t if g != me]
                p_, q_ = (o[0], o[1]) if cross(X[o[0]] - X[me], X[o[1]] - X[me]) > 0 else (o[1], o[0])
                nxt[p_] = (q_, v, t, k)
            start = next(iter(nxt)); p_ = start; ordered = []; nb = []
            while True:
                q_, v, t, k = nxt[p_]
                ordered.append((v, t, k)); nb.append(q_)
                p_ = q_
                if p_ == start: break
            self.cells[i] = ordered; self.nbr.append(nb)

    def vjac(self, t, a_local):
        """d vertex / d X[t[a_local]] for the circumcentre of triangle t (2x2)."""
        p = self.X[t]; v = None
        A = 2*np.array([p[1]-p[0], p[2]-p[0]])
        b = np.array([p[1]@p[1]-p[0]@p[0], p[2]@p[2]-p[0]@p[0]])
        v = np.linalg.solve(A, b); Ai = np.linalg.inv(A)
        if a_local == 0:
            rhs = 2*np.outer(np.ones(2), v - p[0])
        elif a_local == 1:
            rhs = np.zeros((2, 2)); rhs[0] = 2*(p[1] - v)
        else:
            rhs = np.zeros((2, 2)); rhs[1] = 2*(p[2] - v)
        return Ai @ rhs


def cell_energy(P, verts, m, K, need_grad=False, nrm=None):
    """U_i, and optionally dU/dP and dU/dverts (n,2). nrm[f] = outward unit normal of edge f
    (direction to the neighbour); rho_f = (m/V) h_f(c)/h_f(P) stays finite on zero-length edges."""
    n = len(verts); a = verts; b = np.roll(verts, -1, axis=0)
    crs = cross(a, b); V = 0.5*crs.sum()
    J = (crs[:, None]*(a + b)).sum(0)/6.0
    c = J/V
    W = 0.5*cross(a - P, b - P); B = 0.5*cross(a - c, b - c)
    M = m*B/V
    if nrm is None:
        rho = M/W
    else:
        hx = ((a - P)*nrm).sum(1); hc = ((a - c)*nrm).sum(1)
        rho = (m/V)*hc/hx
    e = K*rho**(G-1)/(G-1); Pf = K*rho**G; H = G/(G-1)*K*rho**(G-1)
    U = (M*e).sum()
    if not need_grad:
        return U, dict(V=V, c=c, W=W, B=B, rho=rho, Pf=Pf)
    gP = np.zeros(2); gv = np.zeros((n, 2)); gc = np.zeros(2)
    Hbar_m = (H*M).sum()
    for f in range(n):
        af, bf = a[f], b[f]
        # -P_f dW_f
        gv[f] += -Pf[f]*0.5*R(bf - P)
        gv[(f+1) % n] += -Pf[f]*0.5*L(af - P)
        gP += -Pf[f]*0.5*R(af - bf)
        # (m/V) H_f dB_f
        s = m/V*H[f]
        gv[f] += s*0.5*R(bf - c)
        gv[(f+1) % n] += s*0.5*L(af - c)
        gc += s*0.5*R(af - bf)
    # -(sum H M / V) dV
    dV = np.array([0.5*R(a[(k+1) % n] - a[(k-1) % n]) for k in range(n)])
    gv += -(Hbar_m/V)*dV
    # centroid chain: dc = (dJ - c dV)/V
    for k in range(n):
        ak, an, ap = a[k], a[(k+1) % n], a[(k-1) % n]
        ck, cp = cross(ak, an), cross(ap, ak)
        dJ = (np.outer(ak + an, R(an)) + ck*np.eye(2) + np.outer(ap + ak, L(ap)) + cp*np.eye(2))/6.0
        Jc = (dJ - np.outer(c, dV[k]))/V
        gv[k] += Jc.T @ gc
    return U, dict(gP=gP, gv=gv, V=V, c=c, W=W, B=B, rho=rho, Pf=Pf)


def total(x, Lbox, m, K, grad=False):
    mesh = Mesh(x, Lbox); U = 0.0; F = np.zeros_like(x); per = []
    for i in range(len(x)):
        verts = np.array([v for v, _, _ in mesh.cells[i]])
        nrm = np.array([(mesh.X[j] - x[i])/np.linalg.norm(mesh.X[j] - x[i]) for j in mesh.nbr[i]])
        Ui, d = cell_energy(x[i], verts, m[i], K[i], need_grad=grad, nrm=nrm)
        U += Ui; per.append((Ui, d))
        if grad:
            F[i] -= d['gP']
            for k, (v, t, _) in enumerate(mesh.cells[i]):
                for al in range(3):
                    F[mesh.idx[t[al]]] -= mesh.vjac(t, al).T @ d['gv'][k]
    return U, F, per, mesh


def U0(x, Lbox, m, K):
    """single-state cell energy sum m e(rho_i)"""
    mesh = Mesh(x, Lbox); s = 0.0
    for i in range(len(x)):
        verts = np.array([v for v, _, _ in mesh.cells[i]])
        V = 0.5*cross(verts, np.roll(verts, -1, axis=0)).sum(); rho = m[i]/V
        s += m[i]*K[i]*rho**(G-1)/(G-1)
    return s
