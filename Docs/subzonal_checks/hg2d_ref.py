"""Reference for GFS_HG_AFFINE (design_astra.md sec. 3-5, audit corrections): per pair face
xi = ln(1-u), u = 2 delta.Delta/d^2, delta = c_i - x_i, Delta = x_j - x_i; analytic gradient through
the centroid (vertex circumcentre Jacobians) + explicit terms; 2-hop weighted LSQ affine filter;
Hugoniot coefficient. Checks: grad vs FD, affine fields -> r = 0, energy identity, FD force."""
import numpy as np
from sub2d import Mesh, cross, R as Rv, L as Lv
G = 5/3
def cell_geom(mesh, i):
    x = mesh.x; P0 = x[i]
    cells = mesh.cells[i]; nb = mesh.nbr[i]; n = len(cells)
    a = np.array([v for v, _, _ in cells]) - P0     # corners rel. to x_i
    b = np.roll(a, -1, 0); crs = cross(a, b); V = 0.5*crs.sum()
    J = (crs[:, None]*(a + b)).sum(0)/6; c = J/V
    # centroid Jacobian wrt corners: Jc_k (2x2)
    Jc = []
    for k in range(n):
        ak, an, ap = a[k], a[(k+1) % n], a[(k-1) % n]
        ck, cp = cross(ak, an), cross(ap, ak)
        dJ = (np.outer(ak + an, Rv(an)) + ck*np.eye(2) + np.outer(ap + ak, Lv(ap)) + cp*np.eye(2))/6
        dV = 0.5*Rv(an - ap)
        Jc.append((dJ - np.outer(c, dV))/V)
    # dc/dx_b through the vertices: vertex k = circumcentre of triangle t (indices into mesh.X)
    dc = {}
    for k, (v, t, _) in enumerate(cells):
        for al in range(3):
            g = int(t[al]); Jv = mesh.vjac(t, al)
            dc[g] = dc.get(g, np.zeros((2, 2))) + Jc[k] @ Jv
    return V, c, dc, nb, a
def face_grads(mesh, i):
    V, c, dc, nb, a = cell_geom(mesh, i)
    me = 4*mesh.N + i; out = []
    n = len(nb)
    for k in range(n):
        j = nb[k]; Delta = mesh.X[j] - mesh.X[me]; d2 = Delta @ Delta
        A = np.linalg.norm(a[(k+1) % n] - a[k])
        u = 2*(c @ Delta)/d2
        g = {}
        for bb, Jb in dc.items():
            g[bb] = Jb.T @ (2*Delta/d2)
        g[me] = g.get(me, 0) + (-2*Delta/d2 - 2*c/d2 + 4*(c @ Delta)*Delta/d2**2)
        g[j] = g.get(j, 0) + (2*c/d2 - 4*(c @ Delta)*Delta/d2**2)
        for bb in g: g[bb] = -g[bb]/(1 - u)
        out.append((j, A, np.sqrt(d2), u, g))
    return V, out
def xi_all(x, Lb):
    mesh = Mesh(x, Lb); res = {}
    for i in range(len(x)):
        V, c, dc, nb, a = cell_geom(mesh, i)
        for k, j in enumerate(nb):
            Delta = mesh.X[j] - mesh.X[4*mesh.N + i]
            res[(i, int(mesh.idx[j]))] = np.log(1 - 2*(c @ Delta)/(Delta @ Delta))
    return res

def stencil(mesh, i):
    """2-hop weighted stencil of cell i: dict key=(id, rounded rel pos) -> [w, y(rel pos)]"""
    me = 4*mesh.N + i; x0 = mesh.X[me]; S = {}
    def add(idg, y, w):
        k = (int(idg), tuple(np.round(y, 9)))
        if k in S: S[k][0] += w
        else: S[k] = [w, y.copy()]
    V, c, dc, nb, a = cell_geom(mesh, i); hi = np.sqrt(V); n = len(nb); wsum = 0
    for k, j in enumerate(nb):
        y = mesh.X[j] - x0; A = np.linalg.norm(a[(k+1) % n] - a[k]); w = A/(y @ y); wsum += w
        add(mesh.idx[j], y, w)
        jj = int(mesh.idx[j])   # neighbours of j (central-tile cell jj, shifted by the image offset)
        off = mesh.X[j] - mesh.x[jj]
        cj = mesh.cells[jj]; nbj = mesh.nbr[jj]; m_ = len(nbj); aj = np.array([vv for vv, _, _ in cj]) - mesh.x[jj]
        for kk, q in enumerate(nbj):
            yq = mesh.X[q] - mesh.x[jj]; Aq = np.linalg.norm(aj[(kk+1) % m_] - aj[kk]); wq = Aq/(yq @ yq)
            y2 = yq + off + mesh.x[jj] - x0
            if np.linalg.norm(y2) < 1e-12: continue
            add(mesh.idx[q], y2, hi*w*wq)
    add(i, np.zeros(2), wsum)
    return S
def hg_terms(mesh, i, vel, rho, P, m):
    """per-face filtered operator l_b, rate r, coefficient C, force contributions, heat"""
    me = 4*mesh.N + i; S = stencil(mesh, i)
    keys = list(S.keys()); w = np.array([S[k][0] for k in keys]); y = np.array([S[k][1] for k in keys])
    yb = (w[:, None]*y).sum(0)/w.sum(); yt = y - yb; Sm = (w[:, None, None]*yt[:, :, None]*yt[:, None, :]).sum(0)
    Si = np.linalg.pinv(Sm)
    V, faces = face_grads(mesh, i)
    F = {}; H = 0.0; rates = []
    for (jj, A, d, u, g) in faces:
        # g keyed by tiled index -> (id, rel pos)
        gk = {}; yk = {}
        for bb, gb in g.items():
            k = (int(mesh.idx[bb]), tuple(np.round(mesh.X[bb] - mesh.X[me], 9)))
            gk[k] = gk.get(k, 0) + gb; yk[k] = mesh.X[bb] - mesh.X[me]
        M = sum(np.outer(gb, yk[k]) for k, gb in gk.items())
        l = {k: gk.get(k, 0) - w[q]*(M @ Si @ yt[q]) for q, k in enumerate(keys)}
        for k in gk:
            if k not in l: l[k] = gk[k]
        r = sum(l[k] @ vel[k[0]] for k in l)
        hh = d/2; W = A*d/4; rf = rho*(1 - u); cf = np.sqrt(G*P*(1 - u)**G/rf); b = (G + 1)/4*hh
        C = rf*hh*W*(b*abs(r) + np.sqrt(b*b*r*r + cf*cf))
        for k in l: F[k[0]] = F.get(k[0], 0) - C*r*l[k]
        H += C*r*r; rates.append(r)
    return F, H, np.array(rates)
