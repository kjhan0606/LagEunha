import numpy as np
from scipy.spatial import Voronoi
g = 5.0/3.0; K = 1.0
def e(rho): return K*rho**(g-1)/(g-1)

L = 4.0; n = 4; N = n*n
def cells(x):
    """periodic Voronoi: returns per-cell list of (A_f, h_x, h_c, j) and V_i, c_i."""
    shifts = np.array([[a,b] for a in (-1,0,1) for b in (-1,0,1)])*L
    pts = np.concatenate([x + s for s in shifts])  # index k*N + i, k=4 is center
    vor = Voronoi(pts)
    out = []
    for i in range(N):
        idx = 4*N + i
        reg = vor.regions[vor.point_region[idx]]
        verts = vor.vertices[reg]
        # order by angle around generator
        ang = np.arctan2(verts[:,1]-x[i,1], verts[:,0]-x[i,0]); verts = verts[np.argsort(ang)]
        xs, ys = verts[:,0], verts[:,1]
        cross = xs*np.roll(ys,-1) - np.roll(xs,-1)*ys
        V = 0.5*cross.sum()
        c = np.array([((xs+np.roll(xs,-1))*cross).sum(), ((ys+np.roll(ys,-1))*cross).sum()])/(6*V)
        faces = []
        for (p,q), rv in zip(vor.ridge_points, vor.ridge_vertices):
            if idx in (p,q) and -1 not in rv:
                j = q if p==idx else p
                v1, v2 = vor.vertices[rv]
                A = np.linalg.norm(v2-v1)
                dvec = pts[j]-pts[idx]; d = np.linalg.norm(dvec); ev = dvec/d
                mid = 0.5*(pts[idx]+pts[j])
                h_x = d/2; h_c = ev @ (mid - c)
                faces.append((A, h_x, h_c, j % N))
        out.append((V, c, faces))
    return out
def energies(x, m=1.0):
    Us = []; Phis = []
    for V, c, faces in cells(x):
        rho = m/V
        sumVc = sum(0.5*A*hc for A,_,hc,_ in faces)
        assert abs(sumVc - V) < 1e-9*V
        U = sum(0.5*A*hc*rho * e(rho*hc/hx) for A,hx,hc,_ in faces)
        Us.append(U); Phis.append(U - m*e(rho))
    return np.array(Us), np.array(Phis)
def Utot(x): return energies(x)[0].sum()
def neighbors(x):
    return [frozenset(j for *_, j in f) for _,_,f in cells(x)]

rng = np.random.default_rng(3)
x0 = (np.array([[i,j] for i in range(n) for j in range(n)]) + 0.5)*L/n
# 1) lattice: centroid = generator, Phi = 0
U, Phi = energies(x0)
print("lattice: max|Phi| =", np.abs(Phi).max())
# sheared lattice (uniform shear of lattice stays a lattice) -> Phi still 0
xs = x0.copy(); xs[:,0] += 0.3*xs[:,1]; xs[:,0] %= L
U, Phi = energies(xs); print("sheared lattice: max|Phi| =", np.abs(Phi).max(), " (not exactly periodic-compatible shear; just indicative)")
# 2) random perturbation: Phi >= 0
xr = (x0 + rng.uniform(-0.25,0.25,x0.shape)*L/n) % L
U, Phi = energies(xr)
print("perturbed: Phi min =", Phi.min(), " max =", Phi.max(), " sum Phi / sum m e(rho) =", Phi.sum()/(U.sum()-Phi.sum()))
# 3) momentum conservation: sum_i F_i = 0 via finite differences
def grad(x, h=1e-6):
    G = np.zeros_like(x)
    for a in range(N):
        for k in range(2):
            xp = x.copy(); xp[a,k]+=h; xm = x.copy(); xm[a,k]-=h
            G[a,k] = (Utot(xp)-Utot(xm))/(2*h)
    return G
G = grad(xr); print("sum of forces (should be 0):", -G.sum(axis=0), " |F| typical:", np.abs(G).mean())
# 4) continuity / C1 across a flip: move generator 0 along x, find neighbor-set change
nb0 = neighbors(xr)
s_vals = np.linspace(0, 0.5*L/n, 2001)
flips = []
prev = nb0[0]
Uvals = []
for s in s_vals:
    xt = xr.copy(); xt[0,0] = (xr[0,0]+s) % L
    nb = neighbors(xt)[0]
    Uvals.append(Utot(xt))
    if nb != prev:
        flips.append(s); prev = nb
Uvals = np.array(Uvals)
print("flips found at s =", flips[:4])
if flips:
    s_f = flips[0]
    def U_at(s):
        xt = xr.copy(); xt[0,0] = (xr[0,0]+s) % L; return Utot(xt)
    # bracket the flip precisely by bisection on neighbor set
    lo, hi = s_f - s_vals[1], s_f
    for _ in range(50):
        mid = 0.5*(lo+hi)
        xt = xr.copy(); xt[0,0] = (xr[0,0]+mid) % L
        if neighbors(xt)[0] == nb0[0] if flips[0]==flips[0] else False: lo = mid
        else: hi = mid
    sf = 0.5*(lo+hi)
    eps = 1e-5
    U_left, U_right = U_at(sf-eps), U_at(sf+eps)
    dU_left = (U_at(sf-eps)-U_at(sf-2*eps))/eps
    dU_right = (U_at(sf+2*eps)-U_at(sf+eps))/eps
    print(f"flip at s={sf:.8f}: U jump = {U_right-U_left:.3e} (|U|={U_left:.3e})")
    print(f"  dU/ds left = {dU_left:.6e}, right = {dU_right:.6e}, jump = {dU_right-dU_left:.3e}")
    # compare with the cell-volume Lagrangian gradient continuity
    def Ucell_at(s):
        xt = xr.copy(); xt[0,0] = (xr[0,0]+s) % L
        Us, Ph = energies(xt); return (Us-Ph).sum()
    dl = (Ucell_at(sf-eps)-Ucell_at(sf-2*eps))/eps; dr = (Ucell_at(sf+2*eps)-Ucell_at(sf+eps))/eps
    print(f"  cell-volume U: dU/ds left = {dl:.6e}, right = {dr:.6e}, jump = {dr-dl:.3e}")
    # second derivative jump of cell-volume U (expected: C1 but not C2)
    d2l = (Ucell_at(sf-eps)-2*Ucell_at(sf-2*eps)+Ucell_at(sf-3*eps))/eps**2
    d2r = (Ucell_at(sf+3*eps)-2*Ucell_at(sf+2*eps)+Ucell_at(sf+eps))/eps**2
    print(f"  cell-volume U: d2U left={d2l:.4e} right={d2r:.4e}")
