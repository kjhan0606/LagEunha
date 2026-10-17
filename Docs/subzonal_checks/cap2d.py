"""2D free-surface cell with a rigid cap box (side L = 5 Delta): does the subzonal Phi
pull the edge generator toward vacuum, and how does it compare with the cell-volume
force and with the physical free-surface force P0*Delta?  Finite-difference forces of
the total potential over all cells (Sutherland-Hodgman clipping of the cap square by
the bisector half-planes).  Variants: full (cap fans included), A (perimeter weight),
B (cap fan mass 0), P (reference displacement projected tangentially to the cap).
"""
import numpy as np
G = 5.0/3.0
D = 1.0            # lattice spacing
L = 5.0*D          # cap box side

def clip(poly, n, b):
    """keep points with n.r <= b"""
    out = []
    m = len(poly)
    for k in range(m):
        p, q = poly[k], poly[(k+1) % m]
        fp, fq = n@p - b, n@q - b
        if fp <= 0: out.append(p)
        if fp*fq < 0:
            t = fp/(fp - fq); out.append(p + t*(q - p))
    return out

def cell(i, X):
    xi = X[i]
    poly = [xi + np.array(s)*L/2 for s in [(-1, -1), (1, -1), (1, 1), (-1, 1)]]
    kinds = ['cap']*4
    for j in range(len(X)):
        if j == i: continue
        d = X[j] - xi; dd = d@d
        if dd > 2*L*L: continue
        n = d/np.sqrt(dd); b = n@xi + 0.5*np.sqrt(dd)
        poly = clip(poly, n, b)
    out = []
    for p in poly:
        if not out or np.linalg.norm(p - out[-1]) > 1e-10: out.append(p)
    if len(out) > 1 and np.linalg.norm(out[0] - out[-1]) < 1e-10: out.pop()
    return np.array(out)

def edge_kind(poly, i, X):
    """for each edge (a,b): 'pair' with outward normal to neighbour, or 'cap'"""
    xi = X[i]; m = len(poly); kinds = []
    for k in range(m):
        a, b = poly[k], poly[(k+1) % m]
        mid = 0.5*(a + b); t = b - a
        n = np.array([t[1], -t[0]]); n /= np.linalg.norm(n)   # outward for CCW
        h = (mid - xi)@n
        # cap edge: lies on the box, i.e. |(mid-xi)_x| or |..._y| = L/2 and axis-aligned
        if abs(abs((mid - xi)[0]) - L/2) < 1e-9 and abs(t[0]) < 1e-9 or abs(abs((mid - xi)[1]) - L/2) < 1e-9 and abs(t[1]) < 1e-9:
            kinds.append(('cap', n))
        else:
            kinds.append(('pair', n))
    return kinds

def e_of(rho, K): return K*rho**(G-1)/(G-1)

def Ucell(i, X, m, K, variant):
    poly = cell(i, X); xi = X[i]
    a = poly; b = np.roll(poly, -1, axis=0)
    crs = a[:, 0]*b[:, 1] - a[:, 1]*b[:, 0]
    V = 0.5*crs.sum(); c = (crs[:, None]*(a + b)).sum(0)/(6*V)
    kinds = edge_kind(poly, i, X)
    iscap = np.array([k[0] == 'cap' for k in kinds])
    W = 0.5*((a - xi)[:, 0]*(b - xi)[:, 1] - (a - xi)[:, 1]*(b - xi)[:, 0])
    U0 = m*e_of(m/V, K)
    if variant == 'cell':
        return U0, dict(V=V, c=c, Lcap=0)
    if variant == 'P':
        # project the reference displacement tangentially to the cap normals
        cc = c.copy()
        for (kind, n) in kinds:
            if kind == 'cap':
                cc = xi + (cc - xi) - n*((cc - xi)@n)
        c = cc
    B = 0.5*((a - c)[:, 0]*(b - c)[:, 1] - (a - c)[:, 1]*(b - c)[:, 0])
    if variant == 'B':
        Bp = np.where(iscap, 0, B); M = m*Bp/Bp.sum()
        sel = ~iscap
        U = (M[sel]*e_of(M[sel]/W[sel], K)).sum()
        return U, dict(V=V, c=c, Lcap=0)
    M = m*B/V
    U = (M*e_of(M/W, K)).sum()
    Lcap = sum(np.linalg.norm(b[k] - a[k]) for k in range(len(a)) if iscap[k])
    Lper = sum(np.linalg.norm(b[k] - a[k]) for k in range(len(a)))
    if variant == 'A':
        w = 1 - Lcap/Lper
        U = U0 + w*(U - U0)
    return U, dict(V=V, c=c, Lcap=Lcap, phi=U - U0, U0=U0, W=W, M=M, iscap=iscap)

def Utot(X, m, K, variant, idx=None):
    idx = range(len(X)) if idx is None else idx
    return sum(Ucell(i, X, m, K, variant)[0] for i in idx)

def force(X, i, m, K, variant, h=1e-6):
    F = np.zeros(2)
    # cells affected: i and those within 2L
    near = [j for j in range(len(X)) if np.linalg.norm(X[j] - X[i]) < 2*L + 1e-9]
    for d in range(2):
        Xp = X.copy(); Xp[i, d] += h; Xm = X.copy(); Xm[i, d] -= h
        F[d] = -(Utot(Xp, m, K, variant, near) - Utot(Xm, m, K, variant, near))/(2*h)
    return F

if __name__ == '__main__':
    # lattice: columns x = -4..4, rows y = -4..0 ; row y=0 is the free surface (vacuum above)
    import sys
    hexlat = len(sys.argv) > 1 and sys.argv[1] == 'hex'
    if hexlat:
        ry = np.sqrt(3)/2*D
        pts = [((ix + 0.5*(iy % 2))*D, iy*ry) for ix in range(-5, 6) for iy in range(-5, 1)]
        pts = [(x, y) for (x, y) in pts]
        A0 = D*ry
    else:
        pts = [(ix*D, iy*D) for ix in range(-4, 5) for iy in range(-4, 1)]; A0 = D*D
    X = np.array(pts, float); N = len(X)
    m = 1.0; K = 1.0; rho0 = m/A0; P0 = K*rho0**G
    i0 = [k for k in range(N) if abs(X[k, 0]) < 1e-9 and abs(X[k, 1]) < 1e-9][0]
    print("interior P0 =", P0, " physical surface force P0*Delta =", P0*D)
    for var in ['cell', 'full', 'A', 'B', 'P']:
        U, d = Ucell(i0, X, m, K, var)
        F = force(X, i0, m, K, var)
        extra = ""
        if 'phi' in d: extra = " Phi/U0=%.3f Lcap=%.2f" % (d['phi']/d['U0'], d['Lcap'])
        print("%5s: edge-cell rho/rho0=%.3f  c-x=(%.2f,%.2f)  F/(P0 D)=(%.4f, %.4f)%s"
              % (var, m/d['V']/rho0, *(d['c'] - X[i0]), *(F/(P0*D)), extra))
    # tangential hourglass stiffness at the surface: displace edge generator in x
    print("\n tangential stiffness k_t = d^2U/dx^2 /(P0) (surface cell), and normal k_n:")
    for var in ['cell', 'full', 'P']:
        h = 1e-2
        def U_at(dx, dy):
            Xp = X.copy(); Xp[i0] += (dx, dy)
            near = [j for j in range(N) if np.linalg.norm(X[j] - X[i0]) < 2*L + 1e-9]
            return Utot(Xp, m, K, var, near)
        kt = (U_at(h, 0) - 2*U_at(0, 0) + U_at(-h, 0))/h**2
        kn = (U_at(0, h) - 2*U_at(0, 0) + U_at(0, -h))/h**2
        h2 = h/10
        kt2 = (U_at(h2, 0) - 2*U_at(0, 0) + U_at(-h2, 0))/h2**2
        print("%5s: k_t/P0 = %.4f (h=1e-2) %.4f (h=1e-3)  k_n/P0 = %.4f" % (var, kt/P0, kt2/P0, kn/P0))
    # interior cell for comparison
    i1 = [k for k in range(N) if abs(X[k, 0]) < 1e-9 and abs(X[k, 1] + 2*D) < 1e-9][0]
    for var in ['cell', 'full']:
        h = 1e-3
        def U_at(dx, dy):
            Xp = X.copy(); Xp[i1] += (dx, dy)
            near = [j for j in range(N) if np.linalg.norm(X[j] - X[i1]) < 2*L + 1e-9]
            return Utot(Xp, m, K, var, near)
        kt = (U_at(h, 0) - 2*U_at(0, 0) + U_at(-h, 0))/h**2
        print("interior %5s: k/P0 = %.4f" % (var, kt/P0))
    # continuity of U_i when the cap appears: outer neighbour at distance s crossing L
    print("\n continuity across cap appearance (1 cell, neighbour below at D, neighbour above at s):")
    for s in [4.8, 4.95, 4.999, 5.001, 5.05, 5.2]:
        Xc = np.array([(0, 0), (0, -D), (0, s), (-D, 0), (D, 0), (-D, -D), (D, -D), (-D, s), (D, s)], float)
        U, d = Ucell(0, Xc, m, K, 'full')
        print("  s=%.3f  V=%.4f  Phi=%.5f  Lcap=%.3f" % (s, d['V'], d['phi'], d['Lcap']))
