"""Reference geometry for disc-clipped Voronoi cells (Astra audit 15): cell_i = disc(x_i, R) intersected with the
Voronoi half-planes. Exact area, centroid (with circular-arc moments), arc components, arc subzone
xi = log(B(c)/W), and analytic derivatives with respect to generator positions, checked against FD.
All geometry is in coordinates relative to x_i, so everything is a function of D_j = x_j - x_i."""
import numpy as np

def J(v): return np.array([-v[1], v[0]])
def cross(a, b): return a[0]*b[1] - a[1]*b[0]

def clipped_cell(D, R):
    """D: (k,2) neighbour offsets (all candidates within 2R). Returns boundary as a CCW list of pieces:
    ('seg', p, q, info_p, info_q) or ('arc', th0, th1, j0, j1); info = ('v', j, k) Voronoi vertex of bisectors j,k
    or ('e', j) bisector-circle point. Polygon by clipping a box with half-planes, then circle intersection."""
    # half-plane j: p.D_j <= |D_j|^2/2
    B = 4*R
    poly = [(np.array([-B, -B]), ('box',)), (np.array([B, -B]), ('box',)), (np.array([B, B]), ('box',)), (np.array([-B, B]), ('box',))]
    # store vertices with labels of the two lines meeting there; edge k from poly[k] to poly[k+1] lies on line lab[k]
    lines = ['box0', 'box1', 'box2', 'box3']
    for j, d in enumerate(D):
        c0 = 0.5*(d @ d); newp = []; newl = []
        n = len(poly)
        for k in range(n):
            p, lp = poly[k]; q, lq = poly[(k+1) % n]; L = lines[k]
            fp = p @ d - c0; fq = q @ d - c0
            if fp <= 0: newp.append((p, lp)); newl.append(L)
            if fp*fq < 0:
                t = fp/(fp - fq); x = p + t*(q - p)
                if fp <= 0:   # leaving: new vertex on L and j, next edge on j
                    newp.append((x, ('v', L, j))); newl.append(j)
                else:         # entering: vertex of j and L, next edge on L
                    newp.append((x, ('v', j, L))); newl.append(L)
        poly = newp; lines = newl
        if not poly: return None
    # now intersect convex polygon with circle |p| <= R
    n = len(poly); pieces = []
    inside = [np.hypot(*p) < R for p, _ in poly]
    if all(inside):
        for k in range(n): pieces.append(('seg', poly[k][0], poly[(k+1) % n][0], poly[k][1], poly[(k+1) % n][1], lines[k]))
        return pieces
    # collect boundary as sequence: segments inside and arcs
    events = []   # list of (kind, data) along polygon order
    for k in range(n):
        p, lp = poly[k]; q, lq = poly[(k+1) % n]; L = lines[k]; dpq = q - p
        a = dpq @ dpq; b = 2*(p @ dpq); c = p @ p - R*R; disc = b*b - 4*a*c
        ts = []
        if disc > 0:
            s = np.sqrt(disc); ts = [t for t in ((-b - s)/(2*a), (-b + s)/(2*a)) if 0 < t < 1]
        pts = [(0.0, p, lp, inside[k])] + [(t, p + t*dpq, ('e', L), None) for t in ts] + [(1.0, q, lq, inside[(k+1) % n])]
        # walk sub-segments
        for m in range(len(pts) - 1):
            t0, x0, l0, _ = pts[m]; t1, x1, l1, _ = pts[m+1]
            mid = 0.5*(x0 + x1)
            if np.hypot(*mid) < R: events.append(('seg', x0, x1, l0, l1, L))
    # arcs: connect each segment end that lies on the circle (label 'e') to the next segment start on the circle
    segs = events
    out = []
    for m, sg in enumerate(segs):
        out.append(sg)
        nxt = segs[(m+1) % len(segs)]
        if sg[4][0] == 'e':
            th0 = np.arctan2(sg[2][1], sg[2][0]); th1 = np.arctan2(nxt[1][1], nxt[1][0])
            while th1 <= th0: th1 += 2*np.pi
            out.append(('arc', th0, th1, sg[4][1], nxt[3][1]))
    if not segs:   # whole disc
        out = [('arc', 0.0, 2*np.pi, None, None)]
    return out

def area_moment(pieces, R):
    V = 0.0; M = np.zeros(2)
    for pc in pieces:
        if pc[0] == 'seg':
            p, q = pc[1], pc[2]; cr = cross(p, q); V += 0.5*cr; M += (p + q)*cr/6
        else:
            th0, th1 = pc[1], pc[2]; V += 0.5*R*R*(th1 - th0)
            M += (R**3/3)*np.array([np.sin(th1) - np.sin(th0), np.cos(th0) - np.cos(th1)])
    return V, M

def arc_xi(pieces, R):
    V, M = area_moment(pieces, R); c = M/V; out = []
    for pc in pieces:
        if pc[0] != 'arc': continue
        th0, th1 = pc[1], pc[2]
        W = 0.5*R*R*(th1 - th0); N = R*np.array([np.sin(th1) - np.sin(th0), np.cos(th0) - np.cos(th1)])
        Bc = W - 0.5*(c @ N)
        out.append((pc[3], pc[4], np.log(Bc/W), W, N, Bc))
    return V, c, out

def _pt(lab, D, dD, R, xguess):
    """point and tangent for a boundary label"""
    if lab[0] == 'v':
        j, k = lab[1], lab[2]
        A = np.array([D[j], D[k]]); b = 0.5*np.array([D[j] @ D[j], D[k] @ D[k]]); p = np.linalg.solve(A, b)
        rhs = np.array([(D[j] - p) @ dD[j], (D[k] - p) @ dD[k]]); return p, np.linalg.solve(A, rhs)
    if lab[0] == 'e':
        j = lab[1]; th = np.arctan2(xguess[1], xguess[0])
        # refine th on 2R n.D = |D|^2 near guess
        for _ in range(50):
            n = np.array([np.cos(th), np.sin(th)]); f = 2*R*(n @ D[j]) - D[j] @ D[j]; fp = 2*R*(J(n) @ D[j]); th -= f/fp
        n = np.array([np.cos(th), np.sin(th)]); dth = ((D[j] - R*n) @ dD[j])/(R*(J(n) @ D[j]))
        return R*n, R*J(n)*dth
    raise ValueError(lab)

def _th(lab, D, dD, R, thguess):
    j = lab
    th = thguess
    for _ in range(50):
        n = np.array([np.cos(th), np.sin(th)]); f = 2*R*(n @ D[j]) - D[j] @ D[j]; fp = 2*R*(J(n) @ D[j]); th -= f/fp
    n = np.array([np.cos(th), np.sin(th)]); return th, ((D[j] - R*n) @ dD[j])/(R*(J(n) @ D[j]))

def tangent(pieces, D, dD, R):
    """forward-mode derivative of (V, c, arc xi list) along the neighbour-offset perturbation dD (fixed topology)"""
    V = 0.0; dV = 0.0; M = np.zeros(2); dM = np.zeros(2); arcs = []
    for pc in pieces:
        if pc[0] == 'seg':
            p, dp = _pt(pc[3], D, dD, R, pc[1]); q, dq = _pt(pc[4], D, dD, R, pc[2])
            cr = cross(p, q); dcr = cross(dp, q) + cross(p, dq)
            V += 0.5*cr; dV += 0.5*dcr; M += (p + q)*cr/6; dM += ((dp + dq)*cr + (p + q)*dcr)/6
        else:
            th0, d0 = _th(pc[3], D, dD, R, pc[1]); th1, d1 = _th(pc[4], D, dD, R, pc[2])
            while th1 <= th0: th1 += 2*np.pi
            V += 0.5*R*R*(th1 - th0); dV += 0.5*R*R*(d1 - d0)
            M += (R**3/3)*np.array([np.sin(th1) - np.sin(th0), np.cos(th0) - np.cos(th1)])
            dM += (R**3/3)*np.array([np.cos(th1)*d1 - np.cos(th0)*d0, np.sin(th1)*d1 - np.sin(th0)*d0])
            arcs.append((th0, th1, d0, d1))
    c = M/V; dc = (dM - c*dV)/V; out = []
    for th0, th1, d0, d1 in arcs:
        W = 0.5*R*R*(th1 - th0); dW = 0.5*R*R*(d1 - d0)
        N = R*np.array([np.sin(th1) - np.sin(th0), np.cos(th0) - np.cos(th1)])
        dN = R*np.array([np.cos(th1)*d1 - np.cos(th0)*d0, np.sin(th1)*d1 - np.sin(th0)*d0])
        Bc = W - 0.5*(c @ N); dB = dW - 0.5*(dc @ N + c @ dN)
        out.append((np.log(Bc/W), dB/Bc - dW/W))
    return V, dV, c, dc, out
