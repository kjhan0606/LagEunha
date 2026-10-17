"""Reference for Astra audit 18: Voronoi cells closed by (a) the common material polygon D_Gamma
(inner loop + outer loop through boundary generators) or (b) the per-generator disc B(x_i, Rc).
Exact volume gradient G_ak = dV_a/dx_k by the shape derivative:
  faces:  G_aa += A/d (c_f - x_a),  G_aj -= A/d (c_f - x_j)        (c_f = midpoint of the face piece)
  polygon boundary piece on edge e=(b,b+1), params [s0,s1]:  G_{a,Ib} += L n [(s1-s0) - (s1^2-s0^2)/2],
                                                              G_{a,Ib+1} += L n (s1^2-s0^2)/2
  disc arc (moves rigidly with x_a):  G_aa += integral n ds  (= chord normal * chord length per piece)"""
import numpy as np
from scipy.spatial import Delaunay
import shapely
from shapely.geometry import Polygon, MultiPolygon

def J(v): return np.array([-v[1], v[0]])

class Geom:
    def __init__(self, X, mode, inner=None, outer=None, Rc=None, nseg=256):
        self.X = X; self.mode = mode; N = len(X); self.nseg = nseg
        tri = Delaunay(X); ip, ix = tri.vertex_neighbor_vertices
        self.nb = [ix[ip[a]:ip[a+1]] for a in range(N)]
        self.V = np.zeros(N); self.faces = [[] for _ in range(N)]; self.bnd = [[] for _ in range(N)]
        if mode == 'poly':
            self.inner = inner; self.outer = outer
            D = Polygon(X[outer], holes=[X[inner][::-1]])
            self.edges = []   # (b0 id, b1 id) along loops, oriented with fluid on the left
            for loop in (outer, inner):   # outer CCW, inner CW (as given)
                for q in range(len(loop)): self.edges.append((loop[q], loop[(q+1) % len(loop)]))
        B = 8.0
        for a in range(N):
            cell = Polygon([(-B, -B), (B, -B), (B, B), (-B, B)])
            for j in self.nb[a]:
                dvec = X[j] - X[a]; mid = 0.5*(X[j] + X[a])
                # half-plane (p - mid).dvec <= 0 as big polygon
                t = J(dvec)/np.linalg.norm(dvec)*50; n = dvec/np.linalg.norm(dvec)*50
                hp = Polygon([mid + t, mid - t, mid - t - n, mid + t - n])
                cell = cell.intersection(hp)
            if mode == 'poly': cell = cell.intersection(D)
            else:
                circ = Polygon([X[a] + Rc*np.array([np.cos(th), np.sin(th)]) for th in np.linspace(0, 2*np.pi, nseg, endpoint=False)])
                cell = cell.intersection(circ)
            self.V[a] = cell.area
            polys = [cell] if isinstance(cell, Polygon) else [g for g in getattr(cell, 'geoms', []) if isinstance(g, Polygon)]
            for pg in polys:
                pg = shapely.geometry.polygon.orient(pg, 1.0)
                cs = np.array(pg.exterior.coords)
                for k in range(len(cs) - 1):
                    p, q = cs[k], cs[k+1]; L = np.linalg.norm(q - p)
                    if L < 1e-14: continue
                    self._classify(a, p, q, L, Rc)
    def _classify(self, a, p, q, L, Rc):
        X = self.X; tol = 1e-9
        for j in self.nb[a]:
            dvec = X[j] - X[a]; dn = np.linalg.norm(dvec); e = dvec/dn; mid = 0.5*(X[j] + X[a])
            if abs((p - mid) @ e) < tol and abs((q - mid) @ e) < tol:
                self.faces[a].append((j, L, 0.5*(p + q), dn)); return
        if self.mode == 'poly':
            for ei, (b0, b1) in enumerate(self.edges):
                y0, y1 = X[b0], X[b1]; t = y1 - y0; Le = np.linalg.norm(t); n = -J(t)/Le
                if abs((p - y0) @ n) < tol and abs((q - y0) @ n) < tol:
                    s0 = (p - y0) @ t/Le**2; s1 = (q - y0) @ t/Le**2
                    if s1 < s0: s0, s1 = s1, s0
                    self.bnd[a].append(('e', b0, b1, s0, s1, Le, n)); return
        else:
            rmin = Rc*np.cos(np.pi/self.nseg) - 1e-9
            if np.linalg.norm(p - X[a]) >= rmin and np.linalg.norm(q - X[a]) >= rmin and np.linalg.norm(0.5*(p + q) - X[a]) >= rmin:
                n = -J(q - p)/L   # outward for CCW
                self.bnd[a].append(('arc', L*n)); return
        raise RuntimeError('unclassified segment of cell %d' % a)
    def grad(self):
        """sparse G: dict per cell a -> {k: 2-vector}"""
        G = []
        X = self.X
        for a in range(len(X)):
            g = {}
            for (j, A, c, d) in self.faces[a]:
                g[a] = g.get(a, 0) + A/d*(c - X[a]); g[j] = g.get(j, 0) - A/d*(c - X[j])
            for b in self.bnd[a]:
                if b[0] == 'e':
                    _, b0, b1, s0, s1, Le, n = b
                    g[b0] = g.get(b0, 0) + Le*n*((s1 - s0) - 0.5*(s1*s1 - s0*s0))
                    g[b1] = g.get(b1, 0) + Le*n*0.5*(s1*s1 - s0*s0)
                else:
                    g[a] = g.get(a, 0) + b[1]
            G.append(g)
        return G
