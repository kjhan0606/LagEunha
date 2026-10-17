"""Small convex-polygon reference using rational coordinates; not production geometry."""
from fractions import Fraction as F


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def sub(a, b):
    return (a[0]-b[0], a[1]-b[1])


class Polygon:
    geom_type = 'Polygon'
    interiors = ()

    def __init__(self, vertices):
        p = []
        for v in vertices:
            v = tuple(F(z) for z in v)
            if not p or v != p[-1]: p.append(v)
        if len(p) > 1 and p[0] == p[-1]: p.pop()
        twice = sum((cross(a, b) for a, b in zip(p, p[1:]+p[:1])), F(0))
        if twice < 0: p.reverse(); twice = -twice
        self.vertices = p; self.exact_area = twice/2
        self.area = float(self.exact_area)
        self.is_valid = self.area > 0 and all(cross(sub(p[(i+1)%len(p)], p[i]),
                            sub(p[(i+2)%len(p)], p[(i+1)%len(p)])) >= 0 for i in range(len(p)))

    def intersection(self, other):
        p = self.vertices[:]
        for a, b in zip(other.vertices, other.vertices[1:]+other.vertices[:1]):
            if not p: break
            edge = sub(b, a); out = []
            for s, e in zip(p[-1:]+p[:-1], p):
                ds = cross(edge, sub(s, a)); de = cross(edge, sub(e, a))
                if (ds >= 0) != (de >= 0):
                    t = ds/(ds-de)
                    out.append((s[0]+t*(e[0]-s[0]), s[1]+t*(e[1]-s[1])))
                if de >= 0: out.append(e)
            p = out
        return Polygon(p)


def corner_polygons(x, triangles):
    polygons = []; labels = []
    for t, ids in enumerate(triangles):
        p = [tuple(F(z) for z in x[i]) for i in ids]
        center = tuple(sum(v[k] for v in p)/3 for k in (0, 1))
        for i in range(3):
            mid1 = tuple((p[i][k]+p[(i+1)%3][k])/2 for k in (0, 1))
            mid2 = tuple((p[i][k]+p[(i-1)%3][k])/2 for k in (0, 1))
            polygons.append(Polygon([p[i], mid1, center, mid2])); labels.append((t, ids[i]))
    return polygons, labels
