"""Single-flip variational reference with a vertex-zero internal edge mode.

Four generators remain the only triangulation nodes. The diagonal crossing is
an internal material shape degree of freedom, not a Voronoi generator.
This fixed-reference finite element patch is not a complete remeshing solver.
"""
import numpy as np

GAMMA = 5/3


def cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


def geometry(x, triangles):
    p = x[np.asarray(triangles)]; a, b, c = p[:, 0], p[:, 1], p[:, 2]
    area = cross(b-a, c-a)/2
    if (area <= 0).any() or not np.isfinite(area).all():
        raise RuntimeError('inverted material triangle')
    d = np.stack((c-b, a-c, b-a), axis=1)
    grad = np.stack((-d[..., 1], d[..., 0]), axis=-1)/2
    return area, grad


def consistent_mass(n, triangles, mass):
    result = np.zeros((n, n))
    for ids, m in zip(triangles, mass):
        result[np.ix_(ids, ids)] += m*(np.ones((3, 3))+np.eye(3))/12
    return result


class EnrichedFlip:
    def __init__(self, x=None, rho=(1., .01), pressure=2.78e-5):
        self.x = np.array([[.7, -.1], [.9, -.1], [.9, .1], [.7, .1]]) if x is None else np.array(x, float)
        x = self.x
        # Crossing of old diagonal 0--2 and new diagonal 1--3.
        s, t = np.linalg.solve(np.column_stack((x[2]-x[0], -(x[3]-x[1]))), x[1]-x[0])
        if not (0 < s < 1 and 0 < t < 1): raise ValueError('convex quadrilateral required')
        center = (1-s)*x[0]+s*x[2]
        self.old_lam = np.array([1-s, 0., s, 0.])
        self.new_lam = np.array([0., 1-t, 0., t])
        self.T = np.eye(5); self.T[4, :4] = self.new_lam
        self.Q0 = np.vstack((x, center-self.new_lam@x))
        self.old_tri = [(0, 1, 2), (0, 2, 3)]
        self.tri = [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]
        area, _ = geometry(self.T@self.Q0, self.tri)
        self.area0 = area
        self.mass = area*np.array([rho[0], rho[0], rho[1], rho[1]])
        self.entropy = pressure/(self.mass/area)**GAMMA
        self.pressure = pressure
        self.Mphysical = consistent_mass(5, self.tri, self.mass)
        self.M = self.T.T@self.Mphysical@self.T
        oldarea, _ = geometry(x, self.old_tri)
        self.Mold = consistent_mass(4, self.old_tri, oldarea*np.asarray(rho))

    def encode_velocity(self, velocity):
        # Exact old piecewise-linear trace at the diagonal crossing.
        internal = (self.old_lam-self.new_lam)@velocity
        return np.vstack((velocity, internal))

    def state(self, Q, V):
        y = self.T@Q; velocity = self.T@V
        canonical = self.M@V
        momentum = canonical[:4].sum(0)
        angular = np.sum(cross(Q, canonical))
        kinetic = .5*np.sum(V*canonical)
        area, _ = geometry(y, self.tri)
        internal = self.entropy*self.mass**GAMMA*area**(1-GAMMA)/(GAMMA-1)
        return dict(momentum=momentum, angular=angular, kinetic=kinetic,
                    internal=internal.sum(), augmented_energy=kinetic+internal.sum()+self.pressure*area.sum(),
                    area=area, physical_position=y, physical_velocity=velocity)

    def potential(self, Q):
        area, _ = geometry(self.T@Q, self.tri)
        return np.sum(self.entropy*self.mass**GAMMA*area**(1-GAMMA)/(GAMMA-1)+self.pressure*area)

    def force(self, Q):
        y = self.T@Q; area, grad = geometry(y, self.tri)
        pressure = self.entropy*(self.mass/area)**GAMMA
        physical = np.zeros((5, 2))
        np.add.at(physical, np.asarray(self.tri).ravel(), ((pressure-self.pressure)[:, None, None]*grad).reshape(-1, 2))
        return self.T.T@physical

    def acceleration(self, Q):
        return np.linalg.solve(self.M, self.force(Q))

    def residual_value(self, point, Q, V):
        """Old nodal P1 field minus the new P1 field with unchanged nodal speeds.

        Evaluate directly on the old/new triangles, for the exact trace test.
        Q is only used for the original four positions in this fixed-geometry test.
        """
        x = Q[:4]
        def interpolate(triangles):
            for ids in triangles:
                lam = np.linalg.solve(np.vstack((np.ones(3), x[list(ids)].T)), np.r_[1., point])
                if lam.min() >= -1e-12:
                    return lam@V[list(ids)]
            raise ValueError('point outside patch')
        return interpolate(self.old_tri)-interpolate([(0, 1, 3), (1, 2, 3)])
