"""Sedov-Taylor point explosion, exact self-similar solution.

Uniform ambient density rho0, zero ambient pressure, energy E released at
t=0, geometry j=1 (planar), 2 (cylindrical), 3 (spherical).

Self-similar variables xi=r/R(t), R = xi0 (E t^2/rho0)^(1/(j+2)),
Rdot = delta R/t, delta = 2/(j+2):
    u = Rdot U(xi),  rho = rho0 G(xi),  p = rho0 Rdot^2 P(xi).
Mass, momentum and entropy equations become three linear ODEs for
(G', U', P'), integrated from the strong-shock state at xi=1 inward:
    (U-xi) G' + G U' = -(j-1) G U/xi
    (U-xi) U' + P'/G = -((delta-1)/delta) U
   -gamma (U-xi) G'/G + (U-xi) P'/P = -2 (delta-1)/delta
    U(1)=2/(g+1), G(1)=(g+1)/(g-1), P(1)=2/(g+1).
Energy fixes xi0: 1 = delta^2 xi0^(j+2) * A_j int (G U^2/2 + P/(g-1)) xi^(j-1) dxi.

Checks (selftest): spherical gamma=5/3 -> xi0 = 1.1517, gamma=1.4 -> 1.0328
(standard values, e.g. Kamm & Timmes 2007; Landau & Lifshitz sec. 106).
"""
import numpy as np
from scipy.integrate import solve_ivp, simpson

_AJ = {1: 2.0, 2: 2.0 * np.pi, 3: 4.0 * np.pi}


def _similarity(gamma, j=3, xi_min=1e-4, n=4000):
    g = gamma
    delta = 2.0 / (j + 2)
    a = (delta - 1.0) / delta

    def rhs(xi, y):
        G, U, P = y
        w = U - xi
        # A @ [G', U', P'] = b
        A = np.array([[w, G, 0.0],
                      [0.0, w, 1.0 / G],
                      [-g * w / G, 0.0, w / P]])
        b = np.array([-(j - 1) * G * U / xi, -a * U, -2.0 * a])
        return np.linalg.solve(A, b)

    y0 = [(g + 1) / (g - 1), 2 / (g + 1), 2 / (g + 1)]
    xs = np.concatenate([np.linspace(1.0, 0.2, n // 2, endpoint=False),
                         np.geomspace(0.2, xi_min, n // 2)])
    sol = solve_ivp(rhs, (1.0, xi_min), y0, t_eval=xs, method="LSODA",
                    rtol=1e-10, atol=1e-14)
    if not sol.success:
        raise RuntimeError(sol.message)
    xi = sol.t[::-1]
    G, U, P = sol.y[:, ::-1]
    integrand = (0.5 * G * U ** 2 + P / (g - 1)) * xi ** (j - 1)
    I = _AJ[j] * simpson(integrand, x=xi)
    xi0 = (1.0 / (delta ** 2 * I)) ** (1.0 / (j + 2))
    return xi, G, U, P, xi0, delta


_cache = {}


def sedov_xi0(gamma, j=3):
    key = (round(gamma, 12), j)
    if key not in _cache:
        _cache[key] = _similarity(gamma, j)
    return _cache[key][4]


def shock_radius(t, gamma=5.0 / 3.0, E=1.0, rho0=1.0, j=3):
    return sedov_xi0(gamma, j) * (E * t * t / rho0) ** (1.0 / (j + 2))


def sedov_profile(r, t, gamma=5.0 / 3.0, E=1.0, rho0=1.0, p0=0.0, j=3):
    """rho, v_r, P at radii r (array). Outside the shock: rho0, 0, p0."""
    key = (round(gamma, 12), j)
    if key not in _cache:
        _cache[key] = _similarity(gamma, j)
    xi_t, G, U, P, xi0, delta = _cache[key]
    R = xi0 * (E * t * t / rho0) ** (1.0 / (j + 2))
    Rdot = delta * R / t
    r = np.asarray(r, dtype=float)
    xi = r / R
    rho = np.full_like(r, rho0)
    v = np.zeros_like(r)
    p = np.full_like(r, p0)
    ins = xi <= 1.0
    xc = np.clip(xi[ins], xi_t[0], 1.0)
    rho[ins] = rho0 * np.interp(xc, xi_t, G)
    v[ins] = Rdot * np.interp(xc, xi_t, U)
    p[ins] = rho0 * Rdot ** 2 * np.interp(xc, xi_t, P) + 0 * p0
    return rho, v, p, R


if __name__ == "__main__":
    for g in (5.0 / 3.0, 1.4):
        print("gamma=%.4f  xi0(sph)=%.5f  xi0(cyl)=%.5f  xi0(pl)=%.5f" %
              (g, sedov_xi0(g, 3), sedov_xi0(g, 2), sedov_xi0(g, 1)))
