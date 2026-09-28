"""Noh (1987) implosion, exact solution, planar/cyl/spherical (j=1,2,3).

Cold gas (p=0), rho0, radial velocity -v0 everywhere at t=0. A stagnation
shock moves out at D = (gamma-1) v0 / 2.
  r < D t : rho = rho0 ((gamma+1)/(gamma-1))^j, v = 0,
            p = (gamma-1) rho e with e = v0^2/2  (p = rho0 v0^2 ((g+1)/(g-1))^(j-1) (g+1)/2 ... ),
  r > D t : rho = rho0 (1 + v0 t / r)^(j-1), v = -v0, p = 0.
For gamma=5/3, j=3: post-shock rho = 64, p = 64/3 = 21.33, D = 1/3.
"""
import numpy as np


def noh_profile(r, t, gamma=5.0 / 3.0, rho0=1.0, v0=1.0, j=3, p0=0.0):
    r = np.asarray(r, dtype=float)
    D = 0.5 * (gamma - 1.0) * v0
    rs = D * t
    comp = (gamma + 1.0) / (gamma - 1.0)
    rho = np.where(r < rs, rho0 * comp ** j,
                   rho0 * (1.0 + v0 * t / np.maximum(r, 1e-300)) ** (j - 1))
    v = np.where(r < rs, 0.0, -v0)
    p_post = (gamma - 1.0) * rho0 * comp ** j * 0.5 * v0 ** 2
    p = np.where(r < rs, p_post, p0)
    return rho, v, p, rs


if __name__ == "__main__":
    rho, v, p, rs = noh_profile(np.array([0.1, 1.0]), 2.0)
    print("t=2: rs=%.6f rho_post=%.6f p_post=%.6f rho(r=1)=%.6f" %
          (rs, rho[0], p[0], rho[1]))
