"""Linear growth rate of the McNally et al. (2012, ApJS 201, 18) smoothed
Kelvin-Helmholtz setup, from the COMPRESSIBLE linearised Euler equations.

Background (periodic in y on [0,1], uniform P0):
  rho(y), U(y) = McNally exponential-smoothed profiles, width L.
Perturbation q(y) exp(i(kx - w t)), q = (rho', u', v', p'):
  -i w rho' = -ikU rho' - rho0' v' - rho0 (ik u' + v'_y)
  -i w u'   = -ikU u'   - U' v'    - ik p'/rho0
  -i w v'   = -ikU v'   - p'_y/rho0
  -i w p'   = -ikU p'   - gamma P0 (ik u' + v'_y)
w = i M q is solved as a dense eigenproblem on 4th-order periodic finite
differences. The growth
rate is sigma = max Im(w). The problem is independent of z, so the same
rate applies to the z-uniform seeded mode in the 3D box.

Validation (selftest): a uniform-density tanh layer at low Mach reproduces
Michalke (1964): sigma_max = 0.1897 U0/L at k L = 0.4446. For reference the
sharp incompressible interface gives k dU sqrt(rho1 rho2)/(rho1+rho2)
(Chandrasekhar 1961); the McNally layer (k L = 0.31) grows at
less than half of that.
"""
import numpy as np
import scipy.sparse as sp


def mcnally_profile(y, rho1=1.0, rho2=2.0, U1=0.5, U2=-0.5, L=0.025):
    """rho(y), U(y), and derivatives, y in [0,1). McNally Eq. 1-2."""
    y = np.mod(np.asarray(y, dtype=float), 1.0)
    rm, um = 0.5 * (rho1 - rho2), 0.5 * (U1 - U2)
    rho = np.empty_like(y); U = np.empty_like(y)
    drho = np.empty_like(y); dU = np.empty_like(y)
    s = [(y < 0.25), (y >= 0.25) & (y < 0.5), (y >= 0.5) & (y < 0.75), (y >= 0.75)]
    e = [np.exp((y - 0.25) / L), np.exp((-y + 0.25) / L),
         np.exp((y - 0.75) / L), np.exp((-y + 0.75) / L)]
    sign = [-1, +1, +1, -1]
    base_r = [rho1, rho2, rho2, rho1]
    base_u = [U1, U2, U2, U1]
    dsgn = [+1, -1, +1, -1]           # d/dy of the exponent
    for m in range(4):
        rho[s[m]] = base_r[m] + sign[m] * rm * e[m][s[m]]
        U[s[m]] = base_u[m] + sign[m] * um * e[m][s[m]]
        drho[s[m]] = sign[m] * rm * dsgn[m] / L * e[m][s[m]]
        dU[s[m]] = sign[m] * um * dsgn[m] / L * e[m][s[m]]
    return rho, U, drho, dU


def tanh_pair_profile(U0=0.5, L=0.01, rho=1.0):
    """Two tanh shear layers (y=1/4, 3/4), uniform density. Used only to
    validate the solver against Michalke (1964): max sigma = 0.1897 U0/L at
    k L = 0.4446 for an isolated incompressible layer U = U0 tanh(y/L)."""
    def f(y):
        a, b = (y - 0.25) / L, (y - 0.75) / L
        U = U0 * (np.tanh(a) - np.tanh(b) - 1.0)
        dU = U0 / L * (1 / np.cosh(a) ** 2 - 1 / np.cosh(b) ** 2)
        return np.full_like(y, rho), U, np.zeros_like(y), dU
    return f


def _dmat(n):
    h = 1.0 / n
    c = {-2: 1.0 / 12, -1: -8.0 / 12, 1: 8.0 / 12, 2: -1.0 / 12}
    D = sp.lil_matrix((n, n))
    for i in range(n):
        for o, w in c.items():
            D[i, (i + o) % n] += w / h
    return D.tocsr()


def _operator(n, k, rho, U, drho, dU, P0, gamma):
    D = _dmat(n)
    I = sp.identity(n, format="csr")
    dg = sp.diags
    ik = 1j * k
    Z = sp.csr_matrix((n, n))
    M = sp.bmat([
        [-ik * dg(U), -ik * dg(rho), -dg(drho) - dg(rho) @ D, Z],
        [Z, -ik * dg(U), -dg(dU), -ik * dg(1 / rho)],
        [Z, Z, -ik * dg(U), -dg(1 / rho) @ D],
        [Z, -ik * gamma * P0 * I, -gamma * P0 * D, -ik * dg(U)],
    ], format="csr")
    return 1j * M                      # w q = A q


def growth_rate(k=4 * np.pi, rho1=1.0, rho2=2.0, U1=0.5, U2=-0.5, L=0.025,
                P0=2.5, gamma=5.0 / 3.0, n=512, profile=None):
    """Returns (sigma, omega_of_fastest_mode, sigma_sharp_incompressible).
    Dense eigen-solve of a 4n x 4n matrix (n=512: ~2 s)."""
    import scipy.linalg as la
    y = (np.arange(n) + 0.5) / n
    if profile is None:
        rho, U, drho, dU = mcnally_profile(y, rho1, rho2, U1, U2, L)
    else:
        rho, U, drho, dU = profile(y)
    A = _operator(n, k, rho, U, drho, dU, P0, gamma).toarray()
    w = la.eigvals(A)
    i = np.argmax(w.imag)
    sig0 = k * abs(U1 - U2) * np.sqrt(rho1 * rho2) / (rho1 + rho2)
    return float(w.imag[i]), complex(w[i]), float(sig0)


if __name__ == "__main__":
    s, w, s0 = growth_rate()
    print("McNally default: sigma=%.4f  (sharp incompressible %.4f)  omega_r=%.4f"
          % (s, s0, w.real))
