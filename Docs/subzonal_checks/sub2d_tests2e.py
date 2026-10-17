"""Unit tests (e) flip crossing and (f) adiabatic glass dynamics on the 2D reference."""
import numpy as np
from sub2d import total, U0, Mesh, G

rng = np.random.default_rng(5)
n = 8; Lb = 1.0; h = Lb/n
i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
lat = np.stack([(i.ravel()+0.5)*h, (j.ravel()+0.5)*h], 1)
N = len(lat); m = np.full(N, h*h); K = np.ones(N)

# (e) move generator a along y across the square-lattice flip (x_a.x offset breaks 4-point degeneracy)
x0 = (lat + rng.uniform(-0.02, 0.02, lat.shape)*h) % Lb
a = 27
def nfaces(x):
    M = Mesh(x, Lb)   # topology signature: neighbour sets of all cells
    return hash(tuple(tuple(sorted(int(M.idx[q]) for q in nb)) for nb in M.nbr))
ss = np.linspace(-0.2, 0.2, 161)*h
U = []; F = []; nf = []
for s_ in ss:
    x = x0.copy(); x[a, 0] += s_; x[a, 1] += 0.3*s_
    u, f, _, _ = total(x, Lb, m, K, grad=True); U.append(u); F.append(f[a, 0] + 0.3*f[a, 1]); nf.append(nfaces(x))
U = np.array(U); F = np.array(F); nf = np.array(nf)
dUds = np.gradient(U, ss)
flips = np.where(np.diff(nf) != 0)[0]
print('(e) flips along path at s/h =', np.round(ss[flips]/h, 3))
print('    max |F + dU/ds| / max|F| (smooth parts) = %.1e' % (np.abs(F + dUds)[2:-2].max()/np.abs(F).max()))
for k in flips:
    # bracket each flip with a fine bisection on the face count, then one-sided forces
    lo, hi = ss[k], ss[k+1]
    for _ in range(40):
        mid = 0.5*(lo+hi); x = x0.copy(); x[a, 0] += mid; x[a, 1] += 0.3*mid
        if nfaces(x) == nf[k]: lo = mid
        else: hi = mid
    out = []
    for e in (1e-7, 1e-9):
        for side in (-1, 1):
            sv = (lo if side < 0 else hi) + side*e; x = x0.copy(); x[a, 0] += sv; x[a, 1] += 0.3*sv
            u, f, _, _ = total(x, Lb, m, K, grad=True); out.append((u, f[a, 0] + 0.3*f[a, 1]))
    (uL, fL), (uR, fR) = out[2], out[3]
    print('    flip at s/h=%.6f: U jump %.1e (|U|=%.3e), force left %.6e right %.6e jump/|F| = %.2e'
          % (lo/h, uR-uL, abs(uL), fL, fR, abs(fR-fL)/max(abs(fL), 1e-300)))

