"""Unit tests (e) flip crossing and (f) adiabatic glass dynamics on the 2D reference."""
import numpy as np
from sub2d import total, U0, Mesh, G

rng = np.random.default_rng(5)
n = 8; Lb = 1.0; h = Lb/n
i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
lat = np.stack([(i.ravel()+0.5)*h, (j.ravel()+0.5)*h], 1)
N = len(lat); m = np.full(N, h*h); K = np.ones(N)

# (e) move generator a along y across the square-lattice flip (x_a.x offset breaks 4-point degeneracy)
x0 = (lat + rng.uniform(-0.15, 0.15, lat.shape)*h) % Lb
a = 27
def nfaces(x):
    M = Mesh(x, Lb); return sum(len(c) for c in M.cells)
ss = np.linspace(-0.45, 0.45, 181)*h
U = []; F = []; nf = []
for s_ in ss:
    x = x0.copy(); x[a, 1] += s_
    u, f, _, _ = total(x, Lb, m, K, grad=True); U.append(u); F.append(f[a, 1]); nf.append(nfaces(x))
U = np.array(U); F = np.array(F); nf = np.array(nf)
dUds = np.gradient(U, ss)
flips = np.where(np.diff(nf) != 0)[0]
print('(e) flips along path at s/h =', np.round(ss[flips]/h, 3))
print('    max |F + dU/ds| / max|F| (smooth parts) = %.1e' % (np.abs(F + dUds)[2:-2].max()/np.abs(F).max()))
for k in flips:
    # bracket each flip with a fine bisection on the face count, then one-sided forces
    lo, hi = ss[k], ss[k+1]
    for _ in range(40):
        mid = 0.5*(lo+hi); x = x0.copy(); x[a, 1] += mid
        if nfaces(x) == nf[k]: lo = mid
        else: hi = mid
    out = []
    for e in (1e-7, 1e-9):
        for side in (-1, 1):
            x = x0.copy(); x[a, 1] += (lo if side < 0 else hi) + side*e
            u, f, _, _ = total(x, Lb, m, K, grad=True); out.append((u, f[a, 1]))
    (uL, fL), (uR, fR) = out[2], out[3]
    print('    flip at s/h=%.6f: U jump %.1e (|U|=%.3e), force left %.6e right %.6e jump/|F| = %.2e'
          % (lo/h, uR-uL, abs(uL), fL, fR, abs(fR-fL)/max(abs(fL), 1e-300)))

# (f) adiabatic glass: uniform P, generators off-centre, no dissipation; RK4 energy conservation
x0 = (lat + rng.uniform(-0.2, 0.2, lat.shape)*h) % Lb
mesh = Mesh(x0, Lb)
V0 = np.array([0.5*np.sum(np.cross(np.array([v for v, _, _ in c]), np.roll(np.array([v for v, _, _ in c]), -1, 0))) for c in mesh.cells])
Kg = 1.0/(m/V0)**G          # P = 1 in every cell at t=0
c0 = np.sqrt(G*1.0/(m/V0).mean())
def rhs(x, v):
    _, f, _, _ = total(x % Lb, Lb, m, Kg, grad=True); return v, f/m[:, None]
def energy(x, v):
    return total(x % Lb, Lb, m, Kg)[0] + 0.5*np.sum(m[:, None]*v*v)
def run(dt, T):
    x = x0.copy(); v = np.zeros_like(x); E0 = energy(x, v); ts = []; ph = []
    for k in range(int(round(T/dt))):
        k1 = rhs(x, v); k2 = rhs(x + 0.5*dt*k1[0], v + 0.5*dt*k1[1])
        k3 = rhs(x + 0.5*dt*k2[0], v + 0.5*dt*k2[1]); k4 = rhs(x + dt*k3[0], v + dt*k3[1])
        x = x + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); v = v + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        ts.append((k+1)*dt); ph.append(total(x % Lb, Lb, m, Kg)[0] - U0(x % Lb, Lb, m, Kg))
    return (energy(x, v) - E0)/abs(E0), np.array(ts), np.array(ph), x
T = 2.0*h/c0
print('(f) c0=%.3f h/c0=%.3e  Phi(0)=%.3e' % (c0, h/c0, total(x0, Lb, m, Kg)[0] - U0(x0, Lb, m, Kg)))
prev = None
for dt in (h/c0/10, h/c0/20, h/c0/40):
    dE, ts, ph, xe = run(dt, T)
    print('    dt=%.2e: |dE/E|=%.2e  Phi(T)=%.3e  min/max Phi over run %.2e/%.2e' % (dt, abs(dE), ph[-1], ph.min(), ph.max()))
    if prev is not None: print('      position change vs previous dt: %.2e' % np.abs(((xe - prev + 0.5*Lb) % Lb) - 0.5*Lb).max())
    prev = xe
