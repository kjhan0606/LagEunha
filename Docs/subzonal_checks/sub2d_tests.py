"""Unit tests (b)-(g) of Docs/subzonal_derivation.md section 7 on the 2D reference."""
import numpy as np
from sub2d import total, U0, G

rng = np.random.default_rng(3)
n = 8; Lb = 1.0; h = Lb/n
i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
lat = np.stack([(i.ravel()+0.5)*h, (j.ravel()+0.5)*h], 1)
N = len(lat); m = np.full(N, h*h); K = np.ones(N)

def phi(x, m, K):
    U, F, per, mesh = total(x, Lb, m, K, grad=True)
    return U - U0(x, Lb, m, K), U, F

# (b) lattice and uniformly sheared lattice
p0, _, F0 = phi(lat, m, K)
print('(b) square lattice: Phi=%.2e max|F|=%.2e' % (p0, np.abs(F0).max()))
# periodic-consistent shear: shift row j by s*j*h with s*n integer
for s in (1, 3):
    x = lat.copy(); x[:, 0] = (lat[:, 0] + s*h*(lat[:, 1] - 0.5*h)/(n*h)*n/n) % Lb
    pS, _, FS = phi(x, m, K)
    print('(b) sheared lattice (row shift %d/%d): Phi=%.2e max|F|=%.2e' % (s, n, pS, np.abs(FS).max()))

# (c) perturbed lattice: Phi_i > 0, analytic force vs FD
x = (lat + rng.uniform(-0.25, 0.25, lat.shape)*h) % Lb
mR = m*rng.uniform(0.5, 2.0, N); KR = rng.uniform(0.5, 2.0, N)    # (g) non-uniform mass and entropy
for lab, mm, KK in (('uniform m,K', m, K), ('random m,K', mR, KR)):
    U, F, per, mesh = total(x, Lb, mm, KK, grad=True)
    Phi = U - U0(x, Lb, mm, KK)
    ph_i = []
    for ii, (Ui, d) in enumerate(per):
        ph_i.append(Ui - mm[ii]*KK[ii]*(mm[ii]/d['V'])**(G-1)/(G-1))
    eps = 1e-6; errs = []
    for a in rng.choice(N, 6, replace=False):
        for comp in (0, 1):
            xp = x.copy(); xp[a, comp] += eps; xm = x.copy(); xm[a, comp] -= eps
            fd = -(total(xp, Lb, mm, KK)[0] - total(xm, Lb, mm, KK)[0])/(2*eps)
            errs.append(abs(fd - F[a, comp])/(np.abs(F).max()))
    print('(c,g) %s: Phi=%.3e (Phi/U=%.2e) min Phi_i=%.2e  |F_an - F_fd|/max|F| max=%.1e  (d) |sum F|/sum|F|=%.1e'
          % (lab, Phi, Phi/U, min(ph_i), max(errs), np.hypot(*F.sum(0))/np.hypot(F[:, 0], F[:, 1]).sum()))
