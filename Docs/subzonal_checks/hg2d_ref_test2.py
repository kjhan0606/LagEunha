import numpy as np
from sub2d import Mesh
from hg2d_ref import hg_terms
rng = np.random.default_rng(12); n = 8; Lb = 1.0; h = Lb/n
i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
x = (np.stack([(i.ravel()+.5)*h, (j.ravel()+.5)*h], 1) + rng.uniform(-.25, .25, (n*n, 2))*h) % Lb
N = len(x); mesh = Mesh(x, Lb); rho, P, m = 1.0, 1.0, h*h
def run(vel, lab):
    Ftot = np.zeros((N, 2)); Htot = 0; rmax = 0
    for ii in range(N):
        F, H, r = hg_terms(mesh, ii, vel, rho, P, m); Htot += H; rmax = max(rmax, np.abs(r).max())
        for k, f in F.items(): Ftot[k] += f
    print('%-26s max|r|=%.2e  sumH=%.3e  sum v.F + sumH = %.2e  |sum F|=%.2e' % (lab, rmax, Htot, (Ftot*vel).sum() + Htot, np.abs(Ftot.sum(0)).max()))
    return Ftot, Htot
# periodic-consistent affine fields: translation, and shear/compression only as local tests (periodic box breaks
# global affine fields at the seam; use cells away from the seam via a large-scale check below)
run(np.tile([0.3, -0.7], (N, 1)), 'translation')
xc = x - 0.5
run(rng.normal(size=(N, 2)), 'random')
# affine fields evaluated on unwrapped positions: test cells whose 2-hop stencil does not cross the seam
for lab, Gm in (('uniform shear', np.array([[0, 1.0], [0, 0]])), ('rotation', np.array([[0, -1.0], [1, 0]])), ('compression', -np.eye(2)), ('general affine', np.array([[0.3, 1.1], [-0.4, 0.7]]))):
    vel = xc @ Gm.T; rmax = 0; cnt = 0
    for ii in range(N):
        if np.max(np.abs(x[ii] - 0.5)) > 0.2: continue
        F, H, r = hg_terms(mesh, ii, vel, rho, P, m); rmax = max(rmax, np.abs(r).max()); cnt += 1
    print('%-26s interior cells %d: max|r| = %.2e' % (lab, cnt, rmax))
