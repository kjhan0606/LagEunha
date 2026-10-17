import numpy as np
from sub2d import Mesh
from hg2d_ref import face_grads, xi_all
rng = np.random.default_rng(11); n = 8; Lb = 1.0; h = Lb/n
i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
x = (np.stack([(i.ravel()+.5)*h, (j.ravel()+.5)*h], 1) + rng.uniform(-.25, .25, (n*n, 2))*h) % Lb
N = len(x); mesh = Mesh(x, Lb)
# 1) gradient of xi vs finite difference of xi
v = rng.normal(size=x.shape); eps = 1e-6
xp = xi_all((x + eps*v) % Lb, Lb); xm = xi_all((x - eps*v) % Lb, Lb)
err = 0; big = 0
for ii in range(N):
    V, faces = face_grads(mesh, ii)
    for (jj, A, d, u, g) in faces:
        key = (ii, int(mesh.idx[jj]))
        if key not in xp or key not in xm: continue
        fd = (xp[key] - xm[key])/(2*eps)
        an = sum(gb @ v[int(mesh.idx[bb])] for bb, gb in g.items())
        err = max(err, abs(fd - an)); big = max(big, abs(an))
print('xi_dot analytic vs FD: max err %.2e (max |xi_dot| %.2e)' % (err, big))
# 2) translation invariance: sum_b g_b = 0 (by image identity)
mx = 0
for ii in range(N):
    V, faces = face_grads(mesh, ii)
    for (jj, A, d, u, g) in faces:
        mx = max(mx, np.abs(sum(g.values())).max())
print('sum_b g_b (translation): %.2e' % mx)
