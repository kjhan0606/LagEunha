"""Gate 1 for (b'): affine fields give delta_f = 0; ledger, total force and torque for random velocity."""
import sys; sys.path.insert(0, '/tmp/lageunha-push/Docs/design_audit_scripts'); sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks')
import numpy as np
import audit_dyn2d as a
from audit_hg2d import glass, cell_volumes, face_xi_grad
from sub2d import Mesh
x0, h = glass(8, 0.25, 7); N = len(x0); mesh = Mesh(x0, 1.0); cen = 4*N
allf = {i: face_xi_grad(mesh, i, x0, cen)[0] for i in range(N)}
sl = a.slope_ops(mesh, x0, allf, cen); ops = a.bprime_ops(mesh, x0, allf, sl, cen)
for lab, A in (('translation', np.zeros((2, 2))), ('rotation', np.array([[0, -1], [1, 0.]])), ('shear', np.array([[0, 1], [0, 0.]])),
               ('uniform compression', -np.eye(2)), ('general', np.array([[0.3, -1.2], [0.8, 0.5]]))):
    b0 = np.array([0.4, -0.9])
    mx = max(abs(sum(c @ (b0 + A @ p) for _, p, c in t)) for (_, _, _, _, t) in ops)
    print('%-20s max |delta_f| = %.2e' % (lab, mx))
rng = np.random.default_rng(2); v = rng.normal(size=x0.shape); m = np.full(N, h*h); C = rng.uniform(0.5, 2, len(ops))
F = np.zeros_like(x0); Hs = 0.0
for Cf, (_, _, _, _, t) in zip(C, ops):
    dl = sum(c @ v[p] for p, _, c in t)
    for p, _, c in t: F[p] -= Cf*dl*c
    Hs += Cf*dl*dl
print('faces %d  ledger (sum v.F + H)/H = %.1e  |sum F|/sum|F| = %.1e' % (len(ops), (np.sum(v*F) + Hs)/Hs, np.linalg.norm(F.sum(0))/np.abs(F).sum()))
tq = 0.0
for Cf, (_, _, _, _, t) in zip(C, ops):
    dl = sum(c @ v[p] for p, _, c in t)
    for p, pos, c in t: tq += pos[0]*(-Cf*dl*c[1]) - pos[1]*(-Cf*dl*c[0])
print('torque (unwrapped positions) / sum|F| = %.1e' % (abs(tq)/np.abs(F).sum()))
