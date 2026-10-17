"""Astra design 8 operator tests: affine/quadratic reproduction, cubic order, ledger, pair absorption eta_f and lambda_pair."""
import sys; sys.path.insert(0, '/tmp/lageunha-push/Docs/design_audit_scripts'); sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks')
import numpy as np
import audit_dyn2d as a
from audit_hg2d import glass, face_xi_grad
from sub2d import Mesh
def build(n, seed=7, amp=0.25):
    x0, h = glass(n, amp, seed); N = len(x0); mesh = Mesh(x0, 1.0); cen = 4*N
    allf = {i: face_xi_grad(mesh, i, x0, cen)[0] for i in range(N)}
    return x0, h, N, mesh, cen, allf
x0, h, N, mesh, cen, allf = build(8)
c0 = np.array([0.5, 0.5])
fields = dict(affine=lambda p: np.array([0.3 + 0.37*p[0] - 1.21*p[1], -0.7 + 0.83*p[0] + 0.52*p[1]]),
              quadratic=lambda p: np.array([0.2 + p[0]**2 - 0.7*p[0]*p[1] + 0.4*p[1]**2, -0.1 + 0.6*p[0]**2 + 1.3*p[0]*p[1] - 0.8*p[1]**2]))
for mode in ('lin', 'quad', 'quadx'):
    a.FIT_STATS.update(kmax=0.0, fail=0)
    ops = a.bprime_ops_q(mesh, x0, allf, cen, mode)
    out = []
    for lab, fn in fields.items():
        out.append('%s %.1e' % (lab, max(abs(sum(c @ fn(p) for _, p, c in t)) for (_, _, _, _, t) in ops)))
    # monomial reproduction: each of 6 monomials in each component
    mono = [lambda p: 1.0, lambda p: p[0], lambda p: p[1], lambda p: p[0]**2, lambda p: p[0]*p[1], lambda p: p[1]**2]
    mx = 0.0
    for mfn in mono:
        for comp in (0, 1):
            fn = (lambda p, mfn=mfn, comp=comp: np.array([mfn(p), 0.0]) if comp == 0 else np.array([0.0, mfn(p)]))
            mx = max(mx, max(abs(sum(c @ fn(p) for _, p, c in t)) for (_, _, _, _, t) in ops))
    m = np.full(N, h*h)
    eta = []
    for (i, j, tj, f, t) in ops:
        lb = sum((c @ (np.array([1.0 if q == i else (-1.0 if q == j else 0.0)]*1)[0]*(t[0][2]))) / m[q] for q, _, c in t if q in (i, j))
        ab = 1/m[i] + 1/m[j]
        # l^T M^-1 b with b = +e on i, -e on j (merged over images of i, j)
        e = t[0][2]; lMb = sum((c @ e)/m[q]*(1 if q == i else -1) for q, _, c in t if q in (i, j))
        eta.append(lMb/ab)
    eta = np.array(eta)
    print('%-5s affine/quad: %s | 12 monomials max|delta| %.1e | cond max %.1e fails %d | eta_f min/median %.3f/%.3f'
          % (mode, ', '.join(out), mx, a.FIT_STATS['kmax'], a.FIT_STATS['fail'], eta.min(), np.median(eta)))
# cubic order: smooth periodic field, residual vs resolution
for mode in ('lin', 'quad', 'quadx'):
    rr = []
    for n in (8, 12, 16, 24):
        x0, h, N, mesh, cen, allf = build(n)
        ops = a.bprime_ops_q(mesh, x0, allf, cen, mode)
        k = 2*np.pi
        fn = lambda p: np.array([np.sin(k*p[1]) + 0.5*np.cos(k*(p[0] + p[1])), np.cos(k*p[0])*np.sin(k*p[1])])
        rr.append(np.sqrt(np.mean([sum(c @ fn(p) for _, p, c in t)**2 for (_, _, _, _, t) in ops])))
    rr = np.array(rr); hs = 1/np.array([8, 12, 16, 24])
    print('%-5s smooth periodic rms delta: %s  slopes %s' % (mode, ' '.join('%.2e' % r for r in rr), ' '.join('%.2f' % s for s in np.diff(np.log(rr))/np.diff(np.log(hs)))))
