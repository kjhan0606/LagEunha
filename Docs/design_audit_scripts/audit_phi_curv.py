import sys, numpy as np
sys.path.insert(0,'/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad'); sys.path.insert(0,'/tmp/lageunha-push/Docs/subzonal_checks')
from audit_hg2d import glass, cell_volumes, GAM
from sub2d import Mesh, total, U0
for n in (16, 24, 32):
    x, h = glass(n, 0.2, 5); N = len(x); Lb = 1.0; m = np.full(N, h*h); mesh = Mesh(x, Lb); V = cell_volumes(mesh)
    r = np.hypot(x[:,0]-0.5, x[:,1]-0.5); vphi = np.where(r<0.2, 5*r, np.where(r<0.4, 2-5*r, 0.0))
    P = np.where(r<0.2, 5+12.5*r**2, np.where(r<0.4, 9+12.5*r**2-20*r+4*np.log(r/0.2), 3+4*np.log(2))); K = P/(m/V)**GAM
    v = np.stack([-vphi*(x[:,1]-0.5)/np.maximum(r,1e-12), vphi*(x[:,0]-0.5)/np.maximum(r,1e-12)], 1)
    KE = 0.5*np.sum(m*(v**2).sum(1))
    Phi = lambda xx: total(xx % Lb, Lb, m, K)[0] - U0(xx % Lb, Lb, m, K)
    out = []
    for s in (0.0, 0.01, 0.02, 0.04, 0.08):
        out.append(Phi(x + s*v))   # s = time with frozen velocity
    out = np.array(out); s = np.array([0.0, 0.01, 0.02, 0.04, 0.08])
    # quadratic fit Phi(s) = Phi0 + a s + b s^2  -> curvature 2b
    A = np.column_stack([np.ones(5), s, s*s]); coef = np.linalg.lstsq(A, out, rcond=None)[0]
    print(f"n={n}: KE={KE:.3e} Phi(s)-Phi0 = {np.round(out-out[0], 5)} at s={s}; fit a={coef[1]:+.3e} 2b={2*coef[2]:.3e} -> t*=sqrt(2KE/2b)={np.sqrt(2*KE/(2*coef[2])) if coef[2]>0 else float('nan'):.3f}; Phi(0.08)/KE = {(out[-1]-out[0])/KE:.3f}")
