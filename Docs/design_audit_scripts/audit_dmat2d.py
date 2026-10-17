import sys, numpy as np
sys.path.insert(0,'/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad'); sys.path.insert(0,'/tmp/lageunha-push/Docs/subzonal_checks')
from audit_hg2d import *
from sub2d import Mesh
for n, amp in ((8, 0.25), (12, 0.1)):
    x, h = glass(n, amp, 5); N = len(x); Lb = 1.0; m = np.full(N, h*h); mesh = Mesh(x, Lb); V = cell_volumes(mesh); K = 1.0/(m/V)**GAM
    c0 = np.sqrt(GAM/(m/V).mean()); rng = np.random.default_rng(1)
    for mach in (0.0, 1.0, 10.0):
        v = rng.normal(size=x.shape)*mach*c0/np.sqrt(2) if mach > 0 else np.zeros_like(x)
        F, H, _, ops = hg_eval(x, Lb, m, K, lambda t: v[mesh.idx[t]], mesh=mesh, return_ops=True)
        D = np.zeros((2*N, 2*N))
        for (i, f, r, C, ids, LT, gvec, cf, rhof, hh, W) in ops:
            q = abs(r); b = (GAM+1)/4*hh; s = np.sqrt(b*b*q*q + cf*cf)
            Cp = rhof*hh*W*(2*b*q + s + b*b*q*q/s)
            Lvec = np.zeros(2*N)
            for k, t in enumerate(ids):
                p = mesh.idx[t]; Lvec[2*p:2*p+2] += LT[k]
            D += Cp*np.outer(Lvec, Lvec)
        Mi = np.repeat(1/np.sqrt(m), 2)
        B = Mi[:, None]*D*Mi[None, :]
        lam_row = np.abs(B).sum(1).max(); lam_eig = np.linalg.eigvalsh(B).max()
        print(f"n={n} amp={amp} Mach_rms={mach}: Lambda_D rowsum = {lam_row*h/c0:.2f} c/h, max eig = {lam_eig*h/c0:.2f} c/h  (ratio {lam_row/lam_eig:.2f}); design dt_spec(hg only)=0.25/Lambda = {0.25/lam_row*c0/h:.3f} h/c vs RK4 limit 2.785/eig = {2.785/lam_eig*c0/h:.3f} h/c")
