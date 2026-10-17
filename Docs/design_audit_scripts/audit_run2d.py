import sys, time, numpy as np
sys.path.insert(0, '/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks')
from audit_hg2d import *
from sub2d import total, U0, Mesh


def formula_checks(n=8, amp=0.25, seed=5):
    x, h = glass(n, amp, seed); N = len(x); Lb = 1.0
    m = np.full(N, h*h); mesh = Mesh(x, Lb); V = cell_volumes(mesh); K = 1.0/(m/V)**GAM   # P = 1
    c0 = np.sqrt(GAM*1.0/(m/V).mean()); cen = 4*N
    rng = np.random.default_rng(1)
    print(f"[2D glass n={n} +-{amp}h, N={N}, P=1, c0={c0:.3f}]")
    # (1) xi = ln(h_c/h_x) == ln(B/W); analytic g vs FD
    xi_hw = []; xi_bw = []
    for i in range(N):
        faces, gm = face_xi_grad(mesh, i, x, cen)
        for f in faces:
            xi_hw.append(f['xi']); xi_bw.append(np.log(f['B']/f['W']))
    xi_hw = np.array(xi_hw); xi_bw = np.array(xi_bw)
    print(f"  xi: ln(h_c/h_x) vs ln(B/W): max|diff| = {np.abs(xi_hw-xi_bw).max():.1e}  (rms xi = {np.sqrt((xi_hw**2).mean()):.3f}, max|xi| = {np.abs(xi_hw).max():.3f})")
    w = rng.normal(size=x.shape)
    def xis(xx):
        ms = Mesh(xx % Lb, Lb); out = {}
        for i in range(N):
            faces, gm = face_xi_grad(ms, i, xx % Lb, cen)
            for f in faces:
                out[(i, ms.idx[f['j']])] = f['xi']
        return out
    eps = 1e-6
    xp, xm = xis(x + eps*w), xis(x - eps*w)
    an = {}
    for i in range(N):
        faces, gm = face_xi_grad(mesh, i, x, cen)
        for f in faces:
            an[(i, mesh.idx[f['j']])] = sum(gg @ w[mesh.idx[t]] for t, gg in f['g'].items())
    keys = sorted(an)
    fd = np.array([(xp[k]-xm[k])/(2*eps) for k in keys]); aa = np.array([an[k] for k in keys])
    print(f"  analytic g_if (h-ratio formula, centroid chain + normal rotation) vs FD xidot: max|diff| = {np.abs(fd-aa).max():.1e} (max|xidot| = {np.abs(aa).max():.2f})")
    # translation invariance of g: sum_a g_a = 0 ; rotation: sum_a g_a . (z x (X_a - x_i)) = 0 ? (not required: xi is NOT affine-invariant on disordered mesh)
    s0 = max(np.abs(sum(gg for gg in f['g'].values())).max() for i in range(N) for f in face_xi_grad(mesh, i, x, cen)[0])
    print(f"  sum_a g_if,a (translation) max = {s0:.1e}")
    # (2) affine fields, unwrapped: r_if with filter
    fields = {'translation': lambda X: np.array([0.7, -0.2]),
              'rotation': lambda X: np.array([-(X[1]-0.5), X[0]-0.5]),
              'uniform compression': lambda X: -(X - 0.5),
              'uniform shear v_x=y': lambda X: np.array([X[1]-0.5, 0.0]),
              'general affine': lambda X: np.array([0.3, -0.1]) + np.array([[0.5, 1.2], [-0.7, 0.2]]) @ (X - 0.5)}
    for name, f in fields.items():
        for filt in (True, False):
            F, H, _, ops = hg_eval(x, Lb, m, K, lambda t: f(mesh.X[t]), mesh=mesh, filt=filt, return_ops=True)
            rr = np.array([o[2] for o in ops])
            ref = c0*h*V.sum()   # rho c h sum V |grad v|^2 with rho=1,|grad v|~1
            if filt:
                print(f"  {name:22s}: filtered  max|r| = {np.abs(rr).max():.1e}  sum H = {H.sum():.1e}  max|F| = {np.abs(F).max():.1e}", end='')
            else:
                print(f"   | unfiltered rms r = {np.sqrt((rr**2).mean()):.3e}, H/(rho c h sumV) = {H.sum()/ref:.2e}")
    # (3) random field: identities
    v = rng.normal(size=x.shape)*0.3*c0
    F, H, _, ops = hg_eval(x, Lb, m, K, lambda t: v[mesh.idx[t]], mesh=mesh, return_ops=True)
    torque = 0.0
    for (i, f, r, C, ids, LT, gvec, cf, rhof, hh, W) in ops:
        Xs = np.array([mesh.X[t] for t in ids]) - x[i]
        torque += -C*r*np.sum(Xs[:, 0]*LT[:, 1] - Xs[:, 1]*LT[:, 0])
    print(f"  random v (0.3c): sum F = {np.abs(F.sum(0)).max():.1e} (mean|F| {np.abs(F).mean():.2e}); torque (local frames) = {torque:.1e}; "
          f"sum v.F + sum H = {(v*F).sum()+H.sum():.1e} (sum H = {H.sum():.3e}, min H_i = {H.min():.1e})")
    # (4) F = -dR/dv by FD (Rayleigh function)
    Rt = rayleigh(ops)
    eps = 1e-6; a = 7; errs = []
    for comp in (0, 1):
        vp = v.copy(); vp[a, comp] += eps; vm = v.copy(); vm[a, comp] -= eps
        _, _, _, op_p = hg_eval(x, Lb, m, K, lambda t: vp[mesh.idx[t]], mesh=mesh, return_ops=True)
        _, _, _, op_m = hg_eval(x, Lb, m, K, lambda t: vm[mesh.idx[t]], mesh=mesh, return_ops=True)
        errs.append(-(rayleigh(op_p) - rayleigh(op_m))/(2*eps) - F[a, comp])
    print(f"  F_a vs -dR/dv_a (FD, Hugoniot Rayleigh): |diff| = {np.abs(errs).max():.1e} (|F_a| = {np.linalg.norm(F[a]):.3e}); R = {Rt:.3e}")
    # acoustic branch too
    F2, H2, _, ops2 = hg_eval(x, Lb, m, K, lambda t: v[mesh.idx[t]], mesh=mesh, hug=False, return_ops=True)
    print(f"  acoustic branch: sum v.F + sum H = {(v*F2).sum()+H2.sum():.1e}; H_hug/H_acou at 0.3c random = {H.sum()/H2.sum():.3f}")
    # (5) pair-approach absorption by the filter
    i0 = N//2 + n//2
    faces0, gm0 = face_xi_grad(mesh, i0, x, cen)
    j0 = mesh.idx[faces0[0]['j']]; e = (mesh.X[faces0[0]['j']] - x[i0]); e /= np.linalg.norm(e)
    for name, vf in (("pair approach (i,j) Mach 1", lambda t: (e*0.5*c0 if mesh.idx[t] == i0 else (-e*0.5*c0 if mesh.idx[t] == j0 else np.zeros(2)))),
                     ("single kick i Mach 1", lambda t: (e*c0 if mesh.idx[t] == i0 else np.zeros(2))),
                     ("checkerboard Mach 0.1", lambda t: np.array([0.0, 0.1*c0*(-1.0)**(mesh.idx[t] // n + mesh.idx[t] % n)]))):
        Hf = hg_eval(x, Lb, m, K, vf, mesh=mesh, filt=True)[1]; Hu = hg_eval(x, Lb, m, K, vf, mesh=mesh, filt=False)[1]
        Hfa = hg_eval(x, Lb, m, K, vf, mesh=mesh, filt=True, hug=False)[1]; Hua = hg_eval(x, Lb, m, K, vf, mesh=mesh, filt=False, hug=False)[1]
        print(f"  {name:28s}: sum H filtered/unfiltered = {Hf.sum()/Hu.sum():.3f} (Hugoniot), {Hfa.sum()/Hua.sum():.3f} (acoustic); H_hug/H_acou = {Hf.sum()/Hfa.sum():.2f}")
    return x, m, K, mesh


def gresho(n, amp=0.2, seed=5):
    x, h = glass(n, amp, seed); N = len(x); Lb = 1.0
    m = np.full(N, h*h); mesh = Mesh(x, Lb); V = cell_volumes(mesh)
    # Gresho: P(r) profile, vphi = 5r (r<0.2), 2-5r (0.2<r<0.4), 0
    r = np.hypot(x[:, 0]-0.5, x[:, 1]-0.5)
    vphi = np.where(r < 0.2, 5*r, np.where(r < 0.4, 2-5*r, 0.0))
    P = np.where(r < 0.2, 5 + 12.5*r**2, np.where(r < 0.4, 9 + 12.5*r**2 - 20*r + 4*np.log(r/0.2), 3 + 4*np.log(2)))
    K = P/(m/V)**GAM
    v = np.stack([-vphi*(x[:, 1]-0.5)/np.maximum(r, 1e-12), vphi*(x[:, 0]-0.5)/np.maximum(r, 1e-12)], 1)
    KE = 0.5*np.sum(m*(v**2).sum(1))
    t0 = time.time()
    F, H, _, _ = hg_eval(x, Lb, m, K, lambda t: v[mesh.idx[t]], mesh=mesh, filt=True)
    Fu, Hu, _, _ = hg_eval(x, Lb, m, K, lambda t: v[mesh.idx[t]], mesh=mesh, filt=False)
    Fa, Ha, _, _ = hg_eval(x, Lb, m, K, lambda t: v[mesh.idx[t]], mesh=mesh, filt=True, hug=False)
    # ring breakdown
    ring = [(r < 0.18), (r > 0.22) & (r < 0.38), (r > 0.42)]
    # old centroid Phi: Phi(x), Phidot = -sum v.F_Phi, and curvature v^T K_Phi v (FD along v)
    def Phi(xx):
        return total(xx % Lb, Lb, m, K)[0] - U0(xx % Lb, Lb, m, K)
    Ut, Ft, per, _ = total(x, Lb, m, K, grad=True)
    # cell-volume force
    from cellhg2d import cell_force
    Uc, Fc, _, _ = cell_force(x, Lb, m, K)
    FPhi = Ft - Fc
    eps = 1e-4*h
    Ph0 = Phi(x); Php = Phi(x + eps*v); Phm = Phi(x - eps*v)
    Phidot = (Php - Phm)/(2*eps); Phidd = (Php - 2*Ph0 + Phm)/eps**2
    print(f"GRESHO glass n={n} +-{amp}h: KE={KE:.4e}, U={Uc:.3e}; filtered H_hg = {H.sum():.3e} -> tau_hg = KE/H = {KE/H.sum():.3e} "
          f"(acoustic C: {KE/Ha.sum():.3e}); unfiltered (old sec.10) H = {Hu.sum():.3e} -> tau = {KE/Hu.sum():.3e}; "
          f"rings (rigid/shear/outer) H_filt = {[f'{H[s].sum():.2e}' for s in ring]}, H_unfilt = {[f'{Hu[s].sum():.2e}' for s in ring]}; "
          f"old Phi: Phi0/U = {Ph0/Uc:.2e}, dPhi/dt(frozen v) = {Phidot:+.3e} (= -sum v.F_Phi {-(v*FPhi).sum():+.3e}), d2Phi/dt2 = {Phidd:.3e} -> t*=sqrt(2KE/Phidd) = {np.sqrt(2*KE/Phidd):.3f}; "
          f"sum v.F_hg + sum H = {(v*F).sum()+H.sum():.1e}  [{time.time()-t0:.0f}s]")
    return dict(n=n, KE=KE, H=H.sum(), Hu=Hu.sum(), Ha=Ha.sum(), Phidd=Phidd)


if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'formula':
        formula_checks(8, 0.25)
        formula_checks(12, 0.1, seed=7)
    elif which == 'gresho':
        for n in (16, 24, 32):
            gresho(n)
    elif which == 'gresho48':
        gresho(48)
