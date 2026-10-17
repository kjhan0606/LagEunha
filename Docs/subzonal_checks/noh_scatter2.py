import sys; sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks')
import numpy as np
from chain1d import Chain, g
def noh(variant, beta, N=200, tend=0.6, P0=1e-6, amp=0.1, kappa=1.0, cfl=0.25, seed=1, label=''):
    rng = np.random.default_rng(seed)
    x = (np.arange(2*N)+0.5)/N - 1.0 + amp/N*rng.uniform(-1, 1, 2*N)
    v = -np.sign(x); rho = np.ones(2*N); P = np.full(2*N, P0)
    xp = np.concatenate([[2*x[0]-x[1]], x, [2*x[-1]-x[-2]]]); m = 0.5*(xp[2:]-xp[:-2])
    ch = Chain(x, v, rho, P, m, variant=variant, beta=beta, kappa=kappa, bc='ghost')
    ch.run(tend, cfl=cfl)
    V, rhon, Pn, dL, dR = ch.state()
    xs = ch.t/3; sel = (np.abs(ch.x) < 0.7*xs) & (np.abs(ch.x) > 0.15*xs)
    lnK = np.log(ch.K[sel]); Kex = (4/3)/4**g
    # gap pressure noise (sub): Phat at gaps inside sel
    d = np.diff(ch.x); rg = m[:-1]/d; Ph = 0.5*(ch.K[:-1]+ch.K[1:])*rg**g
    gsel = sel[:-1] & sel[1:]
    print(f"{variant:4s} b={beta} k={kappa} N={N} cfl={cfl} amp={amp} {label}: std lnK={lnK.std():.4f} std lnP0={np.log(Pn[sel]).std():.4f} std lnPhat={np.log(Ph[gsel]).std():.4f} std rho/4={rhon[sel].std()/4:.4f} mingap={ch.mingap_hist[:,1].min():.3f}")
noh('cell', 1.0); noh('sub', 'W')
noh('lin', 'W', kappa=0.5); noh('lin', 'W', kappa=0.25); noh('lin', 'W', kappa=0.1)
noh('sub', 2.0); noh('sub', 'W', cfl=0.1)
noh('cell', 1.0, N=400); noh('sub', 'W', N=400)
noh('cell', 1.0, seed=2); noh('sub', 'W', seed=2)
