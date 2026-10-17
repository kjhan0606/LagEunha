"""1D Noh with disordered initial spacing: post-shock std(ln K) for cell vs subzonal."""
import sys; sys.path.insert(0, '/tmp/lageunha-push/Docs/subzonal_checks')
import numpy as np
from chain1d import Chain, g
def noh(variant, beta, N=200, tend=0.6, P0=1e-6, amp=0.0, massmode='equal', seed=1):
    rng = np.random.default_rng(seed)
    x = (np.arange(2*N)+0.5)/N - 1.0
    x = x + amp/N*rng.uniform(-1, 1, 2*N)
    v = -np.sign(x); rho = np.ones(2*N); P = np.full(2*N, P0)
    if massmode == 'equal': m = np.full(2*N, 1.0/N)
    else:
        xp = np.concatenate([[2*x[0]-x[1]], x, [2*x[-1]-x[-2]]]); m = 0.5*(xp[2:]-xp[:-2])   # rho=1 exactly
    ch = Chain(x, v, rho, P, m, variant=variant, beta=beta, bc='ghost')
    ch.run(tend, cfl=0.25)
    V, rhon, Pn, dL, dR = ch.state()
    xs = ch.t/3
    sel = (np.abs(ch.x) < 0.7*xs) & (np.abs(ch.x) > 0.15*xs)
    lnK = np.log(ch.K[sel]); Kex = (4/3)/4**g
    print(f"{variant:4s} beta={beta} amp={amp} m={massmode}: steps={ch.nstep} dE={ch.dE:.1e} std lnK={lnK.std():.4f} mean lnK-exact={lnK.mean()-np.log(Kex):+.4f} std rho/4={rhon[sel].std()/4:.4f} mingap={ch.mingap_hist[:,1].min():.3f}")
for amp in (0.0, 0.1, 0.25):
    for mm in ('equal', 'rho1'):
        noh('cell', 1.0, amp=amp, massmode=mm)
        noh('sub', 'W', amp=amp, massmode=mm)
        noh('sub', 1.0, amp=amp, massmode=mm)
