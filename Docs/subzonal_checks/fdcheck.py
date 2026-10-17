import numpy as np, sys
sys.path.insert(0, '/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
from chain1d import Chain, g
rng = np.random.default_rng(0)
N = 20
x = np.sort(np.arange(N)*0.1 + rng.uniform(-0.03, 0.03, N))
v = rng.normal(0, 0.1, N); rho = np.ones(N); P = 1+0.3*rng.uniform(size=N); m = np.full(N, 0.1)
for var in ['cell', 'sub', 'lin']:
    ch = Chain(x, v, rho, P, m, variant=var, beta=0.0)
    K = ch.K.copy()
    # set Q to zero by making v uniform
    v0 = np.zeros(N)
    _, a, _, _ = ch.rhs(x, v0, K)
    F = a*m
    Ffd = np.zeros(N); h = 1e-6
    for i in range(N):
        xp = x.copy(); xp[i] += h; xm = x.copy(); xm[i] -= h
        Ffd[i] = -(ch.energy(xp, v0, K)-ch.energy(xm, v0, K))/(2*h)
    print(var, "max|F-Ffd| =", np.abs(F-Ffd).max(), " max|F| =", np.abs(F).max())
    # energy identity with Q: dE/dt = sum v.F + sum K' G  should be 0
    _, a, Kd, _ = ch.rhs(x, v, K)
    G, _, _ = ch.G_and_grad(x)
    # reversible part: sum m v a_rev + sum K dG/dt ; total with Q: sum m v a + sum Kd G + sum K dG/dt
    dGdt = np.zeros(N); hh = 1e-7
    G2, _, _ = ch.G_and_grad(x+hh*v); dGdt = (G2-G)/hh
    print("   dE/dt (should be ~0):", np.sum(m*v*a)+np.sum(Kd*G)+np.sum(K*dGdt))
