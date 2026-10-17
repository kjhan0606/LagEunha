import numpy as np, sys
sys.path.insert(0, '/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
from chain1d import Chain, g
from exact_riemann import exact_riemann

def pairing(x, sel=None):
    d = np.diff(x)
    r = np.abs(np.log(d[1:]/d[:-1]))   # for cell i=1..N-2
    if sel is None: return r.max()
    return r[sel[1:-1]].max() if np.any(sel[1:-1]) else np.nan

def sod(variant, beta, N=200, tend=0.15, kappa=1.0, cfl=0.25, perturb=0.0):
    # equal-mass cells, rho_L=1 on [0,0.5], rho_R=0.125 on [0.5,1.5]  (domain long enough on right)
    mL = 0.5/ (N*0.8); nL = int(0.5/mL); nR = int(1.0*0.125/mL)
    xL = 0.5*(np.arange(nL)+0.5)/nL
    xR = 0.5 + 1.0*(np.arange(nR)+0.5)/nR
    x = np.concatenate([xL, xR]); n = len(x)
    rng = np.random.default_rng(3)
    if perturb > 0:
        dx = np.diff(x).min()
        x = x + perturb*dx*rng.uniform(-1, 1, n)
    rho = np.where(x < 0.5, 1.0, 0.125); P = np.where(x < 0.5, 1.0, 0.1)
    m = np.full(n, mL); v = np.zeros(n)
    ch = Chain(x, v, rho, P, m, variant=variant, beta=beta, kappa=kappa, bc='fixed_ends')
    ch.run(tend, cfl=cfl)
    V, rhon, Pn, dL, dR = ch.state()
    re, ue, pe, ps, us = exact_riemann(1, 0, 1, 0.125, 0, 0.1, g, ch.x, ch.t, 0.5)
    # shock location exact
    cR = np.sqrt(g*0.1/0.125); SR = cR*np.sqrt((g+1)/(2*g)*ps/0.1+(g-1)/(2*g))
    xs = 0.5+SR*ch.t
    # post-shock plateau region: between contact and shock, use inner 60%
    xc = 0.5+us*ch.t
    sel = (ch.x > xc+0.25*(xs-xc)) & (ch.x < xs-0.25*(xs-xc))
    rs = re[sel][0]
    ring = np.std(rhon[sel])/rs; bias = (np.mean(rhon[sel])-rs)/rs
    # left plateau (between rarefaction tail and contact)
    cL = 1.0*np.sqrt(g); ST = us-np.sqrt(g*ps/ (1.0*(ps/1.0)**(1/g)))
    sel2 = (ch.x > ST*ch.t+0.5+0.25*(xc-ST*ch.t-0.5)) & (ch.x < xc-0.25*(xc-ST*ch.t-0.5))
    ring2 = np.std(rhon[sel2])/re[sel2][0]
    w = (ch.x>0.3)&(ch.x<1.1); L1 = np.sum((np.abs(rhon-re)*V)[w])/np.sum(V[w])
    # numerical shock position: last cell with rho > (rs+0.125)/2
    xs_num = ch.x[np.where(rhon > 0.5*(rs+0.125))[0].max()]
    # mingap relative to local mean spacing, over history
    mg = ch.mingap_hist[:, 1].min()
    print(f"SOD {variant:4s} beta={beta}: steps={ch.nstep} dE={ch.dE:.1e} L1rho={L1:.4f} "
          f"post-shock bias={bias:+.4f} ring={ring:.4f} ringL={ring2:.4f} xs_err={xs_num-xs:+.4f} "
          f"pair_post={pairing(ch.x,sel):.3f} pair_left={pairing(ch.x,sel2):.3f} mingap/local={mg:.3f} dtmin={ch.dtmin:.1e}")
    return ch, (re, ue, pe)

def noh(variant, beta, N=200, tend=0.6, P0=1e-6, kappa=1.0, cfl=0.25):
    # symmetric: gas on [-1,1], u = -sign(x), rho=1, cold. Use half domain with mirror by full symmetric domain.
    x = (np.arange(2*N)+0.5)/N - 1.0
    v = -np.sign(x); rho = np.ones(2*N); P = np.full(2*N, P0); m = np.full(2*N, 1.0/N)
    ch = Chain(x, v, rho, P, m, variant=variant, beta=beta, kappa=kappa, bc='ghost')
    ch.run(tend, cfl=cfl)
    V, rhon, Pn, dL, dR = ch.state()
    xs = ch.t/3
    sel = (np.abs(ch.x) < 0.6*xs) & (np.abs(ch.x) > 0.15*xs)
    bias = np.mean(rhon[sel])/4-1; ring = np.std(rhon[sel])/4
    # shock width: cells with 1.3<rho<3.7 on x>0
    wid = np.sum((rhon > 1.3) & (rhon < 3.7) & (ch.x > 0))
    Pbias = np.mean(Pn[sel])/(4/3)-1
    xs_num = ch.x[ch.x > 0][np.where(rhon[ch.x > 0] > 2.5)[0].max()]
    mg = ch.mingap_hist[:, 1].min()
    # wall heating: density at innermost cells
    inner = rhon[np.argsort(np.abs(ch.x))[:4]]
    print(f"NOH {variant:4s} N={N} P0={P0:.0e} beta={beta}: width={wid} steps={ch.nstep} dE={ch.dE:.1e} rho bias={bias:+.4f} ring={ring:.4f} "
          f"P bias={Pbias:+.4f} xs_err={xs_num-xs:+.4f} rho_wall={inner.mean():.3f} "
          f"pair_post={pairing(ch.x,sel):.3f} mingap/local={mg:.3f} dtmin={ch.dtmin:.1e}")
    return ch

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    res = {}
    if which in ('sod', 'all'):
        for var in ['cell', 'sub', 'lin']:
            for beta in [0.0, 1.0]:
                res[('sod', var, beta)] = sod(var, beta)
    if which in ('noh', 'all'):
        for var in ['cell', 'sub', 'lin']:
            for beta in [0.0, 1.0]:
                res[('noh', var, beta)] = noh(var, beta)
    np.save('/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad/res_%s.npy' % which,
            {k: (v[0].x, v[0].v, v[0].K, v[0].m) if isinstance(v, tuple) else (v.x, v.v, v.K, v.m) for k, v in res.items()}, allow_pickle=True)
