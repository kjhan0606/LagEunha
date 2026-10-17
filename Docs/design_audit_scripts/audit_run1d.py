import sys, numpy as np
sys.path.insert(0, '/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
from audit_design1d import Design1D, g, filter_check
from exact_riemann import exact_riemann

SP = '/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad'


def identities():
    rng = np.random.default_rng(0); N = 64; L = 1.0
    x = np.sort((np.arange(N) + 0.5 + rng.uniform(-0.3, 0.3, N))/N*L); v = rng.normal(size=N)
    m = (1 + 0.5*rng.uniform(-1, 1, N))/N; rho = np.ones(N); P = 1 + 0.5*rng.uniform(-1, 1, N)
    ch = Design1D(x, v, rho, P, m, bc='periodic', Lbox=L)
    ch0 = Design1D(x, v, rho, P, m, bc='periodic', Lbox=L, hg=None, qon=False)
    K = ch.K
    dx, dv, dK, a = ch.rhs(x, v, K); dx0, dv0, dK0, a0 = ch0.rhs(x, v, K)
    G = ch.G(x)
    Udot = dK*G + K*(-a['P']*0)      # placeholder; compute U-dot via P Vdot below
    # reversible: sum v.F_rev = sum P_i Vdot_i
    d, dL, dR, V, cL, cR = ch.faces(x)
    ddot = np.roll(v, -1) - v
    Vdot = 0.5*(ddot + np.roll(ddot, 1))
    print("[identities, periodic disordered N=64, random v, nonuniform m,P]")
    print("  sum F_rev = %.1e ; sum v.F_rev - sum P Vdot = %.1e" % ((m*dv0).sum(), (v*m*dv0).sum() - (a0['P']*Vdot).sum()))
    Fh = a['Fh']; H = a['H']
    print("  hg: sum F_hg = %.1e ; sum v.F_hg + sum H = %.1e (sum H = %.3e, min H = %.1e)" % (Fh.sum(), (v*Fh).sum() + H.sum(), H.sum(), H.min()))
    print("  Q : min Q*du = %.1e, sum HQ = %.3e" % ((a['Q']*a['du']).min(), a['HQ'].sum()))
    # total semi-discrete energy rate
    Edot = (v*m*dv).sum() + (dK*G).sum() - (a['P']*Vdot).sum()
    print("  total dE/dt (semi-discrete) = %.1e  [KE+U scale %.2f]" % (Edot, ch.energy(x, v, K)))
    # affine fields: uniform compression / translation -> H_hg = 0
    for name, vv in (("translation", np.full(N, 0.7)), ("uniform compression", -(x - 0.5))):
        chp = Design1D(x, vv, rho, P, m, bc='fixed_ends', hg='hug')
        _, dvp, _, ap = chp.rhs(x, vv, K)
        print("  %-20s: max|F_hg| = %.1e, sum H_hg = %.1e" % (name, np.abs(ap['Fh']).max(), ap['H'].sum()))
    # stationary Nyquist displacement: neutral
    xu = (np.arange(N) + 0.5)/N; xn = xu + 0.2/N*(-1.0)**np.arange(N)
    chn = Design1D(xn, 0*xn, np.ones(N), np.ones(N), np.full(N, 1.0/N), bc='periodic', Lbox=1.0)
    _, dvn, dKn, an = chn.rhs(xn, 0*xn, chn.K)
    print("  stationary Nyquist displacement 0.2D: max|F| = %.1e, max Kdot = %.1e (neutral)" % (np.abs(dvn).max(), np.abs(dKn).max()))
    filter_check()


def sod(N=200, tend=0.15, **kw):
    D = 1.0/N; x = D*(np.arange(int(1.5*N)) + 0.5)
    rho = np.where(x < 0.5, 1.0, 0.125); P = np.where(x < 0.5, 1.0, 0.1); m = rho*D; v = 0*x
    ch = Design1D(x, v, rho, P, m, bc='fixed_ends', **kw); ch.run(tend)
    s = ch.summary()
    d, dL, dR, V, _, _ = ch.faces(ch.x); rhon = ch.m/V; Pn = ch.K*rhon**g
    re, ue, pe, ps, us = exact_riemann(1, 0, 1, 0.125, 0, 0.1, g, ch.x, ch.t, 0.5)
    cR = np.sqrt(g*0.1/0.125); SR = cR*np.sqrt((g+1)/(2*g)*ps/0.1 + (g-1)/(2*g))
    xs = 0.5 + SR*ch.t; xc = 0.5 + us*ch.t
    sel = (ch.x > xc + 0.25*(xs-xc)) & (ch.x < xs - 0.25*(xs-xc))
    rs = re[sel][0]
    w = (ch.x > 0.3) & (ch.x < 1.1)
    L1r = np.sum((np.abs(rhon-re)*V)[w])/np.sum(V[w]); L1u = np.sum((np.abs(ch.v-ue)*V)[w])/np.sum(V[w]); L1p = np.sum((np.abs(Pn-pe)*V)[w])/np.sum(V[w])
    lk = np.log(ch.K[sel]); nyq = 0.5*(lk[1:]-lk[:-1]); sm = 0.5*(lk[1:]+lk[:-1])
    ic = np.argmin(np.abs(ch.x - xc)); dc = d[max(ic-3, 0):ic+3]
    print(f"SOD(m~rho,N={N}) {kw}: steps={s['steps']} dE={s['dE']:+.1e} mingap/h={s['mingap']:.3f} dt first/last 10%={s['dt_first']:.2e}/{s['dt_last']:.2e} dtmin={s['dtmin']:.1e} "
          f"lim={s['lim']} L1(rho,u,P)=({L1r:.4f},{L1u:.4f},{L1p:.4f}) post-shock rho bias={np.mean(rhon[sel])/rs-1:+.4f} P bias={np.mean(Pn[sel])/ps-1:+.4f} "
          f"std lnK plateau={lk.std():.4f}(nyq {nyq.std():.4f}, smooth {sm.std():.4f}) HQ={s['HQ']:.2e} Hhg={s['Hhg']:.2e} gaps@contact/D={np.round(dc/D,2)} {s['aborted'] or ''}")
    return ch


def noh(N=200, tend=0.6, P0=1e-6, amp=0.0, seed=1, **kw):
    rng = np.random.default_rng(seed)
    x = (np.arange(2*N)+0.5)/N - 1.0 + amp/N*rng.uniform(-1, 1, 2*N)
    v = -np.sign(x); rho = np.ones(2*N); P = np.full(2*N, P0)
    xp = np.concatenate([[2*x[0]-x[1]], x, [2*x[-1]-x[-2]]]); m = 0.5*(xp[2:]-xp[:-2])
    ch = Design1D(x, v, rho, P, m, bc='ghost', **kw); ch.run(tend)
    s = ch.summary()
    d, dL, dR, V, _, _ = ch.faces(ch.x); rhon = ch.m/V; Pn = ch.K*rhon**g
    xs = ch.t/3; sel = (np.abs(ch.x) < 0.7*xs) & (np.abs(ch.x) > 0.15*xs)
    Kex = (4/3)/4**g; lk = np.log(ch.K[sel]); nyq = 0.5*(lk[1:]-lk[:-1]); sm = 0.5*(lk[1:]+lk[:-1])
    xr = ch.x[ch.x > 0]; rr = rhon[ch.x > 0]
    xsh = xr[np.where(rr > 2.5)[0].max()] if np.any(rr > 2.5) else np.nan
    wid = np.sum((rr > 1.3) & (rr < 3.7))
    inner = rhon[np.argsort(np.abs(ch.x))[:4]].mean()
    print(f"NOH(N={N},amp={amp},P0={P0:.0e}) {kw}: steps={s['steps']} dE={s['dE']:+.1e} mingap/h={s['mingap']:.3f} dt first/last={s['dt_first']:.2e}/{s['dt_last']:.2e} dtmin={s['dtmin']:.1e} lim={s['lim']} "
          f"<rho>/4={np.mean(rhon[sel])/4:.4f} std rho/4={np.std(rhon[sel])/4:.4f} <P>/(4/3)={np.mean(Pn[sel])/(4/3):.4f} x_sh/x_ex={xsh/xs:.3f} width={wid} rho_wall={inner:.2f} "
          f"std lnK={lk.std():.4f}(nyq {nyq.std():.4f}, smooth {sm.std():.4f}) <lnK>-ex={lk.mean()-np.log(Kex):+.4f} HQ={s['HQ']:.3e} Hhg={s['Hhg']:.3e} {s['aborted'] or ''}")
    return ch


def pair(M, ratio=1.0, N=64, T=20.0, **kw):
    """periodic uniform chain at rest (rho=1,P=1); heavy particle i+1 with mass ratio (cell enlarged to keep rho); pair kicked with relative Mach M."""
    D = 1.0/N; c0 = np.sqrt(g)
    m = np.full(N, D); i = N//2; m[i+1] = ratio*D
    # positions: cells of width m (rho=1): x_k = cumulative
    V = m.copy(); xe = np.concatenate([[0.0], np.cumsum(V)]); x = 0.5*(xe[:-1] + xe[1:]); L = xe[-1]
    v = np.zeros(N); dl = M*c0
    # momentum-free kick: m_i v_i + m_{i+1} v_{i+1} = 0, v_i - v_{i+1} = dl
    v[i] = dl*m[i+1]/(m[i]+m[i+1]); v[i+1] = -dl*m[i]/(m[i]+m[i+1])
    ch = Design1D(x, v, np.ones(N), np.ones(N), m, bc='periodic', Lbox=L, **kw)
    KE0 = np.sum(0.5*m*v**2)
    ch.run(T*D/c0)
    s = ch.summary(); h = ch.hist
    d0 = D*0.5*(1+ratio)   # initial pair gap
    gap_pair = lambda x: (x[i+1]-x[i])
    print(f"PAIR M={M:5g} ratio={ratio:3g} {kw}: steps={s['steps']} {s['aborted'] or 'ok'} mingap/h(all)={s['mingap']:.3f} pair gap end/d0={gap_pair(ch.x)/d0:.3f} "
          f"dt first/last={s['dt_first']:.1e}/{s['dt_last']:.1e} dtmin/dt0={s['dtmin']/h[0,2]:.2e} lim={s['lim']} KE_end/KE0={h[-1,3]/KE0:.3f} dE={s['dE']:+.1e} HQ/KE0={s['HQ']/KE0:.3f} Hhg/KE0={s['Hhg']/KE0:.3f}")
    return ch


def nyquist(eps, N=64, T=8.0, **kw):
    D = 1.0/N; c0 = np.sqrt(g); x = D*(np.arange(N)+0.5); m = np.full(N, D)
    v = eps*c0*(-1.0)**np.arange(N)
    ch = Design1D(x, v, np.ones(N), np.ones(N), m, bc='periodic', Lbox=1.0, **kw)
    KE0 = np.sum(0.5*m*v**2)
    out = []
    for k in range(1, 9):
        ch.run(k*T/8*D/c0 - ch.t if hasattr(ch, 't') else T/8*D/c0)
        oe = np.abs(np.mean(ch.v*(-1.0)**np.arange(N)))/(eps*c0)
        out.append(oe)
        if ch.aborted: break
    s = ch.summary()
    d, *_ = ch.faces(ch.x)
    print(f"NYQ eps={eps} {kw}: |v_oe|/|v_oe0| at t=(1..8)/8*{T}D/c: {np.round(out,4)} mingap/h={s['mingap']:.3f} end gap oe/D={np.abs(d-D).max()/D:.3f} dE={s['dE']:+.1e} HQ/KE0={s['HQ']/KE0:.3f} Hhg/KE0={s['Hhg']/KE0:.3f} {s['aborted'] or ''}")


if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'id':
        identities()
    elif which == 'sod':
        sod(hg=None, qlaw=1, zface=False, zh=False, dtmode='legacy', cfl=0.25)      # baseline: cell + linear-ish TCF Q (reproduce 0.077)
        sod()                                                                        # design
        sod(hg='acou')                                                               # acoustic hg intermediate (step 3)
        sod(hg=None)                                                                 # no hourglass, face Z, QLAW2
        sod(hgscale=2.0)                                                             # Hugoniot with 4h|r| jump identification
        sod(N=400)
        sod(zface=False)
        sod(qlaw=1)
    elif which == 'noh':
        noh(hg=None, qlaw=1, zface=False, zh=False, dtmode='legacy')
        noh(); noh(hg='acou'); noh(hg=None); noh(hgscale=2.0)
        noh(amp=0.1); noh(amp=0.1, hg=None); noh(amp=0.1, hg='acou'); noh(amp=0.1, hg=None, qlaw=1, zface=False, zh=False, dtmode='legacy')
        noh(amp=0.1, seed=2); noh(amp=0.1, N=400)
        noh(P0=1e-10)
    elif which == 'pair':
        for ratio in (1.0, 8.0, 64.0):
            for M in (1, 10, 100, 1000):
                pair(M, ratio)
        for M in (10, 1000):
            pair(M, 1.0, hg=None); pair(M, 1.0, qon=False); pair(M, 1.0, hg='acou')
    elif which == 'nyq':
        for eps in (1e-3, 1.0):
            nyquist(eps); nyquist(eps, qon=False); nyquist(eps, hg=None); nyquist(eps, hg=None, qon=False)
