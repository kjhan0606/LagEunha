import sys, numpy as np; sys.path.insert(0,'/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
from audit_design1d import Design1D, g
def pair_eq(M, ratio, N=64, T=20.0, **kw):
    """equal spacing D, uniform P=1; particle i+1 has mass ratio*D (density ratio, cold); relative kick M*c_light."""
    D=1.0/N; c0=np.sqrt(g); x=D*(np.arange(N)+0.5); m=np.full(N,D); i=N//2; m[i+1]=ratio*D
    rho=m/D; P=np.ones(N); v=np.zeros(N); dl=M*c0
    v[i]=dl*m[i+1]/(m[i]+m[i+1]); v[i+1]=-dl*m[i]/(m[i]+m[i+1])
    ch=Design1D(x,v,rho,P,m,bc='periodic',Lbox=1.0,**kw); KE0=np.sum(0.5*m*v**2)
    ch.run(T*D/c0); s=ch.summary(); h=ch.hist
    d=np.diff(ch.x); k=np.argmin(h[:,1]); 
    print(f"PAIR-EQ M={M:5g} ratio={ratio:3g} {kw}: steps={s['steps']} {s['aborted'] or 'ok'} mingap/h={s['mingap']:.3f} at t={h[k,0]/(D/c0):.2f}D/c; pair gap end/D={(ch.x[i+1]-ch.x[i])/D:.3f} "
          f"gaps around pair/D={np.round(d[i-2:i+3]/D,3)} dtmin/dt0={s['dtmin']/h[0,2]:.2e} dt last/first={s['dt_last']/s['dt_first']:.2f} KE_end/KE0={h[-1,3]/KE0:.3f} dE={s['dE']:+.1e} HQ/KE0={s['HQ']/KE0:.3f} Hhg/KE0={s['Hhg']/KE0:.3f}", flush=True)
w=sys.argv[1]
if w=='8':
    for M in (1,10,100,1000): pair_eq(M,8.0,T=(20.0 if M<=10 else 4.0),dtmode='geom')
elif w=='64':
    for M in (1,10,100,1000): pair_eq(M,64.0,T=(20.0 if M<=10 else 4.0),dtmode='geom')
elif w=='diag':
    # where was the 0.007 min gap in the uniform-density ratio-64 setup?
    from audit_run1d import pair
    import audit_run1d
    ch=pair(1,64.0,T=20.0,dtmode='geom')
    D=1/64; d=np.diff(ch.x); i=32
    print("  gaps/D around pair (i-3..i+3):", np.round(d[i-3:i+4]/D,3), " V_i/D, V_j/D:", np.round(ch.faces(ch.x)[3][i:i+2]/D,3))
