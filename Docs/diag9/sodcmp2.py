import sys, numpy as np
sys.path.insert(0,'/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad'); import lagsnap
src=open('/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad/exactcmp.py').read()
exec(src[:src.index('mode=sys.argv[1]')])
for p in sys.argv[1:]:
    s_=lagsnap.read(p); t,x,vx,rho,P=s_['t'],s_['x'],s_['vx'],s_['den'],s_['pressure']
    s=np.abs(x-0.75)<0.2; ex=riemann(1,0,1,0.125,0,0.1,(x[s]-0.75)/t)
    e=np.mean(np.abs(np.stack([rho[s],vx[s],P[s]],1)-ex),0)
    K=P[s]/rho[s]**G; Kx=ex[:,2]/ex[:,0]**G
    print(f"{p.split('/')[-2]}: t={t:.3f} L1 rho={e[0]:.4f} u={e[1]:.4f} P={e[2]:.4f} max u={vx[s].max():.3f} (exact {ex[:,1].max():.3f})  L1 lnK={np.mean(np.abs(np.log(K/Kx))):.4f}")
