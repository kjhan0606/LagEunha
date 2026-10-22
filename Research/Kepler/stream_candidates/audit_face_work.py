from pathlib import Path
import json,sys
import numpy as np
from audit_snapshot import read_snapshot
root=Path(__file__).resolve().parent
if sys.argv[1]=='INITIAL_GRID':
 ids0=np.arange(256*256,dtype=np.uint64);h=24/256
 x0=(ids0%256+.5)*h;y0=(ids0//256+.5)*h
 r0=np.hypot(x0-12,y0-12);omega=(r0*r0+.047**2)**(-.75)
 p={'u4if':ids0,'x':x0,'y':y0,'vx':-(y0-12)*omega,'vy':(x0-12)*omega,
    'volume':np.full(len(ids0),h*h),'ie':np.full(len(ids0),np.float32(np.float32(1e-6)*np.float32(h*h)/(np.float32(1.666667)-1)))}
 m={'t':0.,'Lx':24.,'Ly':24.}
else:
 p,m=read_snapshot(sys.argv[1])
ncells=len(p['x'])
f=np.loadtxt(sys.argv[2],skiprows=1)
ids=p['u4if']&np.uint64((1<<56)-1)
lookup=np.full(int(ids.max())+1,-1,dtype=np.int64); lookup[ids.astype(int)]=np.arange(len(ids))
if np.any(f[:,6]!=0): raise ValueError('Wall faces require explicit wall kinematics')
i=lookup[f[:,4].astype(int)]; j=lookup[f[:,5].astype(int)]
if (i<0).any() or (j<0).any(): raise ValueError('Face refers to a missing real ID')
x,y=p['x'].astype(float),p['y'].astype(float); vx,vy=p['vx'].astype(float),p['vy'].astype(float)
dx=x[j]-x[i];dy=y[j]-y[i];dx-=m['Lx']*np.rint(dx/m['Lx']);dy-=m['Ly']*np.rint(dy/m['Ly'])
d2=dx*dx+dy*dy
cf_x=.5*(f[:,12]+f[:,14])-x[i]; cf_y=.5*(f[:,13]+f[:,15])-y[i]
cf_x-=m['Lx']*np.rint(cf_x/m['Lx']);cf_y-=m['Ly']*np.rint(cf_y/m['Ly']);sx,sy=f[:,10],f[:,11]
def geo_rate(ux,uy):
 corr=((ux[j]-ux[i])*(cf_x-.5*dx)+(uy[j]-uy[i])*(cf_y-.5*dy))*(dx*sx+dy*sy)/d2
 return .5*((ux[i]+ux[j])*sx+(uy[i]+uy[j])*sy)-corr
geo=geo_rate(vx,vy); work=f[:,8]
r=np.hypot(x-12,y-12); Om=(r*r+.047**2)**(-.75)
circ=geo_rate(-(y-12)*Om,(x-12)*Om)
reduce=lambda val:np.bincount(i,weights=val,minlength=ncells)
V=p['volume'].astype(float);U=p['ie'].astype(float);gamma=1.666667
Pi=(gamma-1)*U/V
Vg,Vw,Vc=reduce(geo),reduce(work),reduce(circ)
Ud=reduce(f[:,9]+f[:,7]*(vx[i]*sx+vy[i]*sy))
energy_closure=Pi*(Vg-Vw)
other=Ud+Pi*Vg-energy_closure
entropy=np.divide(Ud,U,out=np.zeros_like(Ud),where=U>0)+(gamma-1)*Vg/V
closed=np.hypot(reduce(sx),reduce(sy))
# For v=(-y,x)*(r^2+eps^2)^(-3/4), the streamfunction is
# psi=2*(r^2+eps^2)^(1/4). Its exact outward flux along the oriented
# straight face is psi(vertex1)-psi(vertex2); no quadrature is needed.
psi1=2*((f[:,12]-12)**2+(f[:,13]-12)**2+.047**2)**.25
psi2=2*((f[:,14]-12)**2+(f[:,15]-12)**2+.047**2)**.25
exact_div_flux=reduce(psi1-psi2)
out={'time':m['t'],'snapshot':str(Path(sys.argv[1]).resolve()),'faces':len(f),'assumptions':'W=0, Re=infinity, native phase1 face work. Matched full geometry; periodic minimum-image neighbors. Circular-shadow velocity is postprocessing only.','max_surface_closure':float(closed.max()),'face_work_geo_max_abs':float(abs(work-geo).max()),'faces_work_geo_difference_gt_1e-12':int((abs(work-geo)>1e-12).sum()),'radial':[]}
for lo,hi in [(0,.5),(.5,1),(1,2),(2,4),(4,8)]:
 s=(r>=lo)&(r<hi)
 fs=s[i]
 rms=lambda a:float(np.sqrt(np.mean(a[s]**2)))
 out['radial'].append({'Rlo':lo,'Rhi':hi,'n':int(s.sum()),'rms_Vdot_geom_over_V':rms(Vg/V),'rms_work_geometry_mismatch_over_V':rms((Vg-Vw)/V),'rms_circular_geometry_Vdot_over_V':rms(Vc/V),'rms_entropy_rate':rms(entropy),'sum_abs_pressure_volume_mismatch_rate':float(abs(energy_closure[s]).sum()),'sum_abs_other_thermal_rate':float(abs(other[s]).sum()),'sum_thermal_rate':float(Ud[s].sum()),'faces_work_geo_difference_gt_1e-12':int(((abs(work-geo)>1e-12)&fs).sum())})
 out['radial'][-1].update({'max_exact_circular_divergence_flux_over_V':float(np.max(abs(exact_div_flux[s]/V[s]))), 'max_pstar_over_initial_P':float(np.max(f[fs,7])/1e-6),'median_pstar_over_initial_P':float(np.median(f[fs,7])/1e-6),'sum_positive_extra_thermal_rate':float(np.maximum(other[s],0).sum()), 'internal_energy':float(U[s].sum())})
Path(sys.argv[3]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
