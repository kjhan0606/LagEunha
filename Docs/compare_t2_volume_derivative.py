from pathlib import Path
import json,struct,numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_t2_volume_derivative_20261017');baseline=Path('/gpfs/kjhan/LagForce/Codex_t2_onset_rhs_20261017/circular')
def rhs(d):
 a=np.loadtxt(d/'rhs.dat');a=a[np.argsort(a[:,0])];assert a.shape==(36864,25) and np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(36864));return a
base=rhs(baseline);gamma=1.6666669845581055;p0=2.7777796321447566e-5
layout=json.loads(Path('/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json').read_text());raw=Path('/gpfs/kjhan/LagForce/Codex_canonical_initial_20261017/run/khout.000000.dat').read_bytes();_,n,size,_=struct.unpack_from('<iiii',raw);names=['u4if','x','y','den','mass'];dtype=np.dtype({'names':names,'formats':['<u8']+['<f8']*4,'offsets':[layout[k] for k in names],'itemsize':size});ic=np.frombuffer(raw,dtype=dtype,count=n,offset=16);ic=ic[np.argsort(ic['u4if']&layout['id_mask'])];r0=np.hypot(ic['x']-3,ic['y']-3);regions={'initial_inner_background':(ic['den']<.5)&(r0<=.5),'initial_disk_inner_edge':(ic['den']>.5)&(r0<.6),'initial_disk_core':(r0>=.6)&(r0<1.9)}
out=[]
for idx,h in enumerate([.0001,.00005,.000025]):
 for sign in ['plus','minus']:assert (root/f'h{idx}_{sign}'/'exit.code').read_text().strip()=='0'
 plus=rhs(root/f'h{idx}_plus');minus=rhs(root/f'h{idx}_minus');assert np.array_equal(plus[:,7],base[:,7]) and np.array_equal(minus[:,7],base[:,7])
 fd=(plus[:,6]-minus[:,6])/(2*h);pred=base[:,12];err=fd-pred
 pp=base[:,5]/p0*(base[:,6]/plus[:,6])**gamma;pm=base[:,5]/p0*(base[:,6]/minus[:,6])**gamma;pn=base[:,5]/p0;pdot=-gamma*pn*pred/base[:,6]
 report={'h':h,'regions':{}}
 for name,sel in regions.items():
  w=ic['mass'][sel];mean=lambda q:float(np.average(q[sel],weights=w));variance=lambda q:mean((q-mean(q))**2)
  variance_fd=(variance(pp)-variance(pm))/(2*h);variance_rhs=2*mean((pn-mean(pn))*pdot)
  report['regions'][name]={'n':int(sel.sum()),'Vdot_relative_L1_error':float(np.abs(err[sel]).sum()/np.abs(pred[sel]).sum()),'max_Vdot_over_V_abs_error':float(np.max(np.abs(err[sel]/base[sel,6]))),'Vdot_over_V_rms_fd':float(np.sqrt(mean((fd/base[:,6])**2))),'Vdot_over_V_rms_rhs':float(np.sqrt(mean((pred/base[:,6])**2))),'pressure_variance_rate_fd_adiabatic':variance_fd,'pressure_variance_rate_rhs_compression':variance_rhs,'pressure_variance_rate_difference':variance_fd-variance_rhs,'worst_id':int(base[sel,0][np.argmax(np.abs(err[sel]/base[sel,6]))])}
 out.append(report)
result={'scope':'Central derivative of actual rebuilt Voronoi volumes under exact circular kinematics, compared with baseline geometric Vdot. Adiabatic pressure values reconstructed from volumes solely to check the compression chain rule. No energy evolution or long-term stability claim. Fixed original material bins; remote periodic/frozen boundaries excluded.','samples':out}
(root/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
