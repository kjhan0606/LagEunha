from pathlib import Path
import numpy as np,json,struct
root=Path('/gpfs/kjhan/LagForce/Codex_t2_onset_rhs_20261017')
for case in ['t2','circular']:assert (root/case/'exit.code').read_text().strip()=='0'
a=np.loadtxt(root/'t2/rhs.dat');a=a[np.argsort(a[:,0])];b=np.loadtxt(root/'circular/rhs.dat');b=b[np.argsort(b[:,0])]
assert a.shape==b.shape==(36864,25) and np.isfinite(a).all() and np.isfinite(b).all()
assert np.array_equal(a[:,[0,1,2,5,6,7]],b[:,[0,1,2,5,6,7]])
layout=json.loads(Path('/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json').read_text());raw=Path('/gpfs/kjhan/LagForce/Codex_canonical_initial_20261017/run/khout.000000.dat').read_bytes();_,n,size,_=struct.unpack_from('<iiii',raw);names=['u4if','x','y','den','mass'];dtype=np.dtype({'names':names,'formats':['<u8']+['<f8']*4,'offsets':[layout[k] for k in names],'itemsize':size});ic=np.frombuffer(raw,dtype=dtype,count=n,offset=16);ic=ic[np.argsort(ic['u4if']&layout['id_mask'])];r0=np.hypot(ic['x']-3,ic['y']-3)
pdata=json.loads((root/'pressure.json').read_text());gamma=pdata['gamma'];p0=pdata['pressure0'];pn=a[:,5]/p0
circ=-gamma*pn*b[:,12]/a[:,6];pert=-gamma*pn*(a[:,12]-b[:,12])/a[:,6];total=circ+pert
out={}
for name,sel in [('inner_background',(ic['den']<.5)&(r0<=.5)),('disk_inner_edge',(ic['den']>.5)&(r0<.6))]:
 w=ic['mass'][sel];mean=lambda x:float(np.average(x[sel],weights=w));fluct=pn-mean(pn)
 out[name]={'n':int(sel.sum()),'compression_mean_rate':mean(total),'circular_velocity_compression_mean_rate':mean(circ),'velocity_perturbation_compression_mean_rate':mean(pert),'compression_variance_growth':2*mean(fluct*total),'circular_velocity_variance_growth':2*mean(fluct*circ),'velocity_perturbation_variance_growth':2*mean(fluct*pert),'rms_circular_Vdot_over_V':np.sqrt(mean((b[:,12]/a[:,6])**2)),'rms_actual_Vdot_over_V':np.sqrt(mean((a[:,12]/a[:,6])**2))}
 assert abs(out[name]['compression_variance_growth']-pdata['regions'][name]['compression_variance_growth'])<1e-12
result={'scope':'Instantaneous geometric volume-rate split at identical t2 positions, masses, pressures and volumes. Diagnostic velocity replacement only; no time evolution. Material bins defined by t0 IDs. Circular-velocity contribution is a kinematic diagnostic, not proof of the entire instability mechanism.','regions':out}
(root/'volume_split.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
