from pathlib import Path
import numpy as np,json,struct,hashlib
root=Path(__file__).parent;base=Path('/gpfs/kjhan/LagForce');layout=json.loads((base/'Codex_DPQ_backend_parity_20261017/layout.json').read_text())
names=['u4if','x','y','mass','pressure','den'];dt=np.dtype({'names':names,'formats':['<u8']+['<f8']*5,'offsets':[layout[k] for k in names],'itemsize':layout['size']})
def read(p):
 raw=p.read_bytes();magic,n,size,stress=struct.unpack_from('<iiii',raw);assert magic==-1 and size==dt.itemsize and stress==1 and len(raw)==16+n*size+52
 a=np.frombuffer(raw,dtype=dt,count=n,offset=16);ids=a['u4if']&layout['id_mask'];order=np.argsort(ids);assert np.array_equal(ids[order],np.arange(n));a=a[order]
 return a,struct.unpack_from('<d',raw,16+n*size)[0],hashlib.sha256(raw).hexdigest()
initial,_,ih=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');r0=np.hypot(initial['x']-3,initial['y']-3);p0=initial['pressure'][0];assert np.all(initial['pressure']==p0)
disk=initial['den']>.5;core=disk&(r0>=.6)&(r0<1.9);mass=initial['mass'];results=[]
run=base/'Codex_HG_constrained_20261017/omp8_t20'
for filename in ['khout.002305.dat','khout.004114.dat','khout.006248.dat']:
 p=run/filename;a,t,sha=read(p);assert np.array_equal(a['mass'],mass);err=np.abs(a['pressure']/p0-1)
 result={'snapshot':str(p),'sha256':sha,'time':t,'mass_max_abs_change':0.,'core_thresholds':{},'radial_bins':[]}
 for threshold in [.01,.1,1.]:
  sel=core&(err>threshold);result['core_thresholds'][str(threshold)]={'count':int(sel.sum()),'fraction_cells':float(sel.sum()/core.sum()),'fraction_mass':float(mass[sel].sum()/mass[core].sum())}
 for lo in np.arange(.5,2.,.05):
  hi=lo+.05;sel=disk&(r0>=lo)&(r0<hi)
  if not sel.any():continue
  w=mass[sel];result['radial_bins'].append({'r0_lo':float(lo),'r0_hi':float(hi),'count':int(sel.sum()),'pressure_relative_RMS':float(np.sqrt(np.average(err[sel]**2,weights=w))),'pressure_relative_max':float(err[sel].max()),'mass_fraction_error_gt_0p1':float(mass[sel&(err>.1)].sum()/w.sum())})
 results.append(result)
out={'scope':'Fixed initial material membership, original fixed-frozen-HG long run. Absolute pressure deviation from initial uniform P. Core=initial density>0.5 and 0.6<=r0<1.9. Threshold crossings are error diagnostics, not a physical shock definition.','initial_sha256':ih,'pressure0':float(p0),'snapshots':results}
(root/'radial.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([{k:x[k] for k in ['time','core_thresholds']} for x in results],indent=2))
