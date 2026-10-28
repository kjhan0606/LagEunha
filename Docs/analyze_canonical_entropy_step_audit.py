from pathlib import Path
import numpy as np,json,struct,re,hashlib
root=Path(__file__).parent;old=Path('/gpfs/kjhan/LagForce/Codex_early_growth_budget_20261017');layout=json.loads(Path('/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json').read_text())
fields=['x','y','mass','vx','vy','den','pressure','ie','volume','w2','w2ceil','csound'];dtype=np.dtype({'names':['id']+fields,'formats':['<u8']+['<f8']*len(fields),'offsets':[layout['u4if']]+[layout[k] for k in fields],'itemsize':layout['size']})
def read(p):
 raw=p.read_bytes();magic,n,size,stress=struct.unpack_from('<iiii',raw);assert magic==-1 and size==dtype.itemsize and len(raw)==16+n*size+52
 a=np.frombuffer(raw,dtype=dtype,count=n,offset=16);ids=a['id']&layout['id_mask'];order=np.argsort(ids);assert np.array_equal(ids[order],np.arange(n));return a[order],struct.unpack_from('<d',raw,16+n*size)[0]
initial,zero=read(Path('/gpfs/kjhan/LagForce/Codex_canonical_initial_20261017/run/khout.000000.dat'));r0=np.hypot(initial['x']-3,initial['y']-3);rho0=initial['den']
regions={'inner_background':(rho0<.5)&(r0<=.5),'disk_inner_edge':(rho0>.5)&(r0<.6),'disk_core':(rho0>.5)&(r0>=.6)&(r0<1.9)}
out={'scope':'One actual accepted HG/RK4/HG step from existing evolved checkpoints. Per-cell log entropy change minus integrated TCF heat/U at four RK stages and exact HG log1p(dU/U). Residual includes time quadrature, geometric-volume evolution, roundoff and any unaccounted sources; not automatically pure RK truncation. Fixed INITIAL material bins. Frozen cells excluded. No long-time pass inferred.','cases':{}}
for entry in json.loads((root/'manifest.json').read_text()):
 name=entry['case'];d=root/name;assert (d/'exit.code').read_text().strip()=='0';step=entry['step'];start=d/f'khout.{step:06d}.dat';end=d/f'khout.{step+1:06d}.dat';a,t0=read(start);b,t1=read(end);ref,_=read(old/name/end.name)
 parity={k:bool(np.array_equal(b[k],ref[k])) for k in fields};assert all(parity.values()),(name,parity)
 gamma=float(np.float32(float(re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)',(d/'params.dat').read_text())[1])))
 lines=(d/'log').read_text().splitlines();accepted=[x for x in lines if x.startswith('[ACCEPTED_STATE]')];assert len(accepted)==1;state=dict(re.findall(r'(\w+)=([^\s]+)',accepted[0]));assert state['retries']=='0' and state['finite']=='1';dt=float(state['dt']);assert abs(t1-t0-dt)<2e-15
 rk=[x for x in lines if x.startswith('[ACCEPTED_RK_BUDGET]')];assert len(rk)==1;rk=dict(re.findall(r'(\w+)=([^\s]+)',rk[0]));assert rk['stage_mask_min']==rk['stage_mask_max']=='15'
 heat=np.loadtxt(d/'entropy_integral.dat');heat=heat[np.argsort(heat[:,0])];assert np.array_equal(heat[:,0],np.arange(len(a)));heat=heat[:,1]
 assert np.all(a['ie']>0)&np.all(b['ie']>0)&np.all(a['volume']>0)&np.all(b['volume']>0)
 ds=np.log1p((b['ie']-a['ie'])/a['ie'])+(gamma-1)*np.log1p((b['volume']-a['volume'])/a['volume']);res=ds-heat
 record={'start_time':t0,'dt':dt,'diagnostic_on_off_physical_fields_bitwise_equal':parity,'binary':(d/'binary.sha256').read_text().strip(),'start_sha256':hashlib.sha256(start.read_bytes()).hexdigest(),'end_sha256':hashlib.sha256(end.read_bytes()).hexdigest(),'accepted_budget_residual_E':float(rk['residual_E']),'regions':{}}
 for k,s in regions.items():
  assert np.all(np.hypot(a['x'][s]-3,a['y'][s]-3)<2.5)&np.all(np.hypot(b['x'][s]-3,b['y'][s]-3)<2.5)
  w=a['mass'][s];mean=lambda x:float(np.average(x[s],weights=w));rms=lambda x:float(np.sqrt(np.average(x[s]**2,weights=w)))
  record['regions'][k]={'cells':int(s.sum()),'measured_entropy_rate_mean':mean(ds)/dt,'heat_entropy_rate_mean':mean(heat)/dt,'unexplained_entropy_rate_mean':mean(res)/dt,'unexplained_entropy_rate_RMS':rms(res)/dt,'heat_entropy_rate_RMS':rms(heat)/dt,'unexplained_over_heat_RMS':rms(res)/max(rms(heat),1e-300),'entropy_residual_max_abs':float(np.max(abs(res[s]))),'negative_heat_cells':int(np.sum(heat[s]<-1e-18))}
 out['cases'][name]=record
(root/'result_v2.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:{'dt':v['dt'],'inner':v['regions']['inner_background']} for k,v in out['cases'].items()},indent=2))
