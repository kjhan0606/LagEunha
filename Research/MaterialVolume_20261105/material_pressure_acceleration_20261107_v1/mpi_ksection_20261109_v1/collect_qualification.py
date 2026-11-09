from pathlib import Path
import hashlib,json,re,struct,subprocess,os,numpy as np
O=Path(__file__).resolve().parent;R=O.parents[3]
E=R/'Research/MaterialVolume_20261105/material_kh_stop_20261106_v1/extract_ids'
env=os.environ.copy();env.update(json.loads((R/'Research/MaterialVolume_20261105/material_finite_slack_20261107_v1/native_runtime_paths.json').read_text()))
fields=['x','y','vx','vy','ax','ay','mass','U','V','W','pi','inertia','ID'];runs={}
def read(label):
 d=O/label;status=json.loads((d/'status.json').read_text());assert status['returncode']==0
 prefix='cylmat' if list(d.glob('cylmat.*.dat')) else 'rtout' if list(d.glob('rtout.*.dat')) else 'khout'
 files=sorted(d.glob(prefix+'.*.dat'));assert len(files)>1
 f=files[-1];export=d/'state_ids.bin'
 if not export.exists():subprocess.run([str(E),str(f),str(export)],env=env,check=True)
 a=np.fromfile(export,'f8').reshape(-1,13);a=a[np.argsort(a[:,-1])]
 assert np.isfinite(a).all() and (a[:,6:9]>0).all() and len(np.unique(a[:,-1]))==len(a)
 with f.open('rb') as x:x.seek(-52,2);t=struct.unpack('d',x.read(8))[0]
 assert abs(t-status['target'])<1e-13,(label,t,status['target'])
 text=(d/'log').read_text();matches=re.findall(r'\[MATERIAL_SSP_ACCEPT\] t=(\S+).*?retries=(\d+)',text)
 accepted=sorted(set((float(t),int(n)) for t,n in matches));events=sorted(set(re.findall(r'\[CYLINDER_MATERIAL_EVENT\].*',text)))
 runs[label]={'status':status,'output':str(f.relative_to(R)),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'saved_time':t,'cells':len(a),'accepted_steps':len(accepted),'maximum_retries':max((n for _,n in accepted),default=0),'boundary_events':events}
 return a
pairs=[('kh_multi_r1','kh_multi_r4'),('kh_restart_r1','kh_restart_r2'),('cyl_event_r1','cyl_event_r4')]
for c in ['sod_final','shear_final','gresho_final','noh_stages','rt_final','sedov_final','kh_reject_verified']:
 pairs.append((c+'_r1',c+'_r4'))
for n in [2,4]:pairs.append(('cyl_restart_r1','cyl_restart_r'+str(n)))
comparison=[]
for x,y in pairs:
 a,b=read(x),read(y);assert a.shape==b.shape
 error=abs(a-b).max(axis=0);same=np.array_equal(a.view('u8'),b.view('u8'))
 assert np.array_equal(a[:,-1],b[:,-1]) and np.array_equal(a[:,6],b[:,6])
 assert np.allclose(a,b,rtol=1e-11,atol=1e-12),(x,y,error)
 comparison.append({'a':x,'b':y,'bitwise_equal':bool(same),'maximum_absolute_difference':dict(zip(fields,error.tolist()))})
 print(x,y,'bitwise',same,'max',error.max())
assert runs['kh_reject_verified_r1']['maximum_retries']==runs['kh_reject_verified_r4']['maximum_retries']==1
report={'scope':'MPI correctness qualification of unchanged native 2D material SSPRK2; not long-term physics or strong-scaling validation.','fields':fields,'runs':runs,'comparisons':comparison,'builds':{}}
for label in ['build_final','build_reject']:
 p=O/label/'native.exe';report['builds'][label]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
assert report['builds']['build_final']['sha256']!=report['builds']['build_reject']['sha256']
(O/'qualification_report.json').write_text(json.dumps(report,indent=2)+'\n')
