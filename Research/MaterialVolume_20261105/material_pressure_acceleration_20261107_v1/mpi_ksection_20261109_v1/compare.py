from pathlib import Path
import json,subprocess,numpy as np,os,hashlib
O=Path(__file__).resolve().parent;R=O.parents[3];E=R/'Research/MaterialVolume_20261105/material_kh_stop_20261106_v1/extract_ids'
env=os.environ.copy();env.update(json.loads((R/'Research/MaterialVolume_20261105/material_finite_slack_20261107_v1/native_runtime_paths.json').read_text()))
def read(label):
 p=O/label;status=json.loads((p/'status.json').read_text());assert status['returncode']==0,(label,status)
 prefix='cylmat' if list(p.glob('cylmat.*.dat')) else 'rtout' if list(p.glob('rtout.*.dat')) else 'khout'
 files=sorted(p.glob(prefix+'.*.dat'));assert len(files)>1,label
 f=files[-1];out=p/'state_ids.bin'
 if not out.exists():subprocess.run([str(E),str(f),str(out)],env=env,check=True)
 a=np.fromfile(out,dtype='f8').reshape(-1,13);a=a[np.argsort(a[:,-1])]
 assert len(np.unique(a[:,-1]))==len(a) and np.isfinite(a).all() and (a[:,6:9]>0).all()
 return a,f
results=[]
for c in ['kh','rt','sedov','noh','cylinder']:
 label='maps' if c in ['noh','cylinder'] else 'suite';a,f=read(c+'_r1_'+label)
 for n in [1,2,3,4]:
  b,g=read(c+'_r'+str(n)+'_'+label)
  results.append({'case':c,'ranks':n,'cells':len(a),'max_difference':np.abs(b-a).max(axis=0).tolist(),'bitwise_equal':bool(np.array_equal(a.view('u8'),b.view('u8'))),'source':str(g),'sha256':hashlib.sha256(g.read_bytes()).hexdigest()})
  print(c,n,'max',np.max(np.abs(b-a)),'bitwise',results[-1]['bitwise_equal'])
 try:
  b,g=read(c+'_r1_legacy');print(c,'legacy_max',np.max(np.abs(b-a)))
  results.append({'case':c,'reference':'pre_MPI_production','bitwise_equal':bool(np.array_equal(a.view('u8'),b.view('u8'))),'max_difference':np.abs(b-a).max(axis=0).tolist()})
 except FileNotFoundError:pass
(O/'short_comparison.json').write_text(json.dumps(results,indent=2))
