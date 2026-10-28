"""Matched off/on eight-step entropy attribution; not a t0-to-t10 history."""
import json,re,sys,subprocess
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
base=Path('/gpfs/kjhan/LagForce');root=base/'Codex_entropy_multistep_20261018'
old=base/'Codex_mpi_coupled_step_20261017'
paths={}
for case in ['off','on']:
 d=root/case;assert (d/'exit.code').read_text().strip()=='0'
 states=[dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in (d/'log').read_text().splitlines() if line.startswith('[ACCEPTED_STATE]')]
 assert len(states)==8 and all(x['retries']=='0' and x['finite']=='1' for x in states)
 outputs=[p for p in d.iterdir() if p.name.startswith('khout.') and not p.is_symlink()]
 assert len(outputs)==1;paths[case]=outputs[0]
subprocess.run(['python3',str(old/'compare_particles.py'),'--layout',str(base/'Codex_DPQ_backend_parity_20261017/layout.json'),
 '--output',str(root/'comparison.json'),str(paths['off']),str(paths['on'])],check=True)
comparison=json.loads((root/'comparison.json').read_text());assert comparison['all_physical_fields_bitwise_equal']
a,t0=read(root/'on/khout.002305.dat');b,t1=read(paths['on']);ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
assert np.array_equal(a['mass'],b['mass']) and np.array_equal(a['mass'],ic['mass'])
q=np.loadtxt(root/'on/entropy.dat');q=q[np.argsort(q[:,0])]
qoff=np.loadtxt(root/'off/entropy.dat');qoff=qoff[np.argsort(qoff[:,0])]
assert np.array_equal(q[:,0],np.arange(len(a))) and np.array_equal(q[:,1],qoff[:,1])
gamma=float(np.float32(1.666667));ds=np.log(b['ie']/a['ie'])+(gamma-1)*np.log(b['volume']/a['volume'])
r=np.hypot(ic['x']-3,ic['y']-3)
masks={'inner_background':(r<.5)&(ic['den']<.5),'inner_disk_edge':(r<.6)&(ic['den']>.5),'disk_core':(r>=.6)&(r<1.9)&(ic['den']>.5)}
out={'start_time':t0,'end_time':t1,'duration':t1-t0,'start_sha256':sha(root/'on/khout.002305.dat'),'end_sha256':sha(paths['on']),
 'physical_fields_on_off_bitwise_equal':True,'total_entropy_on_off_bitwise_equal':True,
 'total_minus_components_max_abs':float(np.max(abs(q[:,1]-q[:,2]-q[:,3]))),'regions':{},
 'scope':'Same trajectory, eight accepted steps. Initial-ID material bins. Fraction is of modeled log-entropy, not energy or the full prior history. Residual includes quadrature, geometric volume and roundoff.'}
for name,mask in masks.items():
 avg=lambda v:float(np.average(v[mask],weights=ic['mass'][mask]))
 out['regions'][name]={'cells':int(mask.sum()),'measured_ds_mean':avg(ds),'TCF_ds_mean':avg(q[:,2]),'HG_ds_mean':avg(q[:,3]),
  'HG_fraction_of_modeled_entropy':avg(q[:,3])/avg(q[:,1]),'closure_residual_mean':avg(ds-q[:,1]),
  'closure_residual_RMS':avg((ds-q[:,1])**2)**.5,'modeled_ds_RMS':avg(q[:,1]**2)**.5}
with (root/'verification.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
