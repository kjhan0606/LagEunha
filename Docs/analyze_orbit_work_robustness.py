import json,re,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_orbit_work_robustness_20261018'
reference=b/'Codex_orbit_work_20261018/on'
def accepted(d):
 stages=[];hg=[];discarded=0;out=[]
 for line in (d/'log').read_text().splitlines():
  if line.startswith('[RETRY] step rejected:'):
   discarded+=len(stages);stages=[];hg=[]
  elif line.startswith('[ORBIT_WORK_STAGE]'):
   stages.append(dict(re.findall(r'(\w+)=([^\s]+)',line)))
  elif line.startswith('[HG_PERTURBATION]'):
   hg.append(dict(re.findall(r'(\w+)=([^\s]+)',line)))
  elif line.startswith('[ACCEPTED_STATE]'):
   state=dict(re.findall(r'(\w+)=([^\s]+)',line))
   assert state['finite']=='1'
   assert [int(s['stage']) for s in stages]==[1,2,3,4] and len(hg)==8
   out.append({'state':state,'stages':stages,'HG':hg});stages=[];hg=[]
 assert len(out)==1 and not stages and not hg
 return out[0],discarded
ref,tref=read(reference/'khout.000225.dat');base,_=accepted(reference)
results={}
for name in ['omp4','retry']:
 d=root/name;assert (d/'exit.code').read_text().strip()=='0'
 a,t=read(d/'khout.000225.dat');assert t==tref
 equal={f:bool(np.array_equal(a[f],ref[f])) for f in ref.dtype.names}
 assert all(equal.values())
 trial,n=accepted(d);assert int(trial['state']['retries'])==int(name=='retry')
 assert n==(4 if name=='retry' else 0)
 values=['wdt','hydro_minus_TCF','TCF','shear','total']
 diffs={k:max(abs(float(x[k])-float(y[k])) for x,y in zip(trial['stages'],base['stages'])) for k in values}
 assert all(v==0 for v in diffs.values())
 results[name]={'reader_fields_bitwise_equal':equal,'discarded_stage_records':n,'accepted_stage_max_abs_difference':diffs,'retries':int(trial['state']['retries'])}
out={'cases':results,'scope':'MPI1 one accepted step. OMP1 vs4 and forced final rejection .005 then .0025. Parser discards rejected stage/HG rows on retry markers. No production MPI or CUDA force extension and no long-time stability claim.'}
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
