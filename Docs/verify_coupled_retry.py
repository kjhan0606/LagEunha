from pathlib import Path
import json,re,subprocess
root=Path('/gpfs/kjhan/LagForce/Codex_coupled_retry_20261017')
layout='/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json'

def record(log,tag):
 lines=[s for s in log.splitlines() if s.startswith('['+tag+']')]
 assert len(lines)==1,(tag,len(lines))
 return dict(re.findall(r'(\w+)=([^\s]+)',lines[0]))

def snapshot(d):
 paths=[p for p in d.iterdir() if p.name.startswith('khout.') and p.suffix=='.dat' and p.name!='khout.002305.dat']
 assert len(paths)==1,(d,paths)
 return paths[0]

out={'scope':'MPI1, real HG half/RK4/HG half, forced final rejection then half-dt retry versus direct half-dt. CPU force/HG; CPU or CUDA nearest-neighbor backend, OMP1/8. No long-time or multi-rank evolution claim.','cases':{}}
for backend in ['cpu','gpu']:
 for omp in [1,8]:
  name=f'{backend}{omp}';a=root/(name+'_direct');b=root/(name+'_retry')
  for d in [a,b]: assert (d/'exit.code').read_text().strip()=='0'
  report=root/(name+'_particles.json')
  subprocess.run(['python3',str(root/'compare_particles.py'),'--layout',layout,'--output',str(report),str(snapshot(a)),str(snapshot(b))],check=True)
  comparison=json.loads(report.read_text());assert comparison['all_physical_fields_bitwise_equal'],name
  la=(a/'log').read_text();lb=(b/'log').read_text()
  assert '[TEST_FINAL_REJECT]' not in la
  injection=record(lb,'TEST_FINAL_REJECT');assert float(injection['dt'])==.001 and int(injection['after_final_HG'])==1 and float(injection['elapsed_before_accept'])==0
  assert lb.count('[RETRY] step rejected:')==1
  sa=record(la,'ACCEPTED_STATE');sb=record(lb,'ACCEPTED_STATE')
  assert sa.pop('retries')=='0' and sb.pop('retries')=='1'
  assert sa==sb,(name,'accepted state differs',[(k,sa[k],sb.get(k)) for k in sa if sa[k]!=sb.get(k)])
  assert sa['step']=='1' and float(sa['elapsed'])==.0005 and float(sa['dt'])==.0005
  assert record(la,'ACCEPTED_INITIAL')==record(lb,'ACCEPTED_INITIAL')
  assert record(la,'ACCEPTED_RK_BUDGET')==record(lb,'ACCEPTED_RK_BUDGET')
  assert record(la,'ACCEPTED_HG_BOUNDARY')==record(lb,'ACCEPTED_HG_BOUNDARY')
  out['cases'][name]={'physical_fields_bitwise_equal':True,'one_rejected_trial':True,'accepted_ledger_equal_except_retry_count':True,'accepted_RK_budget_equal':True,'accepted_HG_boundary_equal':True,'accepted_steps':1,'accepted_elapsed':.0005,'t':comparison['t'],'dt':comparison['dt']}
(root/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
