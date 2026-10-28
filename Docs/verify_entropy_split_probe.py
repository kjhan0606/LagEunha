import json,subprocess,re
from pathlib import Path
import numpy as np
root=Path('/gpfs/kjhan/LagForce/Codex_entropy_split_20261018')
old=Path('/gpfs/kjhan/LagForce/Codex_mpi_coupled_step_20261017')
result={}
arrays={}
for case in ['off','on','retry']:
 d=root/case;assert (d/'exit.code').read_text().strip()=='0'
 report=d/'same_binary_comparison.json'
 subprocess.run(['python3',str(old/'compare_particles.py'),'--layout','/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json',
 '--output',str(report),str(root/'off/khout.002306.dat'),str(d/'khout.002306.dat')],check=True)
 comp=json.loads(report.read_text());assert comp['all_physical_fields_bitwise_equal']
 a=np.loadtxt(d/'entropy.dat');a=a[np.argsort(a[:,0])];assert np.array_equal(a[:,0],np.arange(36864)) and np.isfinite(a).all()
 arrays[case]=a
 states=[s for s in (d/'log').read_text().splitlines() if s.startswith('[ACCEPTED_STATE]')];assert len(states)==1
 state=dict(re.findall(r'(\w+)=([^\s]+)',states[0]));assert int(state['retries'])==int(case=='retry')
 result[case]={'physical_fields_match_same_binary_off':True,'retries':int(state['retries']),'snapshot_sha256':comp['candidate_sha256']}
assert np.array_equal(arrays['off'][:,1],arrays['on'][:,1])
assert np.array_equal(arrays['on'],arrays['retry'])
a=arrays['on'];closure=np.max(abs(a[:,1]-a[:,2]-a[:,3]));assert closure<1e-14
result.update(total_on_off_bitwise_equal=True,retry_components_bitwise_equal=True,component_closure_max_abs=float(closure),
 scope='One MPI1 split-HG accepted step from t10; on/off physical fields compared within new binary; old binary differs at roundoff level. No integrated historical t0-to-t12 attribution.')
with (root/'verification.json').open('x') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
