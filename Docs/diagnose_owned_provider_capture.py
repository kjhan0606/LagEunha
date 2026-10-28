"""Bounded comparison witness for the immutable owned-provider component."""
import argparse,json,time
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);args=p.parse_args()
base=Path('/gpfs/kjhan/LagForce');run=args.run;m=json.loads((run/'provider_rank0.json').read_text())
root=Path(m['buffer_root']);velocity=base/'Codex_full_disk_fine_second_midpoint_velocity_gravity_20261019';density=base/'Codex_full_disk_fast_density_replay_20261019'
def lines(path):
 with path.open() as f:return [json.loads(s) for s in f]
rows={i:json.loads((root/f'rows_{i}.json').read_text()) for i in (2,3,4)}
fits=lines(density/'cells.jsonl');vel=lines(velocity/'cells.jsonl')
with np.load(velocity/'prepared.npz') as f:original={k:f[k] for k in f.files}
with np.load(json.loads((root/'context.json').read_text())['state']) as f:state={k:f[k] for k in f.files}
report=dict(state_matches={k:bool(np.array_equal(original[k],state[k])) for k in ('state','mass','ids','factors','Psi')},
 density_first_difference=None,velocity_first_difference=None)
for i,(a,b) in enumerate(zip(rows[3],fits)):
 for key in ('id','status','prior','parameter','nodes24'):
  if a[key]!=b[key]:
   report['density_first_difference']=dict(i=i,key=key,actual=a[key],expected=b[key]);break
 if report['density_first_difference']:break
for i,(a,b) in enumerate(zip(rows[4],vel)):
 for key in ('id','status','metadata','diagnostics','generator','gravity','gravity_absolute_sum','moment_defect','points'):
  if a[key]!=b[key]:
   report['velocity_first_difference']=dict(i=i,key=key,actual=a[key],expected=b[key]);break
 if report['velocity_first_difference']:break
ng=m['gravity_nodes'];gn=np.memmap(root/'gravity.bin',mode='r',dtype=np.float64,shape=(ng,7));expected=np.load(velocity/'gravity_nodes.npy',mmap_mode='r')
report.update(actual_nodes_first=gn[:4].tolist(),expected_nodes_first=expected[:4].tolist(),actual_shape=list(gn.shape),expected_shape=list(expected.shape),column_unequal_counts_first262144=np.count_nonzero(gn[:262144]!=expected[:262144],axis=0).tolist(),maximum_absolute_difference_first262144=np.max(abs(gn[:262144]-expected[:262144]),axis=0).tolist())
(run/'difference_witness.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
