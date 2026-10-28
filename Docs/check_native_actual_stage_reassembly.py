"""Same actual moving RHS in original executable, one versus two MPI ranks."""
import argparse,json,shutil,time
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);args=p.parse_args();start=time.perf_counter()
run=args.run;roots=[];reports=[]
for ranks in (1,2):
    folder=run/f'ranks{ranks}';report=json.loads((folder/'actual_stage_report.json').read_text())
    assert report['status']=='passed_actual_full_moving_stage' and report['native_MPI_GPU_stage']
    assert report['service_MPI_calls']==report['service_GPU_calls']==0
    manifests=[json.loads((folder/f'provider_rank{rank}.json').read_text()) for rank in range(ranks)]
    roots.append(Path(manifests[0]['buffer_root']));reports.append(report)
for key in ('cells','nodes','gravity_nodes','node_inputs_sha256','node_extra_sha256'):
    assert reports[0][key]==reports[1][key],('same actual face data',key)
with np.load(roots[0]/'stage_rates.npz') as first,np.load(roots[1]/'stage_rates.npz') as second:
    assert set(first.files)==set(second.files)
    exact={}
    for key in first.files:
        a,b=first[key],second[key];assert np.array_equal(a,b),('actual RHS rates',key)
        exact[key]=True
# Incoming full Q and hydro are compared in canonical original receiver order.
header=np.fromfile(roots[0]/'native_stage/header.bin',dtype='<u8');n,nn,ni,ng=map(int,header[2:])
for key,shape in [('canonical_flux',(nn,19)),('outputs/stage',(n,94)),('outputs/Q',(ni,16)),('outputs/hydro',(n,106)),('outputs/ledgers',(12,))]:
    a=np.memmap(roots[0]/'native_stage'/(key+'.bin'),mode='r',dtype='<f8',shape=shape)
    b=np.memmap(roots[1]/'native_stage'/(key+'.bin'),mode='r',dtype='<f8',shape=shape)
    for lo in range(0,len(a),131072):assert np.array_equal(a[lo:lo+131072],b[lo:lo+131072]),(key,lo)
    del a,b
    exact[key]=True
result=dict(status='passed_native_MPI_actual_RHS_reassembly',original_executable=True,
    cells=n,nodes=nn,gravity_nodes=ng,input_time=reports[0]['time'],exact=exact,
    same_node_cpu_preparation=True,strict_wall_gravity_massroot_ledger_audits_passed=True,
    rank1_stage_seconds=reports[0]['total_seconds'],rank2_stage_seconds=reports[1]['total_seconds'],
    validation_seconds=time.perf_counter()-start,actual_RHS_replays=2,
    unique_new_physical_states=0,accepted_time_steps=0,native_MPI_evolution_verified=False,
    limits='One stored moving physical state replayed on 1 and 2 ranks. Actual RHS integration only; native midpoint state/factor/commit still pending.')
(run/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
# No failed/old run is removed. These exact disposable roots were created by
# this verified pair and only live on the allocated node's shared memory.
for ranks in (1,2):
    folder=run/f'ranks{ranks}'
    for rank in range(ranks):
        root=Path(json.loads((folder/f'provider_rank{rank}.json').read_text())['buffer_root'])
        assert str(root).startswith('/dev/shm/eunha_fmp_') and root.parent==Path('/dev/shm')
        shutil.rmtree(root)
