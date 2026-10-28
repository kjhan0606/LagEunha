"""Compare native MPI-owned preparation to the immutable moving-stage capture.

Runs inside the allocated compute node after the original executable exits.
No new RHS/evolution; on success explicitly release only owned /dev/shm files.
"""
import argparse,json,os,time
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--ranks',type=int,required=True);p.add_argument('--same-host-reference',type=Path)
    args=p.parse_args();started=time.perf_counter()
    base=Path('/gpfs/kjhan/LagForce');run=args.run
    manifests=[json.loads((run/f'provider_rank{rank}.json').read_text()) for rank in range(args.ranks)]
    assert all(m['status']=='passed_owned_cpu_MPI_preparation' and m['rank']==i and m['ranks']==args.ranks for i,m in enumerate(manifests))
    density=base/'Codex_full_disk_fast_density_replay_20261019'
    velocity=base/'Codex_full_disk_fine_second_midpoint_velocity_gravity_20261019'
    def jsonl(path):
        with path.open() as f:return [json.loads(line) for line in f]
    if args.same_host_reference:
        ref=args.same_host_reference;ref_report=json.loads((ref/'report.json').read_text())
        assert ref_report['status']=='passed_original_same_host_cpu_reference' and ref_report['host']==os.uname().nodename
        baseline_bounds=json.loads((ref/'density_bounds.json').read_text())
        baseline_fits=json.loads((ref/'density_rows.json').read_text());baseline_vel=json.loads((ref/'velocity_rows.json').read_text())
        with np.load(ref/'prepared.npz') as f:prepared={k:f[k] for k in f.files}
        original=np.memmap(ref_report['gravity_file'],mode='r',dtype=np.float64,shape=(ref_report['gravity_nodes'],7))
    else:
        baseline_bounds=json.loads((density/'density_bounds.json').read_text())
        baseline_fits=jsonl(density/'cells.jsonl');baseline_vel=jsonl(velocity/'cells.jsonl')
        with np.load(velocity/'prepared.npz') as f:prepared={k:f[k] for k in f.files}
        original=np.load(velocity/'gravity_nodes.npy',mmap_mode='r')
    n=manifests[0]['cells'];assert len(baseline_fits)==len(baseline_vel)==n
    rows={command:[] for command in (2,3,4)};node_count=0;rank_rows=[]
    for m in manifests:
        root=Path(m['buffer_root']);assert str(root).startswith('/dev/shm/eunha_fmp_')
        lo=m['owned_start'];hi=lo+m['owned_count'];assert lo==n*m['rank']//args.ranks and hi==n*(m['rank']+1)//args.ranks
        for command in rows:rows[command]+=json.loads((root/f'rows_{command}.json').read_text())
        ng=m['gravity_nodes'];start=int(prepared['gravity_offsets'][lo]);end=int(prepared['gravity_offsets'][hi]);assert end-start==ng
        gn=np.memmap(root/'gravity.bin',mode='r',dtype=np.float64,shape=(ng,7))
        for off in range(0,ng,262144):
            length=min(262144,ng-off)
            assert np.array_equal(gn[off:off+length],original[start+off:start+off+length]), ('gravity nodes',m['rank'],off)
        del gn;node_count+=ng
        rank_rows.append(dict(rank=m['rank'],lo=lo,hi=hi,gravity_nodes=ng,
            bounds_seconds=m['bounds_seconds'],density_seconds=m['density_seconds'],velocity_seconds=m['velocity_seconds']))
    # Bounds may be the historical raw list or the more recent status wrapper.
    for i,row in enumerate(rows[2]):
        expected=baseline_bounds[i]
        if isinstance(expected,dict) and 'bounds' in expected:expected=expected['bounds']
        assert row['bounds']==expected,('density bounds',i)
    for i,(row,expected) in enumerate(zip(rows[3],baseline_fits)):
        for key in ('id','status','prior','parameter','nodes24'):
            assert row[key]==expected[key],('density reconstruction',i,key)
    for i,(row,expected) in enumerate(zip(rows[4],baseline_vel)):
        for key in ('id','status','metadata','diagnostics','generator','gravity','gravity_absolute_sum','moment_defect','points'):
            assert row[key]==expected[key],('velocity preparation',i,key)
    assert node_count==len(original)
    report=dict(status='passed_native_MPI_owned_provider_bitwise',ranks=args.ranks,cells=n,
        global_density_bounds_bitwise=True,density_priors_parameters_bitwise=True,
        velocity_metadata_moments_gravity_bitwise=True,all_gravity_node_values_bitwise=True,
        gravity_nodes=node_count,rank_rows=rank_rows,max_wall_seconds=manifests[0]['max_wall_seconds'],
        validation_seconds=time.perf_counter()-started,provider_calls=0,evolution_steps=0,
        same_host_reference=str(args.same_host_reference) if args.same_host_reference else None,
        limits='CPU preparation component only; not an accepted MPI time evolution or whole solver speedup')
    (run/'comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
    # Delete only this verified component's explicitly enumerated node-local
    # outputs. Failed candidates stay available for a bounded witness audit.
    for m in manifests:
        root=Path(m['buffer_root'])
        for name in ('gravity.bin','context.json','context_2.json','context_3.json','context_4.json',
                     'rows_2.json','rows_3.json','rows_4.json','result_2.json','result_3.json','result_4.json',
                     'global_bounds.json','global_fits.json'):
            path=root/name
            if path.is_file():path.unlink()
        root.rmdir()
    return 0
if __name__=='__main__':raise SystemExit(main())
