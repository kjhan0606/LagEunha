"""Independent original CPU path on the provider's host, before MPI scaling.

An Ada-host capture is not a bitwise oracle for different CPU/libm dispatch.
This check retains exact equality, using the unchanged original preparation
entry points on the same H200 node, plus all independent q48 density gates.
"""
import argparse,hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base/'Codex_full_disk_fast_modules_v1_20261019'))
import full_disk_stage_memory
import full_disk_density_memory as density
import full_disk_velocity_memory as velocity

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    args.output.mkdir(exist_ok=False);start=time.perf_counter()
    manifest=json.loads((args.run/'provider_rank0.json').read_text());root=Path(manifest['buffer_root'])
    assert os.uname().nodename==manifest['host'] and manifest['ranks']==1
    ctx=json.loads((root/'context.json').read_text())
    with np.load(ctx['state']) as f:initial={k:f[k] for k in f.files}
    native=json.loads((Path(ctx['geometry'])/'full_exact_supports.json').read_text())
    with np.load(Path(ctx['geometry'])/'mesh.npz') as f:mesh={k:f[k] for k in f.files}
    fits=density.prepare(initial,native,mesh,ctx['stage_key'],workers=8,independent_q48=True)
    existing=json.loads((root/'rows_3.json').read_text())
    for i,(row,expected) in enumerate(zip(fits.rows,existing)):
        for key in ('id','status','prior','parameter','nodes24'):
            assert row[key]==expected[key],('original density versus spawned provider',i,key)
    (args.output/'density_report.json').write_text(json.dumps(fits.report,indent=2))
    with (args.output/'density_rows.json').open('x') as f:json.dump(fits.rows,f,separators=(',',':'))
    print(json.dumps(dict(status='same_host_original_density_bitwise',independent_q48=True,seconds=time.perf_counter()-start)),flush=True)
    v=velocity.prepare(initial,native,mesh,fits.rows,ctx['stage_key'],workers=8)
    rows=json.loads((root/'rows_4.json').read_text());prepared=v.prepared
    for key in ('generator','gravity','gravity_absolute_sum','moment_defect'):
        assert np.array_equal(prepared[key],np.array([row[key] for row in rows])),key
    for key in rows[0]['metadata']:
        assert np.array_equal(prepared[key],np.array([row['metadata'][key] for row in rows])),key
    ng=len(v.gravity_nodes);gn=np.memmap(root/'gravity.bin',mode='r',dtype=np.float64,shape=(ng,7))
    for off in range(0,ng,262144):assert np.array_equal(v.gravity_nodes[off:off+262144],gn[off:off+262144]),('original gravity nodes versus spawned provider',off)
    np.savez(args.output/'prepared.npz',**prepared)
    (args.output/'velocity_report.json').write_text(json.dumps(v.report,indent=2))
    (args.output/'velocity_rows.json').write_bytes((root/'rows_4.json').read_bytes())
    (args.output/'density_bounds.json').write_bytes((root/'rows_2.json').read_bytes())
    report=dict(status='passed_original_same_host_cpu_reference',host=os.uname().nodename,
        cells=len(initial['ids']),original_density_and_velocity_bitwise=True,independent_q48_all_passed=True,
        original_gravity_nodes_bitwise=True,gravity_nodes=ng,gravity_file=str(root/'gravity.bin'),
        original_sources={Path(module.__file__).name:hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() for module in (density,velocity)},
        seconds=time.perf_counter()-start,evolution_steps=0,provider_calls=0)
    (args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
    np.show_runtime()
    return 0
if __name__=='__main__':raise SystemExit(main())
