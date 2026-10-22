"""Replay exactly the previously passing four native quarter steps once."""
from pathlib import Path
import hashlib, json, os, subprocess, sys
import numpy as np

root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
base=root/'rk4_diagnostics/fine_cold_interval'
manifest=json.loads((base/'manifest.json').read_text())
ledger=json.loads((base/'ledger/results.json').read_text())
dest=root/'rk4_diagnostics/branch_trace/quarter_native';dest.parent.mkdir(exist_ok=True);dest.mkdir(exist_ok=False)
(dest/'OWNER.txt').write_text('LagForce Codex native branch diagnostic. Same captured interval, four quarter steps. Preserve outputs.\n')
binary=root/'raw_vertex_dedup/v6/eunha2_raw_vertex_dedup.exe'
assert binary.is_file() and (binary.parent/'build_status.txt').read_text().strip()=='BUILD_COMPLETE'
source=Path(manifest['source']);old=base/'quarters'
(dest/'khout.000001.dat').symlink_to(source)
(dest/'params.dat').symlink_to(old/'params.dat')
(dest/'eunha2.exe').symlink_to(binary)
env=os.environ.copy();env.update(json.loads((old/'environment.json').read_text()))
env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='1',I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
    GFS_GPU_FACE_DUMP_CALL='0',GFS_GPU_FACE_DUMP_FIRST='1',GFS_GPU_FACE_DUMP_LAST='16',
    GFS_GPU_FACE_DUMP_IDS=','.join(map(str,ledger['selected_ids'])),GFS_GPU_TESS_DUMP='0')
for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT'):env.pop(key,None)
(dest/'environment.json').write_text(json.dumps({k:v for k,v in env.items() if k.startswith(('GFS_','HYDRO_','KEPLER_','EUNHA_'))},indent=2)+'\n')
code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=dest,env=env).returncode
row={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
     'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
     'exit_code':code,'scope':'instrumented replay, same physics and same four quarter steps, not a long validation',
     'output_policy':'four native checkpoints about156MiB plus selected face/input TSVs; existing outputs preserved'}
if code==0:
    new_files=sorted(f for f in dest.glob('khout.*.dat') if f.name!='khout.000001.dat')
    old_files=sorted(f for f in old.glob('khout.*.dat') if f.name!='khout.000001.dat')
    assert len(new_files)==len(old_files)==4
    comparisons=[]
    for a,b in zip(new_files,old_files):
        p,ma=read_snapshot(a);q,mb=read_snapshot(b)
        ia=np.argsort(p['u4if']);ib=np.argsort(q['u4if'])
        assert np.array_equal(p['u4if'][ia],q['u4if'][ib])
        errors={k:float(abs(p[k][ia].astype(float)-q[k][ib].astype(float)).max()) for k in p.dtype.names if k!='u4if'}
        comparisons.append({'new':str(a),'reference':str(b),'time':ma['t'],
            'same_metadata':ma==mb,'max_absolute_errors':errors,'all_fields_bitwise_equal':all(e==0 for e in errors.values())})
    row['comparisons']=comparisons
    row['same_endpoint']=comparisons[-1]['time']==manifest['target_time']
    row['instrumentation_equivalent']=all(c['same_metadata'] and c['all_fields_bitwise_equal'] for c in comparisons)
    row['trace_calls']=len(list(dest.glob('face_input_meta_r0_call*.json')))
(dest/'manifest.json').write_text(json.dumps(row,indent=2)+'\n')
print(json.dumps(row,indent=2),flush=True)
if code or not row.get('instrumentation_equivalent') or row.get('trace_calls')!=16:
    raise RuntimeError('stop: replay/instrumentation equivalence failed')
