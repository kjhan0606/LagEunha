"""Changed transport law from original t0; preserve the noncorotating trajectory."""
import json,os,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce')
root=b/'Codex_corotating_memory_physical_to_t10_20261018'
proof=json.loads((b/'Codex_corotating_memory_controls_20261018/analysis.json').read_text())
for key in ('off_vs_preserved','restart','retry','parallel_rows'):
    assert all(proof['checks'][key]['fields_bitwise_equal'].values())
    assert proof['checks'][key]['sidecar_bytes_equal']
binary=b/'Codex_corotating_memory_build_20261018/eunha2.exe'
assert sha(binary)==proof['manifest']['binary_sha256']=='36af3ebd200196c7c2eed6dd93a88def07b0909650e972fcc35094fe53af61a3'
source=b/'Codex_canonical_initial_20261017/run/khout.000000.dat'
a,t=read(source);assert t==0
estimate=53*(source.stat().st_size+56+24*len(a));free=shutil.disk_usage(b).free
assert free>2*estimate
print(json.dumps({'target':10.,'dump_interval':.2,'estimated_dump_bytes':estimate,'free_bytes':free}),flush=True)
root.mkdir(exist_ok=False);d=root/'memory_from_t0';d.mkdir()
shutil.copyfile(source,d/'khout.000001.dat');shutil.copyfile(source.parent/'params.dat',d/'params.dat')
template=b/'Codex_reversible_memory_physical_to_t10_20261018/memory_from_t0/environment.json'
case=json.loads(template.read_text())
case.update(GFS_RMEM_COROTATE='1',GFS_RMEM_PARALLEL='1',OMP_NUM_THREADS='4',OMP_STACKSIZE='16M',
            HYDRO_TSTOP='10',GFS_MAX_DT='.0025',EUNHA_DUMP_DT='.2',EUNHA_DUMP_FINAL='1')
assert 'GFS_RMEM_INIT_ZERO' not in case and 'GFS_TEST_REJECT_FINAL_ONCE' not in case
(d/'environment.json').write_text(json.dumps(case,indent=2))
env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(case)
manifest={'target_time':10.,'model':'quadratic residual memory with reconstructed corotational transport',
          'expected_dump_bytes':estimate,'cases':{'memory_from_t0':{'source':str(source),'source_sha256':sha(source),
          'source_step':0,'local_restart_step':1,'start_time':0.,'binary':str(binary),'binary_sha256':sha(binary),
          'dt_cap':.0025,'threads':4,'corotate':True,'parallel_rows':True}},
          'scope':'Changed auxiliary transport model, exact original IC and z0=0. Fixed mass/no steering. Not a stability pass; same environment must be preserved on restart. Existing noncorotating run remains unchanged.'}
with (d/'log').open('x') as stream:
    job=subprocess.Popen(['timeout','86400','mpirun','-np','1',str(binary),'params.dat','1'],
                         cwd=d,env=env,stdout=stream,stderr=subprocess.STDOUT)
    manifest['cases']['memory_from_t0']['pid']=job.pid
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print('corotating_memory_from_t0',job.pid,flush=True)
    code=job.wait()
(d/'exit.code').write_text(str(code));print('exit',code,flush=True);assert code==0
