"""New physical method from preserved original t0 through historical t6-10."""
import json,os,shutil,subprocess,sys
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_reversible_memory_physical_to_t10_20261018'
proof=json.loads((b/'Codex_reversible_memory_controls_20261018/analysis.json').read_text())
assert proof['all_compared_fields_bitwise_equal']
source=b/'Codex_canonical_initial_20261017/run/khout.000000.dat';a,t=read(source);assert t==0
binary=b/'Codex_reversible_memory_integrated_build_20261018/eunha2.exe'
assert sha(binary)=='8614d86c1714151a3ba22602fc8571a8378985ad2ca2a56a1b5e4e0aa061f5d6'
estimate=53*(source.stat().st_size+56+24*len(a));free=shutil.disk_usage(b).free;assert free>2*estimate
print(json.dumps({'target':10.,'dump_interval':.2,'estimated_dump_bytes':estimate,'free_bytes':free}),flush=True)
root.mkdir(exist_ok=False);d=root/'memory_from_t0';d.mkdir()
# icont=0 regenerates IC in KH. Rename an exact copy to index1 for true restart.
shutil.copyfile(source,d/'khout.000001.dat')
shutil.copyfile(source.parent/'params.dat',d/'params.dat')
case=json.loads((b/'Codex_reversible_memory_first_step_20261018/environment.json').read_text())
case.update(HYDRO_TSTOP='10',GFS_MAX_DT='.0025',EUNHA_DUMP_DT='.2',EUNHA_DUMP_FINAL='1',OMP_NUM_THREADS='4')
case.pop('GFS_RMEM_INIT_ZERO',None) # time zero is explicit in loaded snapshot
assert 'GFS_TEST_REJECT_FINAL_ONCE' not in case
(d/'environment.json').write_text(json.dumps(case,indent=2))
env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(case)
with (d/'log').open('x') as stream:
    job=subprocess.Popen(['timeout','86400','mpirun','-np','1',str(binary),'params.dat','1'],cwd=d,env=env,stdout=stream,stderr=subprocess.STDOUT)
    manifest={'target_time':10.,'cases':{'memory_from_t0':{'source':str(source),'source_sha256':sha(source),'source_step':0,'local_restart_step':1,'start_time':0.,'pid':job.pid,'binary':str(binary),'binary_sha256':sha(binary),'dt_cap':.0025,'threads':4}},'scope':'New auxiliary-energy method starts from exact original t0. Existing historical runs and checkpoints preserved. Physical failure interval and later long-time gates remain mandatory.'}
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2));print('memory_from_t0',job.pid,flush=True)
    code=job.wait()
(d/'exit.code').write_text(str(code));print('memory_from_t0 exit',code,flush=True)
assert code==0
