"""One real nonzero-memory RK4 step, old/serial/OpenMP1/OpenMP4 comparison."""
import concurrent.futures,json,os,shutil,subprocess,sys,time
from pathlib import Path
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_memory_parallel_controls_20261018'
source=b/'Codex_reversible_memory_physical_to_t10_20261018/memory_from_t0'
snapshot=source/'khout.000161.dat';side=Path(str(snapshot)+'.rmem');a,t=read(snapshot)
assert abs(t-.4)<1e-12 and side.exists()
old=b/'Codex_reversible_memory_integrated_build_20261018/eunha2.exe'
new=b/'Codex_memory_parallel_build_20261018/eunha2.exe'
assert sha(old)=='8614d86c1714151a3ba22602fc8571a8378985ad2ca2a56a1b5e4e0aa061f5d6'
assert (new.parent/'build.exit').read_text().strip()=='0'
estimate=8*(snapshot.stat().st_size+side.stat().st_size);free=shutil.disk_usage(b).free
assert free>2*estimate
root.mkdir(exist_ok=False);jobs=[];manifest={'source':str(snapshot),'source_sha256':sha(snapshot),
    'source_sidecar_sha256':sha(side),'start_time':t,'expected_dump_and_source_bytes':estimate,'free_bytes':free,'cases':{}}
template=json.loads((source/'environment.json').read_text())
for name,binary,parallel,threads in [('old',old,0,4),('new_serial',new,0,4),('parallel1',new,1,1),('parallel4',new,1,4)]:
    d=root/name;d.mkdir()
    for p in (snapshot,side,source/'params.dat'):shutil.copyfile(p,d/p.name)
    case=dict(template);case.update(GFS_RMEM_PARALLEL=str(parallel),OMP_NUM_THREADS=str(threads),
        OMP_STACKSIZE='16M',HYDRO_TSTOP=str(t+.00125),GFS_MAX_DT='.0025',
        EUNHA_DUMP_DT='1000000',EUNHA_DUMP_FINAL='1')
    assert 'GFS_RMEM_INIT_ZERO' not in case and 'GFS_TEST_REJECT_FINAL_ONCE' not in case
    (d/'environment.json').write_text(json.dumps(case,indent=2))
    env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(case)
    manifest['cases'][name]={'binary':str(binary),'binary_sha256':sha(binary),'threads':threads,'parallel_rows':bool(parallel)}
    jobs.append((name,d,binary,env))
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
def execute(job):
    name,d,binary,env=job;start=time.monotonic()
    with (d/'log').open('x') as stream:
        p=subprocess.Popen(['timeout','900','mpirun','-np','1',str(binary),'params.dat','161'],cwd=d,env=env,stdout=stream,stderr=subprocess.STDOUT)
        (d/'pid').write_text(str(p.pid));print(name,p.pid,flush=True);code=p.wait()
    (d/'exit.code').write_text(str(code));(d/'wall_seconds').write_text(str(time.monotonic()-start))
    print(name,'exit',code,flush=True);return code
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:codes=list(pool.map(execute,jobs))
assert all(code==0 for code in codes),codes
