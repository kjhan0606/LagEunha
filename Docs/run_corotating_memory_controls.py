"""Separate actual-step controls for opt-in corotating memory, not physics validation."""
import json,os,re,shutil,subprocess,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce')
source=b/'Codex_memory_orbit_work_control_20261018/control'
binary=b/'Codex_corotating_memory_build_20261018/eunha2.exe'
root=b/'Codex_corotating_memory_controls_20261018'
assert (binary.parent/'build.exit').read_text().strip()=='0'
snap=source/'khout.000161.dat';side=Path(str(snap)+'.rmem')
estimate=14*(snap.stat().st_size+side.stat().st_size)
free=shutil.disk_usage(b).free;assert free>2*estimate
root.mkdir(exist_ok=False)
manifest={'binary_sha256':sha(binary),'expected_source_and_final_bytes':estimate,'free_bytes':free,
          'scope':'Deliberate nonzero-z t.4 model branches for implementation controls only. No z reset, no long-time physics.', 'cases':{}}
template=json.loads((source/'environment.json').read_text())
print('expected_bytes',estimate,'free_bytes',free,flush=True)
spec=[('off',0,.0025,1,0,0),('on_full',1,.0025,2,0,0),
      ('on_first',1,.0025,1,0,0),('on_resume',1,.0025,1,0,0),
      ('on_direct_half',1,.00125,1,0,0),('on_retry',1,.0025,1,0,1),
      ('on_parallel4',1,.0025,1,1,0)]
records={}
for name,spin,dt,steps,parallel,reject in spec:
    src=root/'on_first' if name=='on_resume' else source
    n=162 if name=='on_resume' else 161
    p=src/f'khout.{n:06d}.dat';s=Path(str(p)+'.rmem');_,t=read(p)
    effective=.00125 if reject else dt
    d=root/name;d.mkdir()
    for f in (p,s,src/'params.dat'):shutil.copyfile(f,d/f.name)
    case=dict(template)
    case.update(GFS_RMEM_COROTATE=str(spin),GFS_RMEM_PARALLEL=str(parallel),OMP_NUM_THREADS='4',
                OMP_STACKSIZE='16M',GFS_MAX_DT=str(dt),HYDRO_TSTOP=str(t+effective*(steps-.5)),
                EUNHA_DUMP_DT='1000000',EUNHA_DUMP_FINAL='1')
    if reject:case['GFS_TEST_REJECT_FINAL_ONCE']='1'
    assert 'GFS_RMEM_INIT_ZERO' not in case
    (d/'environment.json').write_text(json.dumps(case,indent=2))
    manifest['cases'][name]={'source':str(p),'source_sha256':sha(p),'sidecar_sha256':sha(s),
                             'start_time':t,'steps':steps,'corotate':spin,'parallel_rows':parallel}
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2))
    env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(case)
    with (d/'log').open('x') as log:
        proc=subprocess.Popen(['timeout','900','mpirun','-np','1',str(binary),'params.dat',str(n)],
                              cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
        (d/'pid').write_text(str(proc.pid));print(name,proc.pid,flush=True);code=proc.wait()
    (d/'exit.code').write_text(str(code));assert code==0,(name,code)
    lines=(d/'log').read_text().splitlines()
    accepted=[dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in lines if line.startswith('[ACCEPTED_STATE]')]
    assert len(accepted)==steps
    assert all(a['finite']=='1' and float(a['dt'])==effective for a in accepted)
    assert sum(int(a['retries']) for a in accepted)==reject
    assert sum(line.startswith('[TEST_FINAL_REJECT]') for line in lines)==reject
    spinrows=[dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in lines if line.startswith('[RMEM_SPIN_STAGE]')]
    assert bool(spinrows)==bool(spin)
    if spin:
        assert all(float(a['max_spin'])>0 and float(a['dt_limit'])>0 for a in spinrows)
    final=d/f'khout.{n+steps:06d}.dat';a,tf=read(final)
    assert abs(tf-(t+steps*effective))<1e-12
    records[name]={'final':str(final),'time':tf,'accepted':accepted,'spin_records':spinrows}
    print(name,'exit',code,'time',tf,flush=True)

def compare(left,right):
    left=Path(left);right=Path(right);a,ta=read(left);c,tc=read(right);assert ta==tc
    fields={k:bool(np.array_equal(a[k],c[k])) for k in a.dtype.names}
    sides=Path(str(left)+'.rmem').read_bytes()==Path(str(right)+'.rmem').read_bytes()
    assert all(fields.values()) and sides
    return {'time':ta,'fields_bitwise_equal':fields,'sidecar_bytes_equal':sides}
checks={'off_vs_preserved':compare(records['off']['final'],source/'khout.000162.dat'),
        'restart':compare(records['on_full']['final'],records['on_resume']['final']),
        'retry':compare(records['on_direct_half']['final'],records['on_retry']['final']),
        'parallel_rows':compare(records['on_first']['final'],records['on_parallel4']['final'])}
out={'checks':checks,'records':records,'manifest':manifest,
     'scope':'Implementation controls only. Auxiliary transport remains MPI1 host; same restart environment required. Not a Kepler stability result.'}
(root/'analysis.json').write_text(json.dumps(out,indent=2));print(json.dumps(checks,indent=2),flush=True)
