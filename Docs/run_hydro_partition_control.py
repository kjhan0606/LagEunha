"""Validate hydro partition diagnostics without changing the accepted state."""
import json,os,shutil,subprocess,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');old=b/'Codex_memory_parallel_controls_20261018/old'
root=b/'Codex_hydro_partition_control_20261018';binary=b/'Codex_hydro_partition_build_20261018/eunha2.exe'
assert (binary.parent/'build.exit').read_text().strip()=='0'
assert (old/'exit.code').read_text().strip()=='0'
snapshot=old/'khout.000161.dat';side=Path(str(snapshot)+'.rmem')
estimate=2*(snapshot.stat().st_size+side.stat().st_size)
assert shutil.disk_usage(b).free>2*estimate
root.mkdir(exist_ok=False)
for p in (snapshot,side,old/'params.dat'):shutil.copyfile(p,root/p.name)
case=json.loads((old/'environment.json').read_text());case['GFS_TCF_SHARED_ALGEBRA']='1'
(root/'environment.json').write_text(json.dumps(case,indent=2))
manifest={'binary_sha256':sha(binary),'source_sha256':sha(snapshot),'source_sidecar_sha256':sha(side),
    'reference':str(old),'expected_dump_and_source_bytes':estimate}
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))};env.update(case)
with (root/'log').open('x') as stream:
    p=subprocess.Popen(['timeout','900','mpirun','-np','1',str(binary),'params.dat','161'],cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT)
    (root/'pid').write_text(str(p.pid));print('pid',p.pid,flush=True);code=p.wait()
(root/'exit.code').write_text(str(code));assert code==0
left,t=read(old/'khout.000162.dat');right,rt=read(root/'khout.000162.dat');assert t==rt
equal={k:bool(np.array_equal(left[k],right[k])) for k in left.dtype.names}
ztype=[('id','<u8'),('x','<f8'),('y','<f8')]
z1=np.frombuffer((old/'khout.000162.dat.rmem').read_bytes(),dtype=ztype,offset=56)
z2=np.frombuffer((root/'khout.000162.dat.rmem').read_bytes(),dtype=ztype,offset=56)
zequal=bool(np.array_equal(z1,z2))
accepted=[line for line in (root/'log').read_text().splitlines() if line.startswith('[ACCEPTED_STATE]')]
assert len(accepted)==1 and ' finite=1 ' in accepted[0] and ' retries=0 ' in accepted[0]
partition=[line for line in (root/'log').read_text().splitlines() if line.startswith('[ACCEPTED_HYDRO_PARTITION]')]
assert len(partition)==1 and ' phase1=1 ' in partition[0]
result={'partition':partition[0],'time':t,'physical_fields_bitwise_equal':equal,'auxiliary_state_bitwise_equal':zequal,
    'max_abs_difference':{k:float(np.max(abs(left[k]-right[k]))) for k in left.dtype.names if k!='id'},
    'scope':'Hydro partition diagnostic, actual nonzero-memory RK4 step versus preserved old output. State invariance only; not full boundary closure or long-term physics.'}
(root/'analysis.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
assert all(equal.values()) and zequal
