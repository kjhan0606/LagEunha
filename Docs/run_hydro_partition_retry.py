"""Actual final-HG rejection versus direct half-step; preserve all outputs."""
import json, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, '/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read, sha
b = Path('/gpfs/kjhan/LagForce')
source = b/'Codex_hydro_partition_control_20261018'
binary = b/'Codex_hydro_partition_build_20261018/eunha2.exe'
root = b/'Codex_hydro_partition_retry_20261018'
snapshot = source/'khout.000161.dat'
side = Path(str(snapshot)+'.rmem')
_, t = read(snapshot)
estimate = 4*(snapshot.stat().st_size+side.stat().st_size)
free = shutil.disk_usage(b).free
assert free > 2*estimate
root.mkdir(exist_ok=False)
(root/'manifest.json').write_text(json.dumps(dict(binary_sha256=sha(binary),
    source_sha256=sha(snapshot), sidecar_sha256=sha(side), start=t,
    estimated_dump_bytes=estimate, free_bytes=free), indent=2))
records = {}
for name in ('direct', 'retry'):
    d = root/name
    d.mkdir()
    for f in (snapshot, side, source/'params.dat'):
        shutil.copyfile(f, d/f.name)
    case = json.loads((source/'environment.json').read_text())
    case.update(HYDRO_TSTOP=str(t+.000625), GFS_MAX_DT='.00125' if name=='direct' else '.0025')
    if name=='retry': case['GFS_TEST_REJECT_FINAL_ONCE']='1'
    (d/'environment.json').write_text(json.dumps(case, indent=2))
    env = {k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))}
    env.update(case)
    with (d/'log').open('x') as f:
        p = subprocess.Popen(['timeout','900','mpirun','-np','1',str(binary),'params.dat','161'],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
        (d/'pid').write_text(str(p.pid)); print(name,p.pid,flush=True)
        code = p.wait()
    (d/'exit.code').write_text(str(code)); assert code==0
    lines=(d/'log').read_text().splitlines()
    def get(tag):
        hits=[s for s in lines if s.startswith('['+tag+']')]
        assert len(hits)==1
        return dict(re.findall(r'(\w+)=([^\s]+)',hits[0]))
    state=get('ACCEPTED_STATE'); partition=get('ACCEPTED_HYDRO_PARTITION')
    assert state['finite']=='1' and int(state['retries'])==(name=='retry')
    assert float(state['dt'])==.00125
    assert sum(s.startswith('[TEST_FINAL_REJECT]') for s in lines)==(name=='retry')
    records[name]=dict(state=state,partition=partition)
left,lt=read(root/'direct/khout.000162.dat'); right,rt=read(root/'retry/khout.000162.dat')
assert lt==rt and abs(lt-(t+.00125))<1e-12
equal={k:bool(np.array_equal(left[k],right[k])) for k in left.dtype.names}
ztype=[('id','<u8'),('x','<f8'),('y','<f8')]
def z(name): return np.frombuffer((root/name/'khout.000162.dat.rmem').read_bytes(),dtype=ztype,offset=56)
zequal=bool(np.array_equal(z('direct'),z('retry')))
pequal=records['direct']['partition']==records['retry']['partition']
result=dict(time=lt,physical_fields_bitwise_equal=equal,auxiliary_state_bitwise_equal=zequal,
    partition_records_identical=pequal,records=records,
    scope='One forced final-HG rejection compared with direct half-step; not long-term stability or full MPI/CUDA evolution.')
(root/'analysis.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
assert all(equal.values()) and zequal and pequal
