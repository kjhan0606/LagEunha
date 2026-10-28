"""Compare actual production reader fields and auxiliary restart values."""
import json,re,sys,struct
from pathlib import Path
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce')
from run_canonical_av_branches_20261018 import read,sha
b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_reversible_memory_controls_20261018'
manifest=json.loads((root/'manifest.json').read_text());states={};aux={};out={'cases':{},'comparisons':{}}
def sidecar(path,time,ids):
    data=path.read_bytes();header=struct.unpack('<6Q',data[:48]);saved_time=struct.unpack('<d',data[48:56])[0]
    assert header[:3]==(0x524d454d53544154,1,0x0102030405060708)
    assert len(data)==56+24*header[3] and saved_time==time
    records=np.frombuffer(data,dtype=[('id','<u8'),('zx','<f8'),('zy','<f8')],offset=56)
    # Hydro field stores high-bit flags alongside PINDX; sidecar stores PINDX.
    assert np.array_equal(records['id'],ids & np.uint64(72057594037927935))
    assert np.isfinite(records['zx']).all() and np.isfinite(records['zy']).all()
    return records
for name,meta in manifest['cases'].items():
    d=root/name;assert (d/'exit.code').read_text().strip()=='0'
    log=(d/'log').read_text();accepted=[a for a in log.splitlines() if a.startswith('[ACCEPTED_STATE]')]
    assert len(accepted)==meta['expected_steps']
    assert all(' finite=1 ' in line for line in accepted)
    snapshot=d/f"khout.{meta['final_step']:06d}.dat";a,t=read(snapshot);states[name]=(a,t)
    assert sha(Path(meta['source']))==meta['source_sha256']
    if meta['memory_on']:aux[name]=sidecar(Path(str(snapshot)+'.rmem'),t,a['id'])
    else:assert not Path(str(snapshot)+'.rmem').exists()
    out['cases'][name]={'time':t,'snapshot_sha256':sha(snapshot),'accepted_steps':len(accepted),
                        'retries':[int(re.search(r' retries=(\d+)',q).group(1)) for q in accepted],
                        'memory_rows':[q for q in log.splitlines() if q.startswith('[RMEM_ACCEPTED]')]}
first=b/'Codex_reversible_memory_first_step_20261018/khout.000185.dat'
states['omp4']=read(first);aux['omp4']=sidecar(Path(str(first)+'.rmem'),states['omp4'][1],states['omp4'][0]['id'])
all_equal=True
for label,left,right in [('restart','direct_two','resume'),('threads','omp4','omp1'),('retry','omp4','retry'),('default_off','off_old','off_new')]:
    a,t=states[left];c,tc=states[right];assert t==tc
    fields={key:bool(np.array_equal(a[key],c[key])) for key in a.dtype.names}
    maxdiff={key:float(np.max(abs(a[key].astype(float)-c[key].astype(float)))) for key in a.dtype.names}
    row={'fields_bitwise_equal':fields,'field_max_absolute_difference':maxdiff,'time':t}
    all_equal &= all(fields.values())
    if left in aux:
        equal=all(np.array_equal(aux[left][key],aux[right][key]) for key in aux[left].dtype.names)
        row['auxiliary_fields_bitwise_equal']=equal;all_equal &= equal
        row['auxiliary_max_difference']=max(float(np.max(abs(aux[left][key]-aux[right][key]))) for key in ['zx','zy'])
    out['comparisons'][label]=row
assert out['cases']['retry']['retries']==[1]
out['scope']='Actual hydro/gravity/HG/RK4 steps and real sidecar restart. Physical reader fields plus auxiliary values; not whole snapshot byte parity or long-time stability.'
out['all_compared_fields_bitwise_equal']=bool(all_equal)
path=root/'analysis.json'
with path.open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
assert all_equal,'See persisted comparison evidence; do not relaunch unchanged.'
