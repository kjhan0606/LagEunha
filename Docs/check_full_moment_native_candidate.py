"""Native tentative-state update versus an immutable accepted moment update."""
import ctypes as C,json,sys,hashlib,argparse
from pathlib import Path
import numpy as np
b=Path('/gpfs/kjhan/LagForce');build=b/'Codex_full_disk_mpi_bridge_build_v2_20261019'
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False)
assert C.sizeof(C.c_longdouble)==np.dtype(np.longdouble).itemsize
base=np.load(b/'Codex_full_disk_fine_first_accepted_20261019/accepted.npz')
rates=np.load(b/'Codex_full_disk_fine_second_midpoint_stage_20261019/stage0.npz')
factor=np.load(b/'Codex_full_disk_fine_final_endpoint_candidate_20261019/candidate/factor_cuda.npy')
expected=np.load(b/'Codex_full_disk_fine_final_accepted_20261019/accepted.npz')
n=len(base['mass']);lib=C.CDLL(str(build/'libfull_moment_update.so'))
lib.full_moment_build_candidate.argtypes=[C.c_int,C.c_double]+[C.c_void_p]*11
lib.full_moment_build_candidate.restype=C.c_int
inputs=[np.ascontiguousarray(v) for v in (base['state'][:-4].reshape(n,20),base['mass'],rates['stage'],rates['generator'],rates['weight_rates'],base['Psi'],rates['Psi_rates'],factor)]
ptr=lambda a:a.ctypes.data_as(C.c_void_p)
out=[np.full((n,20),-947.),np.full(n,-947.,dtype=np.longdouble),np.full((n,4),-947.)]
status=lib.full_moment_build_candidate(n,.0005,*map(ptr,inputs),*map(ptr,out));assert status==0,status
checks=dict(state=np.array_equal(out[0],expected['state'][:-4].reshape(n,20)),
    Psi=np.array_equal(out[1],expected['Psi']),factors=np.array_equal(out[2],expected['factors'].reshape(n,4)))
assert all(checks.values()),checks
bad=[a.copy() for a in inputs];bad[-1][n//2,62]=9.
sentinels=[np.full_like(a,-947.) for a in out]
code=lib.full_moment_build_candidate(n,.0005,*map(ptr,bad),*map(ptr,sentinels))
assert code!=0 and all(np.all(a==-947.) for a in sentinels)
report=dict(status='passed_native_candidate_update',cells=n,exact_accepted_array_values=checks,
    longdouble_itemsize=np.dtype(np.longdouble).itemsize,longdouble_mantissa_bits=np.finfo(np.longdouble).nmant+1,
    failure_status=code,failure_caller_outputs_unchanged=True,physical_providers=0,evolution_steps=0,
    source_sha256={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in (build/'full_moment_update.c',build/'full_moment_update.h',build/'libfull_moment_update.so')},
    limits='Tentative native state construction only. Native runtime must still enforce all-rank acceptance, geometry/fit/QR gates and transported global ledger.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
