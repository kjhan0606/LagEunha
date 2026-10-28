"""Native/Python state format round trip and initial conservative audit.

No RHS, factor update, accepted step or analytic state replacement.
"""
import argparse,ctypes as C,json,sys
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--library',type=Path,required=True);p.add_argument('--modules',type=Path,required=True);p.add_argument('--geometry',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
sys.path.insert(0,str(args.modules));from full_moment_state_io import load,dump,HEADER
from full_moment_global_audit import admissibility
class Header(C.Structure):
 _fields_=[(k,C.c_uint64) for k in ('magic','version','cells','epoch','run_token','state_elements','ld_bytes','ld_mantissa')]+[('time',C.c_double),('reserved',C.c_uint64*7)]
class State(C.Structure):
 _fields_=[('header',Header)]+[(k,C.c_void_p) for k in ('ids','mass','state','factors','weight_rate_seed','Psi')]+[('mapping',C.c_void_p*6),('bytes',C.c_size_t*6)]
assert C.sizeof(Header)==128
lib=C.CDLL(str(args.library));lib.full_moment_state_load.argtypes=[C.c_char_p,C.POINTER(State)];lib.full_moment_state_load.restype=C.c_int
lib.full_moment_state_write.argtypes=[C.c_char_p,C.POINTER(State)];lib.full_moment_state_write.restype=C.c_int
lib.full_moment_state_publish.argtypes=[C.c_char_p,C.POINTER(State)];lib.full_moment_state_publish.restype=C.c_int
lib.full_moment_state_close.argtypes=[C.POINTER(State)];lib.full_moment_state_close.restype=None
args.output.mkdir(exist_ok=False);base=Path('/gpfs/kjhan/LagForce');ic=base/'Codex_full_curved_kepler_IC_20261019/initial.npz';initial=load(ic)
dump(args.output/'python',initial,1,12345);native=State()
assert lib.full_moment_state_load(str(args.output/'python').encode(),C.byref(native))==0
assert native.header.cells==len(initial['ids']) and native.header.ld_bytes==np.dtype(np.longdouble).itemsize and native.header.ld_mantissa==np.finfo(np.longdouble).nmant+1
assert lib.full_moment_state_write(str(args.output/'native').encode(),C.byref(native))==0
assert lib.full_moment_state_publish(str(args.output/'published').encode(),C.byref(native))==0
assert not (args.output/'published.pending').exists()
published=load(args.output/'published');assert all(np.array_equal(published[k],initial[k]) for k in published)
assert lib.full_moment_state_publish(str(args.output/'published').encode(),C.byref(native))!=0
loaded=load(args.output/'native');exact={key:bool(np.array_equal(loaded[key],initial[key])) for key in loaded};assert all(exact.values()),exact
assert lib.full_moment_state_write(str(args.output/'native').encode(),C.byref(native))!=0
lib.full_moment_state_close(C.byref(native));assert bytes(native)==bytes(C.sizeof(State))
bad=args.output/'bad_header';bad.mkdir();fields=list(HEADER.unpack((args.output/'python/header.bin').read_bytes()));fields[6]=8
(bad/'header.bin').write_bytes(HEADER.pack(*fields));sentinel=State.from_buffer_copy(bytes([0x5A])*C.sizeof(State));before=bytes(sentinel)
status=lib.full_moment_state_load(str(bad).encode(),C.byref(sentinel));assert status and bytes(sentinel)==before
from types import SimpleNamespace
service=SimpleNamespace(ctx=dict(state=str(args.output/'native'),geometry=str(args.geometry),physics_file=str(ic.parent/'report.json')))
request=SimpleNamespace(aux_path=str(ic).encode(),epoch=1,run_token=12345,outer=0,output_dir=str(args.output/'admissibility').encode())
admissibility(service,request)
report=dict(status='passed_native_state_format_and_initial_admissibility',cells=len(initial['ids']),header_bytes=C.sizeof(Header),longdouble_bytes=np.dtype(np.longdouble).itemsize,longdouble_mantissa=np.finfo(np.longdouble).nmant+1,all_fields_exact=exact,durable_atomic_publication_exact=True,overwrite_rejected=True,bad_longdouble_header_no_output_change=True,initial_conservation_audit_passed=True,physical_providers=0,accepted_time_steps=0)
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
