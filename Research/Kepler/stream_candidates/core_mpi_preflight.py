"""One native fine-core step, four MPI ranks versus retained one-rank state."""
from pathlib import Path
import json,os,subprocess,sys
import numpy as np
root=Path(__file__).resolve().parent;sys.path.insert(0,str(root))
from audit_snapshot import read_snapshot
base=root/'core_refinement';d=base/'mpi4_preflight_v4';d.mkdir(exist_ok=False)
source=base/'fine_integrated';reference=base/'native_v4/preflight_fine'
for name in ('params.dat','khout.000001.dat'):(d/name).symlink_to(source/name)
(d/'eunha2.exe').symlink_to(root/'raw_vertex_dedup/v4/eunha2_raw_vertex_dedup.exe')
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',NATIVE_MPI_RANKS='4',I_MPI_PIN='0',I_MPI_HYDRA_BOOTSTRAP='fork',
    KEPLER_RESTART_STEP='1',EUNHA_IC='kepler',HYDRO_TSTOP='1.e-7',EUNHA_DUMP_DT='1.e-7',GFS_EXACT_STOP='1',
    GFS_FAIL_NONPOSITIVE_U='1',GFS_RK4_FAILURE_SNAPSHOT='0',GFS_RK4_AUDIT_STEP_MIN='1',GFS_RK4_AUDIT_STEP_MAX='1',
    GFS_RK4_AUDIT_STAGE='k1_rhs',GFS_GPU_FACE_DUMP_CALL='1',GFS_GPU_TESS_DUMP='0',GFS_GPU_FINAL_VOLUME='1',
    GFS_GPU_RAW_FACES='1',GFS_VELOCITY_LIMITER='0')
for key in ('I_MPI_PMI_LIBRARY','PMI_RANK','PMI_SIZE','PMI_FD','PMI_PORT','GFS_DIAG_DT_MAX'):env.pop(key,None)
code=subprocess.run([sys.executable,str(root/'run_native.py')],cwd=d,env=env).returncode
if code:raise RuntimeError(f'fine MPI4 preflight failed {code}')
subprocess.run([sys.executable,str(root/'audit_face_reciprocity_v2.py'),'--output',str(d/'reciprocity.json'),str(d)],check=True)
a,ma=read_snapshot(reference/'khout.000002.dat');b,mb=read_snapshot(d/'khout.000002.dat')
order=lambda p:np.argsort(p['u4if']&np.uint64((1<<56)-1))
a=a[order(a)];b=b[order(b)]
assert np.array_equal(a['u4if'],b['u4if']) and len(a)==71080
errors={field:float(np.max(abs(a[field].astype(float)-b[field].astype(float)))) for field in a.dtype.names if field!='u4if'}
r={'scope':'one step only; no MPI scaling or long-term stability claim','ranks':4,'exit_code':code,
   'meta':mb,'reference_meta':ma,'all_recorded_fields_bitwise_equal':not any(errors.values()),'state_max_absolute_error':errors}
(d/'comparison.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2),flush=True)
