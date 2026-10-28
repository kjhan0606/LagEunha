"""One approved full-disk midpoint factor update using a saved stage.

No provider evaluation; no accepted step. Existing factor map, independent
Python reference and CPU4/CUDA comparison. Stop and preserve any failure.
"""
import argparse,hashlib,json,os,resource,sys,time,traceback
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_compiled_midpoint_modules_20261018'
sys.path[:0]=[str(modules),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from moment_compiled_backend import Backend
from moment_bridge_audit import factor_reference,audit_factor
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--base',type=Path,required=True);parser.add_argument('--linear',type=Path,required=True);parser.add_argument('--stage',type=Path,required=True);parser.add_argument('--h',type=float,required=True);parser.add_argument('--time',type=float,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.perf_counter()
initial_path=args.base;source=args.stage;report=json.loads((source/'report.json').read_text());assert report['status'] in ('passed_actual_full_initial_stage0','passed_actual_full_moving_stage')
initial=np.load(initial_path);linear=np.load(args.linear);saved=np.load(source/'stage0.npz');state=initial['state'];mass=initial['mass'];ids=initial['ids'];L=initial['factors'];Llinear=linear['factors'];assert np.array_equal(initial['ids'],linear['ids']) and np.array_equal(initial['mass'],linear['mass']);N=len(ids);a=state[:-4].reshape(N,20);raw=np.column_stack([mass,a[:,3:]]);frozen=np.zeros(N,dtype=np.int64)
assert np.array_equal(saved['mass'],mass) and np.array_equal(saved['ids'],ids)
assert saved['initial_state_sha256'].item()==hashlib.sha256(np.ascontiguousarray(linear['state']).view(np.uint8)).hexdigest()
st=saved['stage'];Q=saved['incoming_Q'];offset=saved['incoming_offsets'];A=st[:,56:60].reshape(N,2,2);K=st[:,22:26].reshape(N,2,2);gen=saved['generator'];z=saved['weight_rates'];library=base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so'
backends={'cpu4':Backend(library,'cpu',4),'cuda':Backend(library,'cuda',1)}
preflight=dict(input_Q_bytes=Q.nbytes,output_estimate_bytes=2*N*2000,host_work_estimate_bytes=Q.nbytes+4*1024**3,CPU_threads=4,BLAS_threads=os.environ['OPENBLAS_NUM_THREADS'],factor_provider_calls=0,unique_physical_provider_index=json.loads((source/'manifest.json').read_text()).get('unique_physical_provider_index',1),unique_physical_provider_budget=5,source_stage=str(source),source_runner_sha256=json.loads((source/'manifest.json').read_text())['runner_sha256'])
(args.output/'manifest.json').write_text(json.dumps(dict(scope=__doc__,runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),initial_sha256=hashlib.sha256(initial_path.read_bytes()).hexdigest(),linear_path=str(args.linear),linear_sha256=hashlib.sha256(args.linear.read_bytes()).hexdigest(),stage_path=str(args.stage),preflight=preflight),indent=2))
rows=[]
try:
 for name,h in [('candidate',args.h)]:
  assert 0<h<=.001 and args.time>0
  folder=args.output/name;folder.mkdir();provisional=raw+h*st[:,:18];provisional[:,0]=mass
  np.savez(folder/'input.npz',raw=provisional,h=h)
  inputs=(L,Llinear,A,K,Q,offset,h,provisional,frozen);begin=time.perf_counter();reference=factor_reference(*inputs);outputs={k:b.factor(*inputs) for k,b in backends.items()};audit=audit_factor(inputs,reference,outputs)
  out=outputs['cuda'];nextstate=state.copy();b=nextstate[:-4].reshape(N,20);b[:,:2]=a[:,:2]+h*gen;b[:,2]=a[:,2]+h*z;b[:,3:]=out[:,5:22];nextstate[-4:]=state[-4:]+h*saved['ledgers'].reshape(3,4).sum(0);nextL=out[:,:4].reshape(N,2,2)
  nextPsi=initial['Psi'].astype(np.longdouble)+np.longdouble(h)*saved['Psi_rates'].astype(np.longdouble)
  assert np.isfinite(nextstate).all() and np.isfinite(nextPsi).all() and np.array_equal(out[:,4],mass)
  np.savez(folder/'state.npz',state=nextstate,factors=nextL,mass=mass,ids=ids,Psi=nextPsi,time=args.time,frozen=frozen.astype(bool),psi_constrained_mask=np.ones(N,bool),weight_rate_seed=z)
  np.save(folder/'factor_cuda.npy',out);row=dict(label=name,time=args.time,h=h,seconds=time.perf_counter()-begin,minimum_thermal=float(out[:,51].min()),minimum_factor_diagonal=float(out[:,[0,3]].min()),maximum_factor_orthogonality=float(out[:,58].max()),backend_audit=audit,status='passed_predictor')
  (folder/'report.json').write_text(json.dumps(row,indent=2));rows.append(row);print(json.dumps({k:v for k,v in row.items() if k!='backend_audit'}),flush=True)
 result=dict(status='passed_full_disk_factor_candidate',provider_calls=0,accepted_steps=0,predictors=rows,seconds=time.perf_counter()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,limits='Candidate factor update only; actual moving geometry/density and physical endpoint validation pending. Not a time evolution success or long-term result.')
except Exception as exc:
 result=dict(status='failed_full_disk_factor_candidate',provider_calls=0,accepted_steps=0,completed=rows,error=repr(exc),traceback=traceback.format_exc(),seconds=time.perf_counter()-started)
 (args.output/'report.json').write_text(json.dumps(result,indent=2));raise
(args.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(status=result['status'],seconds=result['seconds'])),flush=True)
