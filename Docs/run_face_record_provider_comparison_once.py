"""One scalar/native face-record comparison; approved state and accumulation batches remain enabled.

Original saved midpoint, original initial seed and stopping rule. Large-region
timers only; no evolution/update/backend-shadow or automatic repeated replay.
"""
import argparse
import ast
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import os
import resource
import sys
import time
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--mode',choices=('scalar','batch'),required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--preflight-only',action='store_true')
args=parser.parse_args()
original=base/'Codex_compiled_midpoint_modules_20261018'
candidate=base/'Codex_face_record_native_modules_v2_20261019'
os.environ['LAGFORCE_NATIVE_EULER_STAR_LIB']=str(candidate/'libnative_euler_star.so')
paths=[candidate,original]
sys.path[:0]=[str(p) for p in paths]+[str(base/'Codex_quadratic_thermal_audit_modules_20261018'),
    str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
import moment_ale_evolve as model
import fixedmass_hllc
import quadratic_velocity_reconstruction as quad
from kinetic_factor_rhs import rhs as original_rhs,ALPHA
from contact_prior_density import ContactPriorDensity,local_contact_priors
from contact_upwind_face_nodes import ContactUpwindNodes
from moment_gravity_source import source
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from kinetic_face_capture import capture as original_capture
import vector_moment_flux
import native_euler_star
import face_record_batch
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read
flags=dict(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',
    LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',
    LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',
    LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',LAGFORCE_VECTOR_KINETIC_ROOT='0',
    LAGFORCE_SHARED_FACE_VELOCITY='1',LAGFORCE_CAPTURE_PRESSURE_RATES='1',
    LAGFORCE_BATCH_FACE_STATE='1',
    LAGFORCE_VECTOR_MOMENT_FLUX='1',
    LAGFORCE_BATCH_FACE_RECORDS='1' if args.mode=='batch' else '0')
os.environ.update(flags)
assert os.environ['LAGFORCE_BATCH_FACE_STATE']=='1'
assert os.environ.get('OPENBLAS_NUM_THREADS')==os.environ.get('OMP_NUM_THREADS')=='1'
assert not os.environ.get('LAGFORCE_ANGULAR_FACE_CELLS','')
model.GAMMA=float(np.float32(1.666667));fixedmass_hllc.GAMMA=model.GAMMA
saved_root=base/'Codex_compiled_midpoint_once_20261018'
old_manifest=json.loads((saved_root/'manifest.json').read_text())
runner=original/'run_compiled_midpoint_once.py'
assert hashlib.sha256(runner.read_bytes()).hexdigest()==old_manifest['runner_sha256']
for name,meta in old_manifest['modules'].items():
    assert hashlib.sha256(Path(meta['path']).read_bytes()).hexdigest()==meta['sha256']
candidate_manifest=json.loads((candidate/'manifest.json').read_text())
for name,digest in candidate_manifest['python_modules'].items():
    assert hashlib.sha256((candidate/(name+'.py')).read_bytes()).hexdigest()==digest
for name,digest in candidate_manifest['native_files'].items():
    assert hashlib.sha256((candidate/name).read_bytes()).hexdigest()==digest
mid=np.load(saved_root/'midpoint/inputs.npz');stage0=np.load(saved_root/'stage0/inputs.npz')
mass=mid['mass'];ids=mid['ids'];frozen=mid['frozen'].astype(bool);n=len(ids)
y=mid['state'].copy();L=mid['factors'].copy();seed=stage0['weight_rates'].copy()
ic,t0=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t0==0.
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
model.geometry=builder;model.boundary_velocity=lambda point:np.zeros(2)
progress={};events=[];timers={};rhs_calls=[]


@contextmanager
def region(name):
    start=time.perf_counter();cpu=time.process_time()
    try:yield
    finally:
        row=timers.setdefault(name,dict(calls=0,wall_seconds=0.,process_cpu_seconds=0.))
        row['calls']+=1;row['wall_seconds']+=time.perf_counter()-start;row['process_cpu_seconds']+=time.process_time()-cpu


def emit(row):
    events.append(row)
    with (args.output/'events.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(json.dumps(row),flush=True)


def sha_array(array):return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def rhs(*values,**keywords):
    with region('rhs_total'):
        dy,diag=original_rhs(*values,**keywords)
    records=diag['audit_face_records'][0]
    row=dict(nodes=len(records),batch=diag.get('batch_face_state'),record_total=diag['face_record_total_timing'],accumulation=diag['moment_accumulation_timing'],native_records=diag['face_record_batch'])
    if row['batch'] is None:
        row.update(scalar_density_calls=sum(1+(r[1]>=0) for r in records),scalar_sample_calls=sum(1+(r[1]>=0) for r in records))
    rhs_calls.append(row)
    return dy,diag


def capture(*values,**keywords):
    with region('diagnostic_capture'):return original_capture(*values,**keywords)


tree=ast.parse(runner.read_text())
decode_fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='decode')
evaluate_fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='evaluate')
cut=next(k for k,node in enumerate(evaluate_fn.body) if ast.get_source_segment(runner.read_text(),node).startswith('pv=np.array'))
evaluate_fn.body=evaluate_fn.body[:cut]+[ast.parse('return dy,diag,z,gravity,moment_defect,len(records),iteration+1').body[0]]
exec(compile(ast.fix_missing_locations(ast.Module(body=[decode_fn,evaluate_fn],type_ignores=[])),str(runner)+'[unchanged provider prefix]','exec'),globals())
old_decode=decode


def decode(*values,**keywords):
    with region('decode_total'):return old_decode(*values,**keywords)


sources={name:dict(path=str(module.__file__),sha256=hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest())
    for name,module in list(sys.modules.items()) if getattr(module,'__file__',None)
    and Path(module.__file__).is_relative_to(base) and Path(module.__file__).suffix=='.py'}
manifest=dict(scope=__doc__,mode=args.mode,flags=flags,time=float(mid['time']),seed='stage0/inputs.npz:weight_rates',
    seed_sha256=sha_array(seed),state_sha256=sha_array(y),factors_sha256=sha_array(L),runtime_sources=sources,
    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    affinity=sorted(os.sched_getaffinity(0)),load_average_start=os.getloadavg(),
    environment={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')},native_star_provenance=native_euler_star.provenance(),candidate_manifest=candidate_manifest)
if args.preflight_only:
    print(json.dumps(dict(status='preflight_only_no_provider',manifest=manifest),indent=2));raise SystemExit(0)
args.output.mkdir(parents=True,exist_ok=False)
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
with region('provider_total'):
    dy,diag,z,gravity,moment_defect,node_count,iterations=evaluate(float(mid['time']),y,L,seed,'midpoint')
expected=np.load(saved_root/'midpoint/reference.npz')
hydro=np.column_stack([diag['mass_rates'],dy[:-4].reshape(-1,20)[:,3:]])
total=hydro.copy();total[~frozen]+=gravity[~frozen]
old_total=expected['original_hydro'].copy();old_total[~frozen]+=expected['gravity'][~frozen]
comparisons={}
for name,value,old in (('weight_rates',z,mid['weight_rates']),('hydro',hydro,expected['original_hydro']),
                       ('gravity',gravity,expected['gravity']),('raw_source',total,old_total)):
    bound=1e-14+1e-10*abs(old);delta=abs(value-old)
    comparisons[name]=dict(bitwise_equal=value.tobytes()==old.tobytes(),maximum_absolute_difference=float(delta.max()),
        maximum_tolerance_fraction=float((delta/bound).max()))
old_iterations=[json.loads(line) for line in (saved_root/'events.jsonl').read_text().splitlines()]
old_iterations=[x for x in old_iterations if x.get('stage')=='midpoint' and 'iteration' in x]
passed=iterations==len(old_iterations) and diag['mass_flux']<2e-13 and all(x['maximum_tolerance_fraction']<=1 for x in comparisons.values())
np.savez(args.output/'replayed.npz',state=y,factors=L,seed=seed,weight_rates=z,hydro=hydro,gravity=gravity,raw_source=total)
report=dict(status='passed' if passed else 'numerical_gate_failed',manifest=manifest,
    iterations=iterations,original_iterations=len(old_iterations),mass_residual=diag['mass_flux'],comparisons=comparisons,
    timers=timers,rhs_calls=rhs_calls,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.,
    load_average_end=os.getloadavg(),
    excluding_face_export={key:timers['provider_total'][key]-timers['diagnostic_capture'][key] for key in ('wall_seconds','process_cpu_seconds')},
    limits='One provider call per mode, no update. Large inclusive regions overlap and must not be summed. Peak RSS includes imports/input/density/diagnostic data. Shared-host wall times are not private CPU allocation. Scalar records retain their original loop: record_total includes scalar star and cannot be split without perturbing that oracle. Native array preparation/star/packing are separate large regions. excluding_face_export subtracts only capture(), retaining gravity/moment-defect audit overhead; it is not a production-step timing.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
assert passed,'Preserved failed numerical result; no fallback or repeat'
