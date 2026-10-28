"""Compare native offsets to every saved non-arc nonflat prior in two stages.

Also replays selected Python rational branches and synthetic geometric cases.
No density refit, RHS or evolution. Eight bounded single-thread workers.
"""
import argparse,hashlib,json,multiprocessing as mp,os,resource,sys,time
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
fast=base/'Codex_full_disk_fast_modules_v1_20261019'
old=base/'Codex_full_disk_arc_roots_modules_v7_20261019'
sys.path[:0]=[str(fast),str(old),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from native_exact_prior_bisection import ExactPriorBisection
from curved_disk_cell import DiskCell
from curved_contact_priors import support_extrema,region_area
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False);start=time.perf_counter()
assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1'
available=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
assert len(os.sched_getaffinity(0))>=8 and 8*16*1024**3<available/2
preflight=dict(workers=8,worker_BLAS_OMP_threads=1,worker_address_space_limit=16*1024**3,available_memory=available)
engine=ExactPriorBisection(fast/'libexact_prior.so')

def python_bisect(cell,normal,lo,hi,target):
    witness=0
    for k in range(64):
        mid=(lo+hi)/2
        area=region_area(cell.clip_halfplane(-normal,-mid))
        if area>target:lo=mid;witness|=1<<k
        else:hi=mid
    return (lo+hi)/2,witness

synthetic=[]
square=DiskCell([(F(-1),F(-1)),(F(1),F(-1)),(F(1),F(1)),(F(-1),F(1))],[0,1,2,3],4.,-1)
tiny=F(1,2**70)
thin=DiskCell([(0,0),(1,0),(1,tiny),(0,tiny)],[0,1,2,3],4.,-2)
for cell in (square,thin):
    for normal in (np.array([1.,0.]),np.array([0.,1.]),np.array([.6,.8])):
        lo,hi=support_extrema(cell,normal)
        for fraction in (0.,np.nextafter(0.,1.),.125,.5,np.nextafter(.5,1.),1.):
            target=fraction*region_area(cell)
            expected=python_bisect(cell,normal,lo,hi,target)
            observed=engine.bisect(cell,normal,lo,hi,target)
            assert expected==observed,(cell.owner,normal,target,expected,observed)
            synthetic.append(dict(owner=cell.owner,normal=normal.tolist(),target=target,offset=observed[0],witness=hex(observed[1])))

def init():
    global engine
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    engine=ExactPriorBisection(fast/'libexact_prior.so')

def worker(i):
    meta=metas[i]
    if meta['flat'] or arcs[i]:return dict(index=i,status='unchanged_shape_dispatch',arc=bool(arcs[i]))
    begin=time.perf_counter();cell=DiskCell.from_native(native[i],2.5)
    assert not cell.arcs
    normal=np.asarray(meta['normal']);lo,hi=support_extrema(cell,normal)
    target=meta['fraction']*areas[i]
    before=time.perf_counter();offset,witness=engine.bisect(cell,normal,lo,hi,target);native_time=time.perf_counter()-before
    row=dict(index=i,id=cell.owner,status='passed' if offset==meta['offset'] else 'failed_offset',
        expected=meta['offset'],observed=offset,witness=hex(witness),native_seconds=native_time,total_seconds=time.perf_counter()-begin)
    if i%509==0:
        before=time.perf_counter();oracle=python_bisect(cell,normal,lo,hi,target);row['python_seconds']=time.perf_counter()-before
        if oracle!=(offset,witness):row['status']='failed_branch';row['python_witness']=hex(oracle[1]);row['python_offset']=oracle[0]
    return row

stages=[('initial','Codex_circular_sparse_support_full_20261019','Codex_full_curved_disk_mesh_20261019','Codex_full_disk_initial_density_20261019'),
 ('moving','Codex_full_disk_fine_second_midpoint_geometry_20261019','Codex_full_disk_fine_second_midpoint_geometry_20261019','Codex_full_disk_fine_second_midpoint_density_20261019')]
reports=[]
for label,gn,mn,dn in stages:
    native=json.loads((base/gn/'full_exact_supports.json').read_text())['cells']
    mesh=np.load(base/mn/'mesh.npz');areas=mesh['moments'][:,0];arcs=mesh['arc_counts']
    metas=[json.loads(l)['prior'] for l in (base/dn/'cells.jsonl').read_text().splitlines()]
    begin=time.perf_counter();rows=[]
    with mp.get_context('fork').Pool(8,initializer=init) as pool:
        for row in pool.imap(worker,range(len(metas)),chunksize=16):
            rows.append(row)
            if row['status'].startswith('failed'):
                (args.output/'failure.json').write_text(json.dumps(row,indent=2));raise RuntimeError(row)
    (args.output/(label+'_cells.jsonl')).write_text(''.join(json.dumps(r)+'\n' for r in rows))
    passed=[r for r in rows if r['status']=='passed']
    report=dict(label=label,all_cells=len(rows),native_offsets_bitwise_equal=len(passed),unchanged_dispatch=len(rows)-len(passed),
        python_branch_witness_oracles=sum('python_seconds' in r for r in passed),
        wall_seconds=time.perf_counter()-begin,native_worker_seconds=sum(r['native_seconds'] for r in passed),
        selected_python_worker_seconds=sum(r.get('python_seconds',0) for r in passed))
    reports.append(report);print(json.dumps(report),flush=True)
result=dict(status='passed_saved_native_exact_prior_bisection',synthetic_cases=synthetic,stages=reports,preflight=preflight,
    source_sha256={str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),fast/'native_exact_prior_bisection.cpp',fast/'native_exact_prior_bisection.py',fast/'libexact_prior.so')},
    provider_calls=0,evolution_steps=0,seconds=time.perf_counter()-start,
    limits='Native bisection checked on all saved non-arc nonflat priors in two states. Local stencil/normal construction and full provider still need integration validation. No tolerance, physics, or quadrature-order change.')
(args.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(status=result['status'],seconds=result['seconds'])),flush=True)
