"""Exact small all-plane oracle and full-registry sparse completion count.

The native output is a box-enclosed power polygon awaiting true disk clip.
No fluid projection, RHS, weight solve, or time evolution occurs here.
"""
import argparse,hashlib,json,os,resource,sys,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_circular_sparse_modules_20261019'
sys.path[:0]=[str(modules),str(base)]
from circular_sparse_geometry import SparsePowerGeometry
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--mode',choices=('fixtures','full'),required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
backend=SparsePowerGeometry(modules/'libcircular_sparse.so');started=time.perf_counter();rows={}
def verify_exact(result,xy,W,ids):
 idrow={int(g):i for i,g in enumerate(ids)};sites=[[F(float(a)),F(float(b))] for a,b in xy];weights=list(map(lambda w:F(float(w)),W));checks=0
 for cell in result['cells']:
  i=idrow[cell['owner']]
  for edge in cell['edges']:
   p=list(map(F,edge['A']));a=[p[k]-sites[i][k] for k in (0,1)]
   for j,site in enumerate(sites):
    dx=site[0]-sites[i][0];dy=site[1]-sites[i][1];bound=(dx*dx+dy*dy+weights[i]-weights[j])/2
    assert dx*a[0]+dy*a[1]<=bound;checks+=1
 return checks

def stats(result):
 names=('tree_visits','pruned_nodes','leaf_checks','exact_power_comparisons','filtered_checks','exact_checks','added_candidates')
 return {name:sum(c[name] for c in result['cells']) for name in names}|dict(vertices=sum(len(c['edges']) for c in result['cells']),candidate_count=sum(len(c['candidate_ids']) for c in result['cells']))

if args.mode=='fixtures':
 grid=np.array([(x,y) for y in range(-3,4) for x in range(-3,4)],float)/4
 # Dyadic weights permit exact encoded gauge-shift comparison.
 weights=((np.arange(len(grid))*7)%13).astype(float)/65536
 for name,xy,W,seed in [('single',np.array([[.125,-.25]]),np.zeros(1),0),('pair',np.array([[-.5,.125],[.5,-.125]]),np.array([1.,-1.])/1024,0),('weighted_no_seed',grid,weights,0),('weighted_sparse_seed',grid,weights,12),('weighted_shift',grid,weights+.125,12)]:
  ids=np.arange(len(xy),dtype=np.int64)*3+100
  result=backend.resolve(xy,W,ids,2.5,args.output/(name+'.json'),seed_count=seed,threads=4)
  checks=verify_exact(result,xy,W,ids);rows[name]=dict(status='passed',exact_all_plane_checks=checks,**stats(result))
  if name=='weighted_no_seed':no_seed=result
  if name=='weighted_sparse_seed':sparse=result
  if name=='weighted_shift':shift=result
  print(json.dumps(dict(fixture=name,**rows[name])),flush=True)
 # Exact geometric sets are independent of seed and constant gauge.
 def signature(result):
  return [[(e['neighbor'],tuple(e['A']),tuple(e['B'])) for e in c['edges']] for c in result['cells']]
 # Polygon cyclic starting index depends on construction; compare edge sets.
 def normalized(result):return [sorted(row) for row in signature(result)]
 assert normalized(no_seed)==normalized(sparse)==normalized(shift)
 try:backend.resolve(np.array([[0.,0.],[0.,0.]]),np.zeros(2),np.arange(2),2.5,args.output/'coincident.json',seed_count=1)
 except RuntimeError:rows['coincident_failure']=dict(status='passed',explicit_failure=True)
 else:raise AssertionError('coincident sites accepted')
else:
 from run_canonical_av_branches_20261018 import read
 ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
 allx=np.column_stack([ic['x']-3,ic['y']-3]);ids=np.flatnonzero(np.linalg.norm(allx,axis=1)<2.5).astype(np.int64);xy=allx[ids];W=np.zeros(len(ids))
 meminfo={}
 for line in Path('/proc/meminfo').read_text().splitlines():
  key,*values=line.split()
  if key in ('MemTotal:','MemAvailable:'):meminfo[key[:-1]]=int(values[0])*1024
 # This stage stores sparse exact polygons only. No full volume/sample arrays.
 estimated_bytes=len(ids)*12*128*128
 preflight=dict(cells=len(ids),registry_sha256=hashlib.sha256(np.ascontiguousarray(xy).tobytes()).hexdigest(),ids_sha256=hashlib.sha256(ids.tobytes()).hexdigest(),buffer_cells=int((np.linalg.norm(xy,axis=1)>2).sum()),estimated_geometry_working_bytes=estimated_bytes,memory=meminfo,affinity_cpus=len(os.sched_getaffinity(0)),omp_threads=4,scope='Only native exact support completion; curved-domain quadrature/node count still pending')
 assert estimated_bytes<meminfo['MemAvailable']/4
 (args.output/'resource_preflight.json').write_text(json.dumps(preflight,indent=2));print(json.dumps(preflight),flush=True)
 result=backend.resolve(xy,W,ids,2.5,args.output/'full_exact_supports.json',seed_count=12,threads=4)
 rows['full']=dict(status='completed_exact_supports_not_yet_disk_mesh',preflight=preflight,**stats(result),artifact_bytes=(args.output/'full_exact_supports.json').stat().st_size)
report=dict(status='passed' if args.mode=='fixtures' else 'completed_exact_supports',scope=__doc__,mode=args.mode,seconds=time.perf_counter()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,results=rows,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'sources.json',modules/'libcircular_sparse.so',modules/'circular_sparse_geometry.py']},limits='Exact rational polygon support certification within an enclosing square; true circle intersection, positive quadrature, IC, fluid RHS and long-time stability remain to be connected.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(status=report['status'],seconds=report['seconds'],peak_RSS_KiB=report['peak_RSS_KiB'])))
