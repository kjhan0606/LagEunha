"""Full true-disk geometry and positive quadrature counts, without fluid RHS.

Four bounded workers, one cell's volume nodes alive per worker. Only compact
cell moments, canonical shared faces and counts are retained, no node dump.
"""
import argparse,hashlib,json,multiprocessing as mp,os,resource,sys,time,traceback
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_circular_curved_modules_v3_20261019';sys.path[:0]=[str(modules)]
from curved_disk_cell import DiskCell,ld,sub,dot,sqrt_enclosure
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
args.output.mkdir(exist_ok=False);started=time.perf_counter();source=base/'Codex_circular_sparse_support_full_20261019/full_exact_supports.json';raw=json.loads(source.read_text());native_cells=raw['cells'];N=len(native_cells)
assert N==19945 and raw['shared_exact_reversed']
assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1' and os.environ['OMP_PROC_BIND']=='false'
input_data=np.load(source.with_suffix('.inputs.npz'));ids=input_data['ids'];sites=input_data['sites'];R=2.5
meminfo={}
for line in Path('/proc/meminfo').read_text().splitlines():
 key,*values=line.split()
 if key in ('MemTotal:','MemAvailable:'):meminfo[key[:-1]]=int(values[0])*1024
preflight=dict(cells=N,workers=4,worker_BLAS_OMP_threads=1,allowed_cpus=sorted(os.sched_getaffinity(0)),memory=meminfo,source_bytes=source.stat().st_size,estimated_additional_output_bytes=N*1000,volume_nodes_retained='one cell per worker, no global node arrays',worker_virtual_memory_limit_bytes=16*1024**3)
assert len(preflight['allowed_cpus'])>=4 and 64*1024**3<meminfo['MemAvailable']/2
(args.output/'resource_preflight.json').write_text(json.dumps(preflight,indent=2));print(json.dumps(preflight),flush=True)

def worker_init():
 resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))

def worker(i):
 try:
  c=DiskCell.from_native(native_cells[i],R);V,q,Q=c.moments();origin=c.origin
  p,w=c.positive_quadrature(24,(.25,.5,2.));pl=p.astype(np.longdouble);wl=w.astype(np.longdouble)
  observed=np.r_[wl.sum(),wl@pl,np.einsum('n,ni,nj->ij',wl,pl,pl).ravel()];target=np.r_[V,q,Q.ravel()]
  error=observed-target;gate=bool(np.all(abs(error)<=1e-14+1e-10*abs(target)))
  face=[];shared_sigs=[]
  for line in c.lines:
   if line.neighbor<0:raise RuntimeError('positive box face survived disk intersection')
   endpoints=(line.A,line.B) if c.owner<line.neighbor else (line.B,line.A)
   digest=hashlib.sha256('|'.join(str(q) for point in endpoints for q in point).encode()).hexdigest()
   shared_sigs.append((min(c.owner,line.neighbor),max(c.owner,line.neighbor),digest))
   if c.owner>line.neighbor:continue
   delta=sub(line.B,line.A);lo,hi=sqrt_enclosure(dot(delta,delta));length=ld((lo+hi)/2)
   mid=np.array([ld((a+b)/2) for a,b in zip(line.A,line.B)]);d=np.array([ld(v) for v in delta]);normal=np.array([d[1],-d[0]])/length
   face.append((c.owner,line.neighbor,*map(float,mid),*map(float,d),float(length),*map(float,normal)))
  return dict(id=c.owner,index=i,status='passed' if gate else 'positive_quadrature_gate_failed',moments=np.asarray(target,dtype=float).tolist(),quadrature_errors=np.asarray(error,dtype=float).tolist(),nodes=len(p),arcs=len(c.arcs),lines=len(c.lines),faces=face,shared_sigs=shared_sigs,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except Exception as exc:return dict(id=native_cells[i]['owner'],index=i,status='failed',error=repr(exc),traceback=traceback.format_exc())

rows=[];faces=[];counts={};max_rss=0;shared_hashes={};shared_counts={}
with mp.get_context('fork').Pool(4,initializer=worker_init) as pool:
 for row in pool.imap(worker,range(N),chunksize=8):
  counts[row['status']]=counts.get(row['status'],0)+1
  for a,b,digest in row.pop('shared_sigs',[]):
   key=(a,b)
   if key in shared_hashes and shared_hashes[key]!=digest:
    witness=dict(pair=key,first=shared_hashes[key],second=digest)
    (args.output/'shared_pair_failure.json').write_text(json.dumps(witness));raise RuntimeError(witness)
   shared_hashes[key]=digest;shared_counts[key]=shared_counts.get(key,0)+1
  faces.extend(row.pop('faces',[]));max_rss=max(max_rss,row.get('peak_RSS_KiB',0));rows.append(row)
  with (args.output/'cells.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
  if len(rows)%500==0:print(json.dumps(dict(cells=len(rows),counts=counts,seconds=time.perf_counter()-started)),flush=True)
assert [r['id'] for r in rows]==ids.tolist()
failures=[r for r in rows if r['status']!='passed']
if failures:
 report=dict(status='failed',counts=counts,failures=failures,seconds=time.perf_counter()-started)
 (args.output/'report.json').write_text(json.dumps(report,indent=2));raise RuntimeError(dict(failed_cells=len(failures),first=failures[0]))
assert all(count==2 for count in shared_counts.values())
faces=np.array(faces);keys=faces[:,:2].astype(np.int64);assert len(np.unique(keys,axis=0))==len(keys)
moments=np.array([r['moments'] for r in rows]);errors=np.array([r['quadrature_errors'] for r in rows]);nodes=np.array([r['nodes'] for r in rows]);arcs=np.array([r['arcs'] for r in rows]);lines=np.array([r['lines'] for r in rows]);assert lines.sum()==2*len(faces)
exact=np.r_[np.pi*R**2,0.,0.,np.pi*R**4/4,0.,0.,np.pi*R**4/4];total=moments.astype(np.longdouble).sum(0);partition_error=total-exact
bound=(128+4*N)*np.finfo(float).eps*np.maximum(abs(moments).sum(0),1.)
assert np.all(abs(partition_error)<=bound)
np.savez(args.output/'mesh.npz',ids=ids,sites=sites,weights=np.zeros(N),radius=R,moments=moments,quad24_errors=errors,quad24_counts=nodes,arc_counts=arcs,line_counts=lines,face_keys=keys,faces=faces[:,2:])
report=dict(status='passed_full_curved_geometry_not_fluid_IC',scope=__doc__,cells=N,canonical_shared_faces=len(faces),shared_circle_clipped_endpoint_representatives_exact_reversed=True,wall_cells=int((arcs>0).sum()),wall_arcs=int(arcs.sum()),positive_quad24_nodes=int(nodes.sum()),maximum_cell_quad24_nodes=int(nodes.max()),quadrature_maximum_component_errors=abs(errors).max(0).tolist(),partition_moments_error=np.asarray(partition_error,dtype=float).tolist(),partition_moments_roundoff_bound=bound.tolist(),worker_maximum_RSS_KiB=max_rss,parent_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,seconds=time.perf_counter()-started,preflight=preflight,output_bytes=sum((args.output/f).stat().st_size for f in ('mesh.npz','cells.jsonl','resource_preflight.json')),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),modules/'curved_disk_cell.py',source)},limits='Geometry and volume quadrature only, order24 on each cell. No conservative fluid IC, independent order48 IC, closed root, actual dynamic wall flux, or time evolution. Full cell-level checks and failures preserved; no full volume nodes stored.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k not in ('preflight','source_sha256')},indent=2))
