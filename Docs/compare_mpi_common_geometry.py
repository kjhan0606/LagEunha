from pathlib import Path
import numpy as np,json,sys,re,math
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path("/gpfs/kjhan/LagForce/Codex_mpi_common_geometry_20261017")

def read(n):
 a=np.concatenate([np.loadtxt(root/f'np{n}/raw_faces.dat.common.rank{r:04d}') for r in range(n)])
 assert a.shape==(221159,7) and np.isfinite(a).all() and np.all(a[:,4]>=0)
 order=np.lexsort((a[:,3],a[:,2],a[:,1],a[:,0]));a=a[order];assert len(np.unique(a[:,:4],axis=0))==len(a)
 return a
base=read(1);out={'scope':'Production shf_common after installing received records in diagnostic cache; no force or evolution','ranks':{}}
for n in [1,2,4]:
 a=read(n);assert np.array_equal(a[:,:4],base[:,:4]);diff=np.abs(a[:,4:]-base[:,4:]);assert np.max(diff)==0
 log=(root/f'np{n}/log').read_text();assert log.count('[MPI_COMMON_GEOMETRY]')==n
 out['ranks'][str(n)]={'faces':len(a),'max_abs_length_centroid_difference':np.max(diff,axis=0).tolist(),'zero_length_faces':int(np.count_nonzero(a[:,4]==0)),'missing_reciprocal':sum(map(int,re.findall(r'missing_reciprocal=(\d+)',log)))}
# Independently reconstruct lengths from the prior unmodified single-rank raw dump.
raw=np.loadtxt('/gpfs/kjhan/LagForce/Codex_shf_raw_geometry_20261017/audit/raw_faces.dat');lookup={}
for row in raw:lookup.setdefault((int(row[0]),int(row[1])),[]).append(row)
max_error=0.;missing=0
for row in base:
 i,j=map(int,row[:2]);dx,dy=row[2:4];d=math.hypot(dx,dy);tx,ty=-dy/d,dx/d;mx,my=.5*dx,.5*dy
 own=[g for g in lookup[(i,j)] if g[4]==dx and g[5]==dy];assert len(own)==1
 ax,ay,bx,by=own[0][6:10];s1=(ax-mx)*tx+(ay-my)*ty;s2=(bx-mx)*tx+(by-my)*ty;lo,hi=min(s1,s2),max(s1,s2);found=False
 for g in lookup.get((j,i),[]):
  if abs(g[4]+dx)+abs(g[5]+dy)>1e-9*d:continue
  s3=(g[6]+dx-mx)*tx+(g[7]+dy-my)*ty;s4=(g[8]+dx-mx)*tx+(g[9]+dy-my)*ty;lo=max(lo,min(s3,s4));hi=min(hi,max(s3,s4));found=True;break
 expected=max(hi-lo,0) if found else 0.;missing+=not found;max_error=max(max_error,abs(expected-row[4]))
assert max_error<1e-12
out['independent_reconstruction']={'faces':len(base),'maximum_length_abs_error':max_error,'missing_reciprocal':missing,'note':'Missing reciprocal faces are counted, not silently treated as a physics pass; compare with original serial behavior.'}
(root/'common_comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
