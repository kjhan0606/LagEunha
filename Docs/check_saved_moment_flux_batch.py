"""Saved stage0/half accumulation only; no fit, RHS, root or evolution."""
import argparse,ast,hashlib,json,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
basepath=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
modules=basepath/'Codex_moment_flux_batch_modules_20261018'
original=basepath/'Codex_compiled_midpoint_modules_20261018'
sys.path[:0]=[str(modules),str(original),str(basepath/'Codex_coupled_relaxation_evolution_20261018'),str(basepath)]
from vector_moment_flux import accumulate
from kinetic_factor_rhs import ALPHA
source=original/'kinetic_factor_rhs.py'
tree=ast.parse(source.read_text());fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='rhs')
branch=next(x for x in fn.body if isinstance(x,ast.If) and ast.unparse(x.test)=="backend_name() == 'batched'")
prefix=ast.parse('''def scalar(records,values,states,mass,gamma):
 N=len(states)
 dr=np.zeros((N,18));ledger=np.zeros(4);F=np.zeros((N,2));H=np.zeros((N,2,2));K=np.zeros((N,2,2));incoming=[[] for _ in range(N)];area_rates=np.zeros(N)
''').body[0]
prefix.body+=branch.orelse+ast.parse('return dr,ledger,F,H,K,incoming,area_rates').body
base=SimpleNamespace(GAMMA=float(np.float32(1.666667)),affine_velocity=lambda s,p:s['v']+s['B']@(p-s['c']))
exec(compile(ast.fix_missing_locations(ast.Module(body=[prefix],type_ignores=[])),str(source)+'[scalar accumulation only]','exec'))
args.output.mkdir(exist_ok=False)
manifest=dict(scope=__doc__,provider_calls=0,fits=0,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),source,modules/'vector_moment_flux.py')})
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
results={};passed=True
for label in ('stage0','midpoint'):
 start=time.perf_counter();root=basepath/'Codex_compiled_midpoint_once_20261018'/label
 saved=np.load(root/'inputs.npz');ref=np.load(root/'reference.npz');nodes=saved['node_inputs'];owners=saved['owners'];extra=saved['extra'];raw=saved['cells'][:,:18];mass=saved['mass']
 states=[dict(c=r[1:3]/r[0],v=r[7:9]/r[0],B=B) for r,B in zip(raw,ref['B'])]
 x=saved['state'][:-4].reshape(-1,20)[:,:2];gen=saved['generator'];z=saved['weight_rates']
 records=[];values=[]
 for v,(i,j),ext in zip(nodes,owners,extra):
  point=v[:2];normal=v[2:4]
  if j>=0:
   d=np.linalg.norm(x[j]-x[i]);wn0=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/d
   wn=wn0+(z[i]-z[j])/(2*d)
  else:d=0.;wn0=0.;wn=0.
  records.append((i,j,v[4],wn0,d,v[5:13].reshape(4,2),v[13:21].reshape(4,2),1.,1.,normal,point,(base.GAMMA-1)*v[29],(base.GAMMA-1)*v[30],v[31],v[32],0.,ext[:8].reshape(4,2),ext[8:16].reshape(4,2)))
  values.append((v[21:25],v[25:29],wn))
 oracle=scalar(records,values,states,mass,base.GAMMA);candidate=accumulate(records,values,states,mass,base.GAMMA)
 comparisons={}
 for name,k in (('dr',0),('exterior_ledger',1),('F',2),('H',3),('K',4),('area_rates',6)):
  old=oracle[k];new=candidate[k];delta=abs(new-old);bound=1e-14+1e-10*abs(old)
  comparisons[name]=dict(bitwise_equal=new.tobytes()==old.tobytes(),max_absolute_difference=float(delta.max(initial=0)),maximum_tolerance_fraction=float((delta/bound).max(initial=0)))
  passed &= bool(np.all(delta<=bound))
 q=[];qo=[];qexact=True
 for old,new in zip(oracle[5],candidate[5]):
  a=np.concatenate(old,axis=1) if old else np.zeros((2,0));b=np.concatenate(new,axis=1) if new else np.zeros((2,0))
  q.append(b);qo.append(a);qexact &= a.tobytes()==b.tobytes()
 comparisons['ordered_incoming']=dict(bitwise_equal=qexact,sample_columns=sum(v.shape[1] for v in q));passed &= qexact
 # Existing saved CSR oracle verifies incoming ordering, affine conventions,
 # F/H/K and raw face transport independently of the extracted scalar loop.
 csrQ=ref['Q'][:,:8].reshape(-1,4,2);flat=np.concatenate([a.T.reshape(-1,4,2) for a in q])
 delta=abs(flat-csrQ);bound=1e-14+1e-10*abs(csrQ)
 comparisons['saved_CSR_incoming']=dict(maximum_absolute_difference=float(delta.max(initial=0)),maximum_tolerance_fraction=float((delta/bound).max(initial=0)));passed &= bool(np.all(delta<=bound))
 face_target=ref['unmasked_hydro']-ref['source']
 for name,new,old in (('saved_face_dr',candidate[0],face_target),('saved_F',candidate[2],ref['F']),('saved_H',candidate[3],ref['H']),('saved_K',candidate[4],ref['K_unmasked'])):
  delta=abs(new-old);bound=1e-14+1e-10*abs(old)
  comparisons[name]=dict(maximum_absolute_difference=float(delta.max(initial=0)),maximum_tolerance_fraction=float((delta/bound).max(initial=0)));passed &= bool(np.all(delta<=bound))
 np.savez(args.output/(label+'.npz'),dr=candidate[0],ledger=candidate[1],F=candidate[2],H=candidate[3],K=candidate[4],area_rates=candidate[6],ordered_Q=flat)
 results[label]=dict(nodes=len(records),seconds=time.perf_counter()-start,comparisons=comparisons)
 print(json.dumps(results[label]),flush=True)
report=dict(status='passed' if passed else 'numerical_gate_failed',manifest=manifest,stages=results,
 limits='Saved captures lack density and mesh-normal speed. Test uses rho=1 with stored internal energy converted to pressure, and reconstructs wn from saved geometry/rates. Same records feed scalar/candidate; independent CSR comparisons retain existing 1e-14+1e-10 component gate. Ordered incoming scalar/candidate must be bitwise. No provider execution.')
(args.output/'report.json').write_text(json.dumps(report,indent=2))
assert passed,'Preserved failing audit; no fallback or automatic replay'
