"""Small closed nonlinear root fixtures; no full fluid provider/evolution."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_circular_curved_modules_v3_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_face_record_native_modules_v2_20261019'),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from closed_weight_root import solve
from kinetic_factor_rhs import ALPHA
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
def face(i,j,speed):
 vl=np.array([speed,0.])+.1*ALPHA;vr=vl.copy();return (i,j,1.,0.,1.,vl,vr,1.,1.2,np.array([1.,0.]),np.array([.5,0.]),.01,.01,.01,0.,0.,ALPHA,ALPHA)
rows={}
for name,n,records in [('pair',2,[face(0,1,.3)]),('chain',3,[face(0,1,.3),face(1,2,-.2)])]:
 z,values,norm,diag=solve(records,n,return_diagnostics=True);assert norm<2e-13 and abs(z.sum())<1e-14
 shifted=solve(records,n,initial=z+.125,return_diagnostics=True)
 assert np.max(abs(shifted[0]-z))<1e-13
 rows[name]=dict(status='passed',z=z.tolist(),norm=norm,diagnostics=diag)
# A stationary mirror wall has no mass leakage, including kinetic samples.
a=list(face(0,-1,.3));a[4]=0.;a[6]=a[5].copy();a[6][:,0]*=-1;a[8]=a[7]
z,values,norm,diag=solve([tuple(a)],1,return_diagnostics=True);assert norm<2e-13 and np.array_equal(z,np.zeros(1));rows['mirror_wall']=dict(status='passed',norm=norm,diagnostics=diag)
for name,n,records in [('disconnected',3,[face(0,1,.3)]),('leaking_wall',1,[face(0,-1,.3)])]:
 try:solve(records,n)
 except RuntimeError as exc:rows[name]=dict(status='passed',explicit_failure=repr(exc),fallback=False)
 else:raise AssertionError(name+' accepted')
report=dict(status='passed',scope=__doc__,results=rows,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),modules/'closed_weight_root.py')},limits='Synthetic graph/mirror tests only; full actual circular domain massroot, dynamic boundary ledgers and two-dt evolution remain unverified.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
