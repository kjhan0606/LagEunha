"""Actual native/kinetic/moment shared-face cancellation on true-disk cells."""
import argparse,hashlib,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_psi_coupled_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ['LAGFORCE_NATIVE_EULER_STAR_LIB']=str(modules/'libnative_euler_star.so')
from circular_power_reference import cells
from face_record_batch import build
from kinetic_factor_rhs import ALPHA,evaluate_records
from vector_moment_flux import accumulate
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
x=np.array([[-.8,-.7],[.9,-.6],[.8,.85],[-.7,.9],[0.,0.]]);W=np.array([.03,-.01,.02,.04,-.02]);N=len(x)
geometry=cells(x,W,2.5);points=[];normals=[];weights=[];owners=[]
for i,c in enumerate(geometry):
 for line in c.lines:
  if i>line.neighbor:continue
  p,n,w=line.quadrature();points.extend(p);normals.extend(n);weights.extend(w);owners.extend([(i,line.neighbor)]*len(w))
p=np.array(points);n=np.array(normals);owners=np.array(owners);ds=np.array(weights);nn=len(p);ids=owners
means=np.stack([np.column_stack([-.4*p[:,1],.4*p[:,0]])+np.array([.01,-.02])*side for side in (0,1)],axis=1)
samples=means[:,:,None,:]+.03*ALPHA[None,None,:,:]
cache=dict(points=p,normals=n,owners=owners,rho=.01*(1+.05*ids),means=means,samples=samples,coordinates=np.broadcast_to(ALPHA,(nn,2,4,2)).copy(),record_scale=ds,record_distance=np.linalg.norm(x[ids[:,1]]-x[ids[:,0]],axis=1))
states=[dict(P=.001*(1+.02*i),c=site.copy(),v=np.array([-.4*site[1],.4*site[0]]),B=np.zeros((2,2))) for i,site in enumerate(x)]
gen=np.array([s['v'] for s in states]);records,stats,_=build(cache,[],states,x,gen,np.zeros(N,bool),x,np.zeros((N,2)),lambda p:np.zeros(2),5/3,ALPHA)
r,J,values=evaluate_records(records,N,np.linspace(-.01,.01,N));dr,ledger,F,H,K,incoming,area=accumulate(records,values,states,np.ones(N),5/3)
roundoff=(128+8*nn)*np.finfo(float).eps
moment_bound=roundoff*np.maximum(abs(dr).sum(0),1e-30);moment_sum=dr.sum(0)
assert np.all(abs(moment_sum)<=moment_bound)
assert np.array_equal(ledger,np.zeros(4))
assert abs(area.sum())<=roundoff*abs(area).sum()
# Potential flux evaluated once per actual shared node and scattered in the
# same alternating owner order, with the very same transport weights.
flux=np.array([rec[2]*(sum(v[0])+sum(v[1]))*(-1/np.linalg.norm(rec[10])) for rec,v in zip(records,values)])
phi_rates=np.zeros(N)
for (i,j),v in zip(owners,flux):phi_rates[i]-=v;phi_rates[j]+=v
assert abs(phi_rates.sum())<=roundoff*abs(flux).sum()
np.savez(args.output/'audit.npz',owners=owners,points=p,normal=n,ds=ds,raw_flux_rates=dr,ledger=ledger,area_rates=area,Phi_flux=flux,Phi_rates=phi_rates)
report=dict(status='passed',scope=__doc__,nodes=nn,moment_sum=moment_sum.tolist(),moment_sum_bound=moment_bound.tolist(),area_sum=float(area.sum()),Phi_rate_sum=float(phi_rates.sum()),provider_calls=0,evolution_steps=0,limits='Static non-equilibrium arbitrary internal face states; no nonlinear mass solve or long-time dynamics.',source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),modules/'circular_power_reference.py',modules/'vector_moment_flux.py']})
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
