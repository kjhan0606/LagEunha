"""Static true-arc disk boundary fixtures; no domain evolution or IC fit."""
import argparse,hashlib,json,os,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');modules=base/'Codex_psi_coupled_modules_20261019'
sys.path[:0]=[str(modules),str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_quadratic_thermal_audit_modules_20261018'),str(base/'Codex_quadratic_split_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
os.environ['LAGFORCE_NATIVE_EULER_STAR_LIB']=str(modules/'libnative_euler_star.so')
from circular_power_reference import cells,mean_zero_gauge
from face_record_batch import build
from kinetic_factor_rhs import ALPHA,evaluate_records
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
R=2.5;results={};saved={};eps=np.finfo(float).eps
fixtures=[('single',np.array([[.1,-.2]]),np.zeros(1)),('pair',np.array([[-.6,.1],[.6,-.1]]),np.array([.03,-.02])),('weighted',np.array([[-.8,-.7],[.9,-.6],[.8,.85],[-.7,.9],[0.,0.]]),np.array([.03,-.01,.02,.04,-.02]))]
for name,x,W in fixtures:
 geometry=cells(x,W,R);moments=[c.moments() for c in geometry]
 V=sum(m[0] for m in moments);q=sum((m[1] for m in moments),np.zeros(2));Q=sum((m[2] for m in moments),np.zeros((2,2)))
 expected=np.pi*R**4/4*np.eye(2);bound=4096*eps*max(np.pi*R**4,1.)
 assert abs(V-np.pi*R**2)<=bound and max(abs(q))<=bound and np.max(abs(Q-expected))<=bound
 divergence=[];shared=[];wall=[]
 for i,cell in enumerate(geometry):
  points,normals,ds=cell.quadrature();area=moments[i][0]
  zeroth=ds@normals;first=np.einsum('n,ni,nj->ij',ds,points,normals)-area*np.eye(2)
  assert max(abs(zeroth))<=bound and np.max(abs(first))<=bound
  divergence.append([float(max(abs(zeroth))),float(np.max(abs(first)))])
  for line in cell.lines:
   j=line.neighbor
   if j<i:continue
   other=[line2 for line2 in geometry[j].lines if line2.neighbor==i];assert len(other)==1
   p,n,w=line.quadrature();p2,n2,w2=other[0].quadrature()
   err=max(float(np.max(abs(p-p2[::-1]))),float(np.max(abs(n+n2[::-1]))),float(np.max(abs(w-w2[::-1]))))
   assert err<=bound;shared.append(err)
  wall.extend(a.quadrature() for a in cell.arcs)
 points=np.concatenate([a[0] for a in wall]);normal=np.concatenate([a[1] for a in wall]);ds=np.concatenate([a[2] for a in wall]);nn=len(points)
 tangent=np.column_stack([-normal[:,1],normal[:,0]]);mean=.4*tangent
 samples=mean[:,None,:]+.03*ALPHA[None,:,:]
 cache=dict(points=points,normals=normal,owners=np.column_stack([np.zeros(nn,int),-np.ones(nn,int)]),rho=np.column_stack([np.full(nn,.01),np.full(nn,np.nan)]),means=np.stack([mean,np.full_like(mean,np.nan)],axis=1),samples=np.stack([samples,np.full_like(samples,np.nan)],axis=1),coordinates=np.broadcast_to(ALPHA,(nn,2,4,2)).copy(),record_scale=ds,record_distance=np.zeros(nn))
 states=[dict(P=.001,c=np.zeros(2),v=np.zeros(2),B=np.zeros((2,2)))];gamma=5/3
 records,stats,audit=build(cache,[],states,np.zeros((1,2)),np.zeros((1,2)),np.array([False]),np.zeros((1,2)),np.zeros((1,2)),lambda p:np.zeros(2),gamma,ALPHA)
 residual,J,values=evaluate_records(records,1,np.zeros(1));flux=[];magnitudes=[]
 for rec,(wl,wr,wn) in zip(records,values):
  _,_,scale,_,_,vl,vr,rl,rr,n,p,pl,pr,ps,sm,*_=rec
  fm=sum(wl)+sum(wr);fp=wl@vl+wr@vr+ps*n
  el=.5*np.sum(vl*vl,axis=1)+pl/((gamma-1)*rl);er=.5*np.sum(vr*vr,axis=1)+pr/((gamma-1)*rr)
  fe=wl@el+wr@er+ps*sm;phi=-1/R;torque=p[0]*fp[1]-p[1]*fp[0]
  flux.append(scale*np.array([fm,fe,phi*fm,torque]))
  magnitudes.append(scale*np.array([sum(abs(wl))+sum(abs(wr)),abs(wl)@abs(el)+abs(wr)@abs(er)+abs(ps*sm),abs(phi)*(sum(abs(wl))+sum(abs(wr))),np.linalg.norm(p)*np.linalg.norm(fp)]))
 flux=np.array(flux);magnitudes=np.array(magnitudes);flux_bound=512*eps*np.maximum(magnitudes,ds[:,None]*1e-30)
 assert np.all(abs(flux)<=flux_bound)
 results[name]=dict(status='passed',cells=len(x),arcs=sum(len(c.arcs) for c in geometry),wall_nodes=nn,area_error=float(V-np.pi*R**2),q_error=q.tolist(),Q_error=(Q-expected).tolist(),divergence_errors=divergence,maximum_shared_pair_error=max(shared,default=0.),wall_maximum_absolute_flux=np.max(abs(flux),axis=0).tolist(),wall_total_flux=flux.sum(0).tolist(),wall_flux_order=['mass','hydro_E','Phi','origin_torque'],wall_flux_bounds_sum=flux_bound.sum(0).tolist())
 saved[name+'_wall_points']=points;saved[name+'_wall_normals']=normal;saved[name+'_wall_ds']=ds;saved[name+'_wall_flux']=flux
# Constant W shifts preserve all true cell moments.
x,W=fixtures[-1][1:];a=cells(x,W,R);b=cells(x,W+.125,R)
gauge_error=max(np.max(abs(np.r_[aa.moments()[0],aa.moments()[1],aa.moments()[2].ravel()]-np.r_[bb.moments()[0],bb.moments()[1],bb.moments()[2].ravel()])) for aa,bb in zip(a,b));assert gauge_error<4096*eps*np.pi*R**4
J=np.array([[1.,-1.,0.],[-1.,2.,-1.],[0.,-1.,1.]])
r=np.array([.1,-.3,.2]);z=mean_zero_gauge(J,r);assert max(abs(J@z+r))<64*eps and abs(z.sum())<64*eps
try:mean_zero_gauge(J,np.array([.1,.2,.3]))
except ValueError:pass
else:raise AssertionError('incompatible closed system accepted')
try:cells(np.array([[0.,0.],[1.,0.]]),np.array([-100.,0.]),R)
except ValueError:pass
else:raise AssertionError('empty cell accepted')
results['gauge_and_explicit_failures']=dict(status='passed',moment_gauge_error=float(gauge_error),gauge_solution=z.tolist(),no_projection=True)
np.savez(args.output/'fixtures.npz',**saved)
report=dict(status='passed',scope=__doc__,radius=R,results=results,boundary='stationary true circular slip wall; mirror normal velocity only; ds=R dtheta',future_massroot='all fluid dynamic; explicit sum(z)=0 gauge, reject incompatible total residual; no frozen anchor',future_initialization='integrate steady field on actual disk-clipped cells, including dynamic buffer r=2..2.5; no reuse of old boundary mass',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),modules/'circular_power_reference.py',modules/'face_record_batch.py',modules/'kinetic_factor_rhs.py']},limits='Small all-constraint static reference. No full domain build, interior IC, constrained fit, RHS evolution, or wall-pressure positivity proof.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
