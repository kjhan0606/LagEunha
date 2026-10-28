"""Hybrid positive density with affine mean velocity and retained residual covariance."""
from pathlib import Path
import ast,importlib.util,json,time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
base=Path('/gpfs/kjhan/LagForce');out=base/'Codex_hybrid_affine_path_20261018';out.mkdir(exist_ok=False)
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
# Reuse only pure function definitions, without executing the diagnostic mains.
for filename,names in [('audit_quadratic_density_reconstruction.py',('feature','minimum')),('audit_remap_true_spatial.py',('recover',))]:
 tree=ast.parse((base/filename).read_text());functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names];assert len(functions)==len(names);exec(compile(ast.Module(body=functions,type_ignores=[]),filename,'exec'),globals())
mod=load('mom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geo=load('geo',base/'Codex_weighted_compression_geometry_20261017/audit.py');rd=load('rd',base/'run_canonical_av_branches_20261018.py')
a,t=rd.read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat');x=np.column_stack((a['x']-3,a['y']-3));ids=np.flatnonzero(np.linalg.norm(x,axis=1)<.7);m=a['mass'][ids];vel=np.column_stack((a['vx'][ids],a['vy'][ids]));thermal=a['ie'][ids];old=mod.polygons(x,np.zeros(len(x)),ids,geo.clip);target=np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
def integrate(xy,dv,rho,v,Rspec,ud):
 dm=dv*rho;M=dm.sum();sp=dm@xy;Q=np.einsum('n,ni,nj->ij',dm,xy,xy);D=np.einsum('n,ni,nj->ij',dm,xy,v);Pi=np.einsum('n,ni,nj->ij',dm,v,v)+M*Rspec
 return np.r_[M,sp,Q[0,0],Q[0,1],Q[1,1],dm@v,D.ravel(),Pi[0,0],Pi[0,1],Pi[1,1],ud*dv.sum()+.5*np.trace(Pi),0.]
q0=[]
for k,p in enumerate(old):
 xy,dv=mod.quadrature(p,4);V=dv.sum();q0.append(integrate(xy,dv,np.full(len(xy),m[k]/V),np.broadcast_to(vel[k],xy.shape),np.zeros((2,2)),thermal[k]/V))
q0=np.array(q0)
def fit_positive(p,origin,h,desired):
 V=mod.area(p);theta=np.zeros(5)
 def grid(order):
  xy,dv=mod.quadrature(p,order);return feature((xy-origin)/h)[:,1:],np.log(dv/V)
 def stats(b,F,lw):
  l=lw+F@b;Z=logsumexp(l);prob=np.exp(l-Z);mu=prob@F;d=F-mu
  return mu,(d.T*prob)@d,Z
 for order,check in ((16,32),(32,48),(48,64)):
  F,lw=grid(order);sol=least_squares(lambda b:stats(b,F,lw)[0]-desired,theta,jac=lambda b:stats(b,F,lw)[1],xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=200);theta=sol.x;mu,_,Z=stats(theta,F,lw);Fc,lwc=grid(check);muc,_,Zc=stats(theta,Fc,lwc);ratio=np.exp(Zc-Z);err=float(abs(np.r_[ratio,ratio*muc]-np.r_[1.,desired]).max())
  if err<1e-11:
   b=np.r_[-Z,theta];low=np.exp(minimum((p-origin)/h,b));return b,float(low),order,err
 raise RuntimeError('Positive fit independent quadrature gate failed')
def reconstruct(polys,q):
 S,C,G,R=recover(q);mass=q[:,0];c=q[:,1:3]/mass[:,None];v=q[:,6:8]/mass[:,None];U=q[:,15]-.5*(q[:,12]+q[:,14]);ev=np.linalg.eigvalsh(R)[:,0];tol=1e-12*np.max(abs(q[:,12:15]),axis=1)+1e-20;fail=[];data=[];lowest=float('inf')
 for k,p in enumerate(polys):
  xy,dv=mod.quadrature(p,4);V=dv.sum();h=np.sqrt(V);origin=x[ids[k]];F=feature((xy-origin)/h);gram=F.T@(dv[:,None]*F)/V
  Q=np.array([[q[k,3],q[k,4]],[q[k,4],q[k,5]]])/mass[k];mean=(c[k]-origin)/h;second=(Q-np.outer(origin,c[k])-np.outer(c[k],origin)+np.outer(origin,origin))/h**2;rhs=np.array([1,*mean,second[0,0],second[0,1],second[1,1]]);b=np.linalg.solve(gram,rhs);low=minimum((p-origin)/h,b);exporder=0;fiterr=0.
  if low<=0:b,low,exporder,fiterr=fit_positive(p,origin,h,rhs[1:])
  lowest=min(lowest,low)
  if low<=0 or ev[k]<-tol[k] or U[k]<=0:fail.append({'id':int(ids[k]),'min_density_over_mean':low,'min_R_eigenvalue':float(ev[k]),'R_tolerance':float(tol[k]),'U':float(U[k])})
  data.append((V,h,origin,b,mass[k]/V,c[k],v[k],G[k],R[k]/mass[k],U[k]/V,exporder,fiterr))
 return data,fail,lowest,float(ev.min())
def advance(source,dest,q,data):
 lo=np.array([p.min(0) for p in source]);hi=np.array([p.max(0) for p in source]);post=np.zeros_like(q);covered=np.zeros(len(ids));newcovered=np.zeros(len(ids))
 for j,p in enumerate(dest):
  cand=np.flatnonzero(np.all(hi>=p.min(0),axis=1)&np.all(lo<=p.max(0),axis=1))
  for i in cand:
   poly=source[i].copy()
   for aa,bb in zip(p,np.roll(p,-1,axis=0)):
    e=bb-aa;n=np.array([e[1],-e[0]]);poly=geo.clip(poly,n,n@aa)
    if len(poly)<3:break
   if mod.area(poly)<=1e-18:continue
   V,h,origin,b,rhobar,c,v,G,res,ud,exporder,fiterr=data[i]
   def piece(order):
    xy,dv=mod.quadrature(poly,order);f=feature((xy-origin)/h)@b;rho=rhobar*(np.exp(f) if exporder else f);vv=v+(xy-c)@G.T
    return integrate(xy,dv,rho,vv,res,ud),dv.sum()
   contribution,av=piece(max(4,exporder))
   if exporder:
    higher,av2=piece(exporder+16);assert np.max(abs(higher-contribution))<1e-12;contribution=higher
   post[j]+=contribution;covered[i]+=av;newcovered[j]+=av
 err=max(float(abs(covered/np.array([d[0] for d in data])-1).max()),float(abs(newcovered/np.array([mod.area(p) for p in dest])-1).max()));assert err<1e-8
 assert np.max(abs(post.sum(0)-q.sum(0)))<1e-10
 return post,err
rows=[]
for steps in (4,8):
 start=time.monotonic();q=q0.copy();source=old;history=[];passed=True;path=list(np.linspace(0,.25,steps+1)[1:])+list(np.linspace(.25,0,steps+1)[1:])
 for leg in range(len(path)+1):
  data,fail,low,eig=reconstruct(source,q);h={'state_index':leg,'min_density_over_mean':low,'min_R_eigenvalue':eig,'failures':fail,'exponential_cells':sum(d[-2]>0 for d in data),'max_fit_check_error':max(d[-1] for d in data)};history.append(h)
  if fail:passed=False;break
  if leg==len(path):break
  dest=mod.polygons(x,path[leg]*target,ids,geo.clip);q,err=advance(source,dest,q,data);source=dest;h['next_coverage_error']=err
  print(json.dumps({'N':steps,'completed_transition':leg+1,'weight':float(path[leg])}),flush=True)
 row={'subdivisions':steps,'passed_full_path':passed,'history':history,'total_moment_drift':(q.sum(0)-q0.sum(0)).tolist(),'seconds':time.monotonic()-start}
 if passed:
  v=q[:,6:8]/q[:,0,None];internal=.5*(q[:,12]+q[:,14])-.5*q[:,0]*np.sum(v*v,axis=1);row.update(internal_K_over_initial_U=float(internal.sum()/thermal.sum()),velocity_change_mass_RMS=float(np.sqrt(np.average(np.sum((v-vel)**2,axis=1),weights=m))),represented_L_change=float(np.sum(q[:,9]-q[:,10])-np.sum(q0[:,9]-q0[:,10])))
 rows.append(row);np.savez(out/f'steps{steps}.npz',ids=ids,state=q)
 print(json.dumps({'N':steps,'passed':passed,'last_state':history[-1]['state_index'],'failed_cells':len(history[-1]['failures']),'min_density':history[-1]['min_density_over_mean'],'min_R_eig':history[-1]['min_R_eigenvalue']}),flush=True)
(out/'result.json').write_text(json.dumps({'scope':__doc__,'rows':rows,'limits':['Uniform initial donor lift, not original point-Lz/Phi equivalent.','Residual covariance is transported separately; no physical pressure closure.','No clipping, no density floor, no energy addition.','No physical time, forces, MPI or backend parity claim.']},indent=2)+'\n')
