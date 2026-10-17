"""Implicit mass-constrained weight path on a manufactured material shear box.
Solve capacities at every queried time; differentiate the constraint for w'.
Exact affine fluid fields still supply fluxes: not general hydro evolution.
"""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.integrate import quad_vec
from positive_transport_patch import capacity,edges
from spacetime_ale_mass import D0,velocity,endpoint,boundary_rows,P0,GAMMA


def flux_rate(x,dx,w,z,T,t,cells,brows):
    flux=np.zeros((len(x),5));boundary=np.zeros(5);top=[];max_face=0.
    q,a=np.polynomial.legendre.leggauss(3);q=(q+1)/2;a=a/2
    for i,cell in enumerate(cells):
        p=np.asarray(cell.exterior.coords[:-1]);delta=x-x[i];bound=np.sum(delta*delta,axis=1)+w[i]-w
        for left,right in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(right-left)
            if length<1e-13:continue
            n=np.array([right[1]-left[1],left[0]-right[0]])/length
            err=abs(2*delta@((left+right)/2-x[i])-bound);err[i]=np.inf;j=int(np.argmin(err))
            if err[j]>1e-9:j=-1
            if j>=0 and i>=j:continue
            r=left+q[:,None]*(right-left);v=velocity(r);vn=v@n
            if j>=0:
                speed=(x[j]@dx[j]-x[i]@dx[i]-r@(dx[j]-dx[i])+.5*(z[i]-z[j]))/np.linalg.norm(x[j]-x[i]);top.append((i,j))
            else:
                mid=(left+right)/2;k=min(range(4),key=lambda k:abs(brows[k][0](t)*mid[0]+brows[k][1](t)*mid[1]-brows[k][2](t)))
                row=brows[k];nn=np.array([row[0](t),row[1](t)])
                speed=(row[2].deriv()(t)-r@np.array([row[0].deriv()(t),row[1].deriv()(t)]))/np.linalg.norm(nn)
            fm=T*vn-speed;fp=v*fm[:,None]+T*P0*n
            energy=P0/(GAMMA-1)+.5*np.sum(v*v,axis=1);fe=energy*fm+T*P0*vn
            fl=r[:,0]*fp[:,1]-r[:,1]*fp[:,0]
            f=length*(a@np.column_stack([fm,fp,fe,fl]));flux[i]+=f
            if j>=0:flux[j]-=f;max_face=max(max_face,abs(float(f[0])))
            else:boundary+=f
    return flux,boundary,tuple(sorted(top)),max_face


class PathSolution:
    def __init__(self,T,seed):
        self.T=T;self.x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
        if seed:self.x+=np.random.default_rng(seed).uniform(-.025,.025,(4,2))
        self.dx=T*velocity(self.x);self.w0=np.array([.01,0.,0.,0.]);self.initial=endpoint(self.x,self.w0,0.)
        self.target=self.initial[:,0];self.brows=boundary_rows(T);self.cache={}
    def at(self,t):
        t=float(t)
        if t in self.cache:return self.cache[t]
        x=self.x+t*self.dx;domain=D0+t*self.T*velocity(D0)
        w,cells,it,error=capacity(x,self.w0,domain,self.target,[SimpleNamespace(rho=1.)])
        zero=np.zeros(4);base=flux_rate(x,self.dx,w,zero,self.T,t,cells,self.brows)[0][:,0]
        L=np.zeros((4,4))
        for i,j,length in edges(x,w,cells):
            c=length/(2*np.linalg.norm(x[j]-x[i]));L[i,i]+=c;L[j,j]+=c;L[i,j]-=c;L[j,i]-=c
        z=np.r_[np.linalg.solve(L[:-1,:-1],base[:-1]),0.]
        f,b,top,maxface=flux_rate(x,self.dx,w,z,self.T,t,cells,self.brows)
        record=dict(w=w,z=z,flux=f,boundary=b,topology=top,mass_error=float(error),
                    instantaneous_mass_flux=float(np.max(abs(f[:,0]))),max_face=maxface)
        self.cache[t]=record;return record
    def breaks(self,count):
        grid=np.linspace(0,1,count);points=[]
        for lo,hi in zip(grid,grid[1:]):
            old=self.at(lo)['topology']
            if old==self.at(hi)['topology']:continue
            for _ in range(32):
                mid=(lo+hi)/2
                if self.at(mid)['topology']==old:lo=mid
                else:hi=mid
            points.append((lo+hi)/2)
        return points
    def integral(self,points,tol):
        def fn(t):
            a=self.at(t);return np.r_[a['flux'].ravel(),a['boundary']]
        value,error,info=quad_vec(fn,0.,1.,epsabs=tol,epsrel=tol,points=points,full_output=True,limit=300)
        if not info.success:raise RuntimeError('adaptive time quadrature failed: '+info.message)
        return value[:20].reshape(4,5),value[20:],float(error)


def run(T,seed):
    path=PathSolution(T,seed);b17=path.breaks(17);b33=path.breaks(33)
    if len(b17)!=len(b33) or (b17 and max(abs(np.array(b17)-b33))>1e-8):raise RuntimeError('event sampling audit failed')
    coarse,_,_=path.integral(b33,1e-9);flux,boundary,est=path.integral(b33,1e-12)
    last=path.at(1);final=endpoint(path.x+path.dx,last['w'],T)
    derivative_error=0.
    for t in (.125,.375,.625,.875):
        h=1e-5
        fd=(path.at(t+h)['w']-path.at(t-h)['w'])/(2*h)
        derivative_error=max(derivative_error,float(np.max(abs(fd-path.at(t)['z']))))
    records=list(path.cache.values())
    result=dict(duration=T,seed=seed,events=b33,evaluated_times=len(records),
        max_sampled_mass_error=max(r['mass_error'] for r in records),
        max_instantaneous_mass_flux=max(r['instantaneous_mass_flux'] for r in records),
        max_face_mass_rate=max(r['max_face'] for r in records),
        integrated_mass_error=float(np.max(abs(flux[:,0]))),
        ale_identity_error=float(np.max(abs(final-path.initial+flux))),
        boundary_ledger_error=float(np.max(abs(flux.sum(0)-boundary))),
        quadrature_refinement_error=float(np.max(abs(flux-coarse))),quadrature_error_estimate=est,
        weight_derivative_fd_error=derivative_error)
    result['pass_gate']=bool(result['max_sampled_mass_error']<1e-11 and result['max_instantaneous_mass_flux']<1e-12 and result['ale_identity_error']<2e-11 and result['weight_derivative_fd_error']<2e-8)
    tag=f'T{T:g}_seed{seed}'
    Path(tag+'.path.json').write_text(json.dumps([dict(t=t,w=r['w'].tolist(),z=r['z'].tolist(),mass_error=r['mass_error'],mass_flux=r['instantaneous_mass_flux']) for t,r in sorted(path.cache.items())],indent=2))
    Path(tag+'.checkpoint.json').write_text(json.dumps(dict(initial=path.initial.tolist(),final=final.tolist(),flux=flux.tolist(),boundary=boundary.tolist(),diagnostics=result),indent=2))
    return result

if __name__=='__main__':
    for T,seed in [(.02,0),(.1,0),(.3,0),(.1,731),(.3,731),(.3,732)]:
        r=run(T,seed);print(json.dumps(r),flush=True)
        if not r['pass_gate']:raise RuntimeError('constrained path failed')
