"""Frozen numerical midpoint field sampled on the true quadratic face path.
Audit only: does not solve integrated mass constraints or advance hydro states.
"""
import json,time
import numpy as np
from quadratic_geometry_gcl import QuadraticPath
from moving_power_spacetime import events
from kinetic_factor_rhs import ALPHA,transport,rhs
from factor_moment_evolve import advance_factors_sylvester
import moment_ale_evolve as base

class FluxAudit:
    def __init__(self,kind,n,h):
        self.path=QuadraticPath(kind,n,h,rhs_function=rhs,advance_function=advance_factors_sylvester);p=self.path;a=p.mid_state[:-4].reshape(p.N,20)
        polys,edges=base.geometry(a[:,:2],a[:,2],h/2)
        self.states=[base.decode(poly,np.r_[m,row[3:]]) for poly,m,row in zip(polys,p.mass,a)]
        centers=np.array([s['density'].origin for s in self.states]);neighbors=[set() for _ in a]
        for i,j,*_ in edges:
            if j>=0:neighbors[i].add(j);neighbors[j].add(i)
        self.grad=[]
        for i,ids in enumerate(neighbors):
            ids=list(ids);dr=centers[ids]-centers[i];length=np.linalg.norm(dr,axis=1)
            dp=np.array([self.states[j]['P']-self.states[i]['P'] for j in ids])
            g,_,rank,_=np.linalg.lstsq(dr/length[:,None],dp/length,rcond=None)
            if rank<2:raise RuntimeError('pressure stencil rank deficient')
            self.grad.append(g)
    def rate(self,t,weight_rates=None,return_jacobian=False):
        p=self.path;x,w=p.state(t);polys,edges=base.geometry(x,w,t*p.h)
        result=np.zeros((p.N,4));boundary=np.zeros(4);jac=np.zeros((p.N,p.N))
        gx,gw=np.polynomial.legendre.leggauss(4);gx=(gx+1)/2;gw=gw/2
        def field(i,point):
            s=self.states[i];v=base.velocity(s,point);rho=float(s['density'].value(point)[0]);pressure=s['P']+self.grad[i]@(point-s['density'].origin)
            if pressure<=0:raise RuntimeError('predicted pressure nonpositive')
            samples=v+ALPHA@p.mid_factors[i].T/np.sqrt(p.mass[i])
            return v,rho,pressure,samples
        for i,j,a,b,n,length in edges:
            if j>=0:row=p.row(i,j)
            else:
                mid=(a+b)/2;row=min(p.boundary,key=lambda r:abs(r[0](t)*mid[0]+r[1](t)*mid[1]-r[2](t)))
            normal=np.array([row[0](t),row[1](t)])
            for z,q in zip(gx,gw):
                point=a+z*(b-a)
                # Polynomial derivatives use normalized path time; convert to physical time.
                wn=(row[2].deriv()(t)-point@np.array([row[0].deriv()(t),row[1].deriv()(t)]))/np.linalg.norm(normal)/p.h
                if weight_rates is not None and j>=0:
                    wn+=(weight_rates[i]-weight_rates[j])/(2*np.linalg.norm(x[j]-x[i])*p.h)
                vl,rl,pl,sl=field(i,point)
                if j>=0:vr,rr,pr,sr=field(j,point)
                else:
                    vr=vl-2*(vl@n-wn)*n;rr=rl;pr=pl;sr=sl-2*(sl@n-wn)[:,None]*n
                left,right,effective_rho=transport(sl,sr,rl,rr,n,wn)
                if j>=0:
                    c=length*q*effective_rho/(2*np.linalg.norm(x[j]-x[i]))
                    jac[i,i]-=c;jac[j,j]-=c;jac[i,j]+=c;jac[j,i]+=c
                fan=base.fan(base.state(rl,vl,pl),base.state(rr,vr,pr),n);sm=fan[0][1];ps=float((fan[2][1][1:3]-fan[1][1][1:3]*sm)@n)
                fm=sum(left)+sum(right);fp=left@sl+right@sr+ps*n
                fe=left@(.5*np.sum(sl*sl,axis=1)+pl/((base.GAMMA-1)*rl))+right@(.5*np.sum(sr*sr,axis=1)+pr/((base.GAMMA-1)*rr))+ps*sm
                f=p.h*length*q*np.r_[fm,fp,fe];result[i]+=f
                if j>=0:result[j]-=f
                else:boundary+=f
        value=np.r_[result.ravel(),boundary]
        return (value,jac) if return_jacobian else value
    def integral(self,ts,order):
        gx,gw=np.polynomial.legendre.leggauss(order);total=np.zeros(4*self.path.N+4)
        for a,b in zip(ts,ts[1:]):
            for z,q in zip(gx,gw):total+=(b-a)*q/2*self.rate((a+b)/2+(b-a)*z/2)
        return total

def run(kind,n,h):
    started=time.monotonic();a=FluxAudit(kind,n,h);p=a.path;x,w=p.state(0)
    ts,_=events(x,np.zeros_like(x),w,np.zeros_like(w),constraint_builder=p.constraints)
    values=[a.integral(ts,k) for k in (8,16,32)];f=values[-1][:-4].reshape(p.N,4);boundary=values[-1][-4:]
    refinement=[float(np.max(abs(x-y))) for x,y in zip(values[:-1],values[1:])]
    r=dict(kind=kind,n=n,duration=h,status='completed',seconds=time.monotonic()-started,events=ts.tolist(),
        quadrature_differences=refinement,integrated_relative_mass_defect=float(max(abs(f[:,0])/p.mass)),boundary_mass=float(boundary[0]),ledger_error=float(max(abs(f.sum(0)-boundary))))
    r['quadrature_gate']=refinement[-1]<1e-9;r['conservation_gate']=r['ledger_error']<1e-11
    r['fixed_mass_gate']=r['integrated_relative_mass_defect']<1e-10
    return r
if __name__=='__main__':
    for kind,n,h in [('shear',2,.01),('shear',3,.01),('shear',2,.05),('wave',2,.01),('wave',3,.01),('wave',2,.05)]:
        try:r=run(kind,n,h)
        except Exception as e:r=dict(kind=kind,n=n,duration=h,status='failed',error=repr(e))
        print(json.dumps(r),flush=True)
