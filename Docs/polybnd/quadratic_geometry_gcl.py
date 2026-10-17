"""Geometry-only audit along a quadratic path from two numerical RHS stages.
No hydro state is corrected. No reference flow is supplied to interior fluxes.
"""
import json,time
import numpy as np
from numpy.polynomial import Polynomial as Poly
from adaptive_moment_batch import setup,D0
from compatible_moment_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors
from moving_power_spacetime import events
from power_spacetime_geometry import basis,moments
from weighted_mass_patch import power_cells

class QuadraticPath:
    def __init__(self,kind,n,h,rhs_function=None,advance_function=None):
        stage_rhs=rhs if rhs_function is None else rhs_function
        stage_advance=advance_factors if advance_function is None else advance_function
        mass,y=setup(kind,n,.01);L=[np.zeros((2,2)) for _ in mass]
        dy,d=stage_rhs(0.,y,mass,linear_pressure=True,covariance_factors=L,return_covariance_terms=True)
        hmL=stage_advance(L,d,h/2);hm=rebuild(y+h*dy/2,mass,hmL)
        dh,_=stage_rhs(h/2,hm,mass,linear_pressure=True,covariance_factors=hmL,return_covariance_terms=True)
        a=y[:-4].reshape(-1,20);f=dy[:-4].reshape(-1,20);g=dh[:-4].reshape(-1,20)
        self.mass=mass.copy();self.mid_state=hm.copy();self.mid_factors=hmL
        self.N=len(mass);self.kind=kind;self.h=h
        self.x=[[Poly([a[i,k],h*f[i,k],h*(g[i,k]-f[i,k])]) for k in range(2)] for i in range(self.N)]
        self.w=[Poly([a[i,2],h*f[i,2],h*(g[i,2]-f[i,2])]) for i in range(self.N)]
        self.boundary=[[Poly([-1]),Poly([0,1.1*h]),Poly([0,-.2*h])],
                       [Poly([1]),Poly([0,-1.1*h]),Poly([1,.2*h])],
                       [Poly([0]),Poly([-1]),Poly([0])],[Poly([0]),Poly([1]),Poly([1])]] if kind=='shear' else [[Poly([a]),Poly([b]),Poly([c])] for a,b,c in [(-1,0,0),(1,0,1),(0,-1,0),(0,1,1)]]
    def domain(self,t):
        p=D0.copy()
        if self.kind=='shear':p[:,0]+=self.h*t*(.2+1.1*p[:,1])
        return p
    def state(self,t):return np.array([[p(t) for p in x] for x in self.x]),np.array([p(t) for p in self.w])
    def row(self,i,j):
        return [2*(self.x[j][k]-self.x[i][k]) for k in range(2)]+[sum(p*p for p in self.x[j])-sum(p*p for p in self.x[i])+self.w[i]-self.w[j]]
    def constraints(self,x,v,w,dw,i):return [self.row(i,j) for j in range(self.N) if j!=i]+self.boundary
    def rate(self,t):
        x,w=self.state(t);cells=power_cells(x,w,self.domain(t));rate=np.zeros((self.N,6));boundary=np.zeros(6)
        q,qw=np.polynomial.legendre.leggauss(3);q=(q+1)/2;qw=qw/2
        for i,c in enumerate(cells):
            if c.area<=0:raise RuntimeError('empty cell on path')
            p=np.array(c.exterior.coords[:-1])
            candidates=[(j,self.row(i,j)) for j in range(self.N) if j!=i]+[(-1,row) for row in self.boundary]
            for a,b in zip(p,np.roll(p,-1,axis=0)):
                length=np.linalg.norm(b-a)
                if length<1e-13:continue
                mid=(a+b)/2
                j,row=min(candidates,key=lambda item:abs(sum(item[1][k](t)*mid[k] for k in range(2))-item[1][2](t)))
                residual=abs(sum(row[k](t)*mid[k] for k in range(2))-row[2](t))
                if residual>1e-9:raise RuntimeError('unclassified face')
                if j>=0 and i>=j:continue
                points=a+q[:,None]*(b-a);normal=np.array([row[0](t),row[1](t)])
                speed=(row[2].deriv()(t)-points@np.array([row[0].deriv()(t),row[1].deriv()(t)]))/np.linalg.norm(normal)
                exchange=length*((qw*speed)@basis(points));rate[i]+=exchange
                if j>=0:rate[j]-=exchange
                else:boundary+=exchange
        return rate,boundary
    def integral(self,ts,order):
        q,qw=np.polynomial.legendre.leggauss(order);total=np.zeros((self.N,6));boundary=np.zeros(6)
        for a,b in zip(ts,ts[1:]):
            for z,w in zip(q,qw):
                r,s=self.rate((a+b)/2+(b-a)*z/2);total+=(b-a)*w/2*r;boundary+=(b-a)*w/2*s
        return total,boundary

def run(kind,n,h):
    started=time.monotonic();path=QuadraticPath(kind,n,h);x,w=path.state(0.)
    ts,_=events(x,np.zeros_like(x),w,np.zeros_like(w),constraint_builder=path.constraints)
    low,_=path.integral(ts,8);high,boundary=path.integral(ts,16)
    x1,w1=path.state(1.);actual=moments(power_cells(x1,w1,path.domain(1.)))-moments(power_cells(x,w,path.domain(0.)))
    error=np.max(abs(high-actual),axis=0);mid,_=path.rate(.5)
    result=dict(kind=kind,n=n,duration=h,status='completed',events=ts.tolist(),seconds=time.monotonic()-started,
                area_error=float(error[0]),first_moment_error=float(max(error[1:3])),second_moment_error=float(max(error[3:])),
                refinement_error=float(np.max(abs(high-low))),boundary_ledger_error=float(np.max(abs(high.sum(0)-boundary))),
                midpoint_geometry_error=float(np.max(abs(mid-actual))))
    result['pass_gate']=max(result[k] for k in ['area_error','first_moment_error','second_moment_error','refinement_error','boundary_ledger_error'])<1e-10
    return result
if __name__=='__main__':
    for kind in ('shear','wave'):
        for n in (2,3):
            for h in (.01,.05):
                try:r=run(kind,n,h)
                except Exception as e:r=dict(kind=kind,n=n,duration=h,status='failed',error=repr(e),pass_gate=False)
                print(json.dumps(r),flush=True)
