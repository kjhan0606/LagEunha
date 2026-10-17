"""Manufactured affine Euler flux gate on a moving material box.
Weights solve integrated ALE mass, not endpoint areas. Endpoint integrals audit it.
Uniform rho/P exact shear supplies flux states: not a general hydro closure.
"""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scipy.optimize import root
from weighted_mass_patch import power_cells
from moving_power_spacetime import constraints,events
from positive_transport_patch import rule

GAMMA=5/3;P0=.6
D0=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
def velocity(x):return np.column_stack([.2+1.1*x[:,1],np.zeros(len(x))])


def boundary_rows(T):
    dv=T*velocity(D0);rows=[]
    for i in range(4):
        j=(i+1)%4
        a=[Poly([D0[i,k],dv[i,k]]) for k in range(2)]
        delta=[Poly([D0[j,k]-D0[i,k],dv[j,k]-dv[i,k]]) for k in range(2)]
        n=[delta[1],-delta[0]];rows.append(n+[n[0]*a[0]+n[1]*a[1]])
    return rows


def endpoint(x,w,T):
    cells=power_cells(x,w,D0+T*velocity(D0));states=[]
    for cell in cells:
        if cell.area<=0:raise ValueError('empty cell')
        r,q=rule(cell);v=velocity(r);E=P0/(GAMMA-1)+.5*np.sum(v*v,axis=1)
        L=r[:,0]*v[:,1]-r[:,1]*v[:,0]
        states.append(q@np.column_stack([np.ones(len(q)),v,E,L]))
    return np.array(states)


def integrate(x0,w0,dw,T,order=8,split=True):
    dx=T*velocity(x0);brows=boundary_rows(T)
    def builder(x,v,w,d,i):return constraints(x,v,w,d,i)[:-4]+brows
    ts,_=events(x0,dx,w0,dw,builder) if split else (np.array([0.,1.]),0)
    qs,ws=np.polynomial.legendre.leggauss(3);qs=(qs+1)/2;ws=ws/2
    qt,wt=np.polynomial.legendre.leggauss(order);total=np.zeros((len(x0),5));boundary=np.zeros(5)
    edge_mass={};tops=[]
    def rate(t):
        x=x0+t*dx;w=w0+t*dw;domain=D0+t*T*velocity(D0)
        cells=power_cells(x,w,domain);flux=np.zeros_like(total);boundflux=np.zeros(5);byedge={};top=[]
        for i,cell in enumerate(cells):
            if cell.area<=0:raise ValueError('empty cell')
            p=np.asarray(cell.exterior.coords[:-1]);delta=x-x[i];rhs=np.sum(delta*delta,axis=1)+w[i]-w
            for a,b in zip(p,np.roll(p,-1,axis=0)):
                length=np.linalg.norm(b-a)
                if length<1e-13:continue
                # power_cells preserves CCW orientation from the input domain.
                n=np.array([b[1]-a[1],a[0]-b[0]])/length
                err=abs(2*delta@((a+b)/2-x[i])-rhs);err[i]=np.inf;j=int(np.argmin(err))
                if err[j]>1e-9:j=-1
                if j>=0 and i>=j:continue
                r=a+qs[:,None]*(b-a);v=velocity(r);vn=v@n
                if j>=0:
                    speed=(x[j]@dx[j]-x[i]@dx[i]-r@(dx[j]-dx[i])+.5*(dw[i]-dw[j]))/np.linalg.norm(x[j]-x[i])
                    top.append((i,j))
                else:
                    # Geometric boundary speed from its implicit moving line.
                    mid=(a+b)/2
                    k=min(range(4),key=lambda k:abs(brows[k][0](t)*mid[0]+brows[k][1](t)*mid[1]-brows[k][2](t)))
                    row=brows[k];normal=np.array([row[0](t),row[1](t)])
                    speed=(row[2].deriv()(t)-r@np.array([row[0].deriv()(t),row[1].deriv()(t)]))/np.linalg.norm(normal)
                fm=T*vn-speed
                fp=v*fm[:,None]+T*P0*n
                E=P0/(GAMMA-1)+.5*np.sum(v*v,axis=1)
                fe=E*fm+T*P0*vn
                fl=r[:,0]*fp[:,1]-r[:,1]*fp[:,0]
                f=length*(ws@np.column_stack([fm,fp,fe,fl]));flux[i]+=f
                if j>=0:flux[j]-=f;byedge[(i,j)]=float(f[0])
                else:boundflux+=f
        return flux,boundflux,byedge,sorted(top)
    for left,right in zip(ts,ts[1:]):
        mid=(left+right)/2;half=(right-left)/2;tops.append(rate(mid)[3])
        for z,a in zip(qt,wt):
            f,b,edges,_=rate(mid+half*z);q=half*a;total+=q*f;boundary+=q*b
            for edge,value in edges.items():edge_mass[edge]=edge_mass.get(edge,0.)+q*value
    return total,boundary,dict(events=ts.tolist(),topologies=tops,max_face_mass=float(max(map(abs,edge_mass.values()),default=0.)))


def run(T,seed,order=8):
    x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
    if seed:x+=np.random.default_rng(seed).uniform(-.025,.025,(4,2))
    w=np.array([.01,0.,0.,0.]);initial=endpoint(x,w,0.)
    calls=0
    def residual(z):
        nonlocal calls;calls+=1
        return integrate(x,w,np.r_[z,0.],T,order)[0][:-1,0]
    solution=root(residual,np.zeros(3),tol=1e-9)
    dw=np.r_[solution.x,0.];flux,boundary,audit=integrate(x,w,dw,T,order)
    final=endpoint(x+T*velocity(x),w+dw,T)
    defect=final-initial+flux;mass=float(np.max(abs(flux[:,0])))
    identity=float(np.max(abs(defect)));relative=float(np.max(abs(final[:,0]/initial[:,0]-1)))
    finer=integrate(x,w,dw,T,16)[0];unsplit=integrate(x,w,dw,T,16,False)[0]
    refinement=float(np.max(abs(finer-flux)));unsplit_error=float(np.max(abs(final-initial+unsplit)))
    # Reconstruct positivity from independently integrated endpoint raw moments.
    cells=power_cells(x+T*velocity(x),w+dw,D0+T*velocity(D0))
    geometric_mass=np.array([c.area for c in cells]);bulkP=(GAMMA-1)*(final[:,3]-.5*np.sum(final[:,1:3]**2,axis=1)/final[:,0])/geometric_mass
    result=dict(duration=T,seed=seed,order=order,newton_success=bool(solution.success),newton_message=str(solution.message),evaluations=calls,
        mass_flux_error=mass,endpoint_relative_mass_error=relative,full_ale_identity_error=identity,
        global_boundary_ledger_error=float(np.max(abs(flux.sum(0)-boundary))),
        quadrature_refinement_error=refinement,unsplit_identity_error=unsplit_error,
        min_endpoint_mass=float(min(final[:,0])),min_pressure_reference=P0,
        naive_mean_velocity_pressure_excess=float(max(bulkP-P0)),**audit)
    result['pass_gate']=bool(mass<2e-12 and relative<2e-11 and identity<2e-12 and refinement<2e-12)
    tag=f'T{T:g}_seed{seed}'
    Path(tag+'.checkpoint.json').write_text(json.dumps(dict(sites_initial=x.tolist(),weights_initial=w.tolist(),weight_increment=dw.tolist(),initial=initial.tolist(),final=final.tolist(),integrated_flux=flux.tolist(),boundary=boundary.tolist(),diagnostics=result),indent=2))
    return result

if __name__=='__main__':
    for T,seed in [(.02,0),(.1,0),(.3,0),(.1,731),(.3,731),(.3,732)]:
        r=run(T,seed);print(json.dumps(r),flush=True)
        if not r['pass_gate']:raise RuntimeError('ALE mass gate failed')
