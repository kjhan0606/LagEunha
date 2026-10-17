"""Positive fragment transport gate; NOT a general pressure-force Euler solver.
Affine free flight is exact Euler only for the tested constant-P nilpotent shear.
Floating convex overlap geometry is audited, never mass-renormalized.
"""
from dataclasses import dataclass, replace
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
from weighted_mass_patch import power_cells

GAMMA = 5/3

@dataclass
class Piece:
    polygon: object
    rho: float
    u: float
    b: np.ndarray
    B: np.ndarray


def rule(poly):
    """Positive triangle rule, exact through polynomial degree two."""
    if poly.geom_type != 'Polygon' or not poly.is_valid or poly.interiors:
        raise ValueError('convex simple polygon required')
    p = np.asarray(orient(poly, 1).exterior.coords[:-1]); points=[]; weights=[]
    bary = np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])
    for k in range(1,len(p)-1):
        tri=p[[0,k,k+1]]; d=tri[1]-tri[0]; e=tri[2]-tri[0]
        area=(d[0]*e[1]-d[1]*e[0])/2
        if area < -1e-16: raise RuntimeError('nonconvex polygon')
        if area <= 0: continue
        points.extend(bary@tri); weights.extend([area/3]*3)
    return np.asarray(points),np.asarray(weights)


def samples(pieces):
    xs=[];vs=[];ms=[];us=[]
    for f in pieces:
        if f.rho<=0 or f.u<=0: raise RuntimeError('nonpositive distribution')
        x,w=rule(f.polygon);xs.extend(x);vs.extend(f.b+x@f.B.T)
        ms.extend(f.rho*w);us.extend(f.u*w)
    return np.asarray(xs),np.asarray(vs),np.asarray(ms),np.asarray(us)


def raw(pieces):
    x,v,m,u=samples(pieces)
    return np.r_[m.sum(),m@x,(x.T@(m[:,None]*x)).ravel(),m@v,
                 (x.T@(m[:,None]*v)).ravel(),(v.T@(m[:,None]*v)).ravel(),u.sum()]


def invariant(r):
    return np.r_[r[0],r[7:9],r[10]-r[11],.5*(r[13]+r[16]),r[17]]


def diagnose(pieces):
    x,v,m,u=samples(pieces);mass=m.sum();c=m@x/mass;mean=m@v/mass
    dx=x-c;dv=v-mean;S=dx.T@(m[:,None]*dx);C=dx.T@(m[:,None]*dv)
    B=np.linalg.solve(S,C).T;res=dv-dx@B.T
    # The Gram matrix is the actual residual covariance, with no eigenvalue clipping.
    R=res.T@(m[:,None]*res)/mass
    subtract=(dv.T@(m[:,None]*dv)-B@S@B.T)/mass
    area=sum(f.polygon.area for f in pieces)
    return dict(mass=mass,P=(GAMMA-1)*u.sum()/area,
                min_R=float(np.linalg.eigvalsh(R).min()),
                max_R=float(np.linalg.eigvalsh(R).max()),
                subtraction_min_R=float(np.linalg.eigvalsh((subtract+subtract.T)/2).min()))


def move(f,dt):
    F=np.eye(2)+dt*f.B;J=np.linalg.det(F)
    if J<=0: raise RuntimeError('material map inverted')
    shift=dt*f.b;B=f.B@np.linalg.inv(F);b=f.b-B@shift
    p=np.asarray(f.polygon.exterior.coords[:-1])@F.T+shift
    return Piece(Polygon(p),f.rho/J,f.u/J,b,B)


def same_field(a,b):
    return a.rho==b.rho and a.u==b.u and np.array_equal(a.b,b.b) and np.array_equal(a.B,b.B)


def remap(donors,receivers):
    covered=np.zeros(len(donors));result=[];closure=0.;merge_error=0.;piece_count=0
    for cell in receivers:
        pieces=[]
        for j,f in enumerate(donors):
            p=cell.intersection(f.polygon)
            if p.area==0: continue
            if p.geom_type!='Polygon': raise RuntimeError('unexpected overlap topology')
            covered[j]+=p.area;pieces.append(replace(f,polygon=p))
        area=sum(f.polygon.area for f in pieces)
        closure=max(closure,abs(area/cell.area-1));piece_count+=len(pieces)
        if not pieces: raise RuntimeError('empty receiver')
        # Only identical physical fields can be merged. No affine projection here.
        if all(same_field(pieces[0],f) for f in pieces):
            merged=[replace(pieces[0],polygon=cell)]
            merge_error=max(merge_error,float(np.max(abs(raw(merged)-raw(pieces)))))
            result.append(merged)
        else: result.append(pieces)
    closure=max(closure,max(abs(covered[i]/f.polygon.area-1) for i,f in enumerate(donors)))
    error=float(np.max(abs(sum((raw(p) for p in result))-raw(donors))))
    if closure>2e-11 or error>2e-11 or merge_error>2e-11:
        raise RuntimeError(f'overlap conservation failed {closure=} {error=} {merge_error=}')
    return result,dict(closure=closure,remap_error=error,merge_error=merge_error,pieces=piece_count)


def edges(sites,w,cells):
    result=[]
    for i,cell in enumerate(cells):
        p=np.asarray(cell.exterior.coords[:-1]);delta=sites-sites[i]
        bound=np.sum(delta*delta,axis=1)+w[i]-w
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(b-a)
            if length<1e-14: continue
            r=abs(2*delta@((a+b)/2-sites[i])-bound);r[i]=np.inf;j=int(np.argmin(r))
            if r[j]<1e-10 and i<j: result.append((i,j,length))
    return result


def capacity(sites,w,domain,target,donors):
    # This gate uses a uniform transported density. Its integral over a receiver
    # is rho*area because the transported donor partition covers the domain.
    # Nonuniform density needs overlap-mass Newton and is deliberately rejected.
    rho=donors[0].rho
    if any(f.rho!=rho for f in donors): raise ValueError('uniform density capacity gate only')
    w=w.copy();w-=w[-1];N=len(w)
    for iteration in range(25):
        cells=power_cells(sites,w,domain);m=rho*np.array([c.area for c in cells])
        err=np.max(abs(m/target-1))
        if err<5e-13: return w,cells,iteration,err
        if min(m)<=0: raise RuntimeError('empty capacity cell')
        L=np.zeros((N,N))
        for i,j,length in edges(sites,w,cells):
            a=rho*length/(2*np.linalg.norm(sites[i]-sites[j]))
            L[i,i]+=a;L[j,j]+=a;L[i,j]-=a;L[j,i]-=a
        step=np.r_[np.linalg.solve(L[:-1,:-1],target[:-1]-m[:-1]),0.]
        for k in range(24):
            trial=w+step*2.**(-k);cc=power_cells(sites,trial,domain)
            mm=rho*np.array([c.area for c in cc])
            if min(mm)>0 and np.max(abs(mm/target-1))<err: w=trial;break
        else: raise RuntimeError('capacity line search failed')
    raise RuntimeError('capacity did not converge')


def run(n,dt,end):
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack([xx.ravel(),yy.ravel()])+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    w=np.zeros(n*n);cells=power_cells(x,w,domain)
    b0=np.array([.2,0.]);B0=np.array([[0.,1.1],[0.,0.]])
    state=[[Piece(c,1.,.6/(GAMMA-1),b0.copy(),B0.copy())] for c in cells]
    target=np.array([raw(f)[0] for f in state]);initial=sum((raw(f) for f in state));iv=invariant(initial)
    history=[];topology={(i,j) for i,j,l in edges(x,w,cells)};events=[]
    steps=round(end/dt)
    if abs(steps*dt-end)>1e-12: raise ValueError('end must be multiple of dt')
    for step in range(steps):
        flat=[f for group in state for f in group];field=flat[0]
        if not all(same_field(field,f) for f in flat): raise RuntimeError('affine patch boundary/generator gate only')
        # Point velocities come from the retained state; no centroid steering.
        x=x+dt*np.array([group[0].b+group[0].B@point for group,point in zip(state,x)])
        domain=domain@(np.eye(2)+dt*field.B).T+dt*field.b
        donors=[move(f,dt) for f in flat]
        flight_error=float(np.max(abs(invariant(raw(donors))-invariant(raw(flat)))))
        w,cells,it,cap=capacity(x,w,domain,target,donors)
        state,audit=remap(donors,cells)
        diag=[diagnose(f) for f in state];total=sum((raw(f) for f in state))
        # Analytic reference is used only below, never for state updates.
        reference=[raw([Piece(c,1.,.6/(GAMMA-1),b0,B0)]) for c in cells]
        newtop={(i,j) for i,j,l in edges(x,w,cells)}
        record=dict(t=(step+1)*dt,P_error=max(abs(d['P']-.6) for d in diag),
            mass_error=float(np.max(abs(np.array([d['mass'] for d in diag])/target-1))),
            min_R=min(d['min_R'] for d in diag),max_R=max(d['max_R'] for d in diag),
            subtraction_min_R=min(d['subtraction_min_R'] for d in diag),
            invariant_error=float(np.max(abs(invariant(total)-iv))),flight_error=flight_error,
            reference_error=max(float(np.max(abs(raw(f)-r))) for f,r in zip(state,reference)),
            capacity_error=float(cap),newton_iterations=it,**audit)
        if newtop!=topology: events.append(record['t'])
        topology=newtop;history.append(record)
    summary=dict(n=n,dt=dt,end=end,steps=steps,events=events)
    for key in ['P_error','mass_error','max_R','invariant_error','flight_error','reference_error','capacity_error','closure','remap_error','merge_error','pieces','newton_iterations']:
        summary[key]=max(r[key] for r in history)
    summary['min_R']=min(r['min_R'] for r in history)
    summary['subtraction_min_R']=min(r['subtraction_min_R'] for r in history)
    summary['pass_gate']=all(summary[k]<2e-11 for k in ['P_error','mass_error','invariant_error','reference_error','closure','remap_error','merge_error']) and summary['min_R']>=-1e-25 and summary['max_R']<1e-25
    tag=f'n{n}_dt{dt:g}_end{end:g}'
    Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
    Path(tag+'.checkpoint.json').write_text(json.dumps(dict(t=end,sites=x.tolist(),weights=w.tolist(),domain=domain.tolist(),target=target.tolist(),state=[[dict(vertices=list(f.polygon.exterior.coords),rho=f.rho,u=f.u,b=f.b.tolist(),B=f.B.tolist()) for f in group] for group in state]),indent=2))
    return summary

if __name__=='__main__':
    for n,dt,end in [(3,.005,.02),(4,.005,.02),(4,.0025,.02),(4,.00125,.02),(4,.000625,.02),(4,.0003125,.02),(4,.005,.1)]:
        result=run(n,dt,end);print(json.dumps(result),flush=True)
        if not result['pass_gate']: raise RuntimeError('transport gate failed')
