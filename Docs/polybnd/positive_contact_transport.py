"""Nonuniform-density contact transport gate. Uniform P, common affine shear.
Mass capacities integrate the transported piecewise density. No pressure kick.
"""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon, MultiPoint
from shapely.ops import unary_union
from positive_transport_patch import (Piece,GAMMA,raw,invariant,diagnose,move,
                                      same_field,remap,edges)
from weighted_mass_patch import power_cells,clip_halfplane


def compress_materials(pieces):
    """Merge identical fields only when their union is convex, with audits.
    This patch has two convex material supports. Reject nonconvex supports;
    a general solver must retain them, never fill their holes with a hull.
    """
    groups=[]
    for f in pieces:
        for group in groups:
            if same_field(group[0],f): group.append(f);break
        else: groups.append([f])
    result=[];area_error=0.;moment_error=0.
    for group in groups:
        union=unary_union([f.polygon for f in group])
        hull=MultiPoint([p for f in group for p in f.polygon.exterior.coords[:-1]]).convex_hull
        defect=hull.symmetric_difference(union).area
        overlap=abs(sum(f.polygon.area for f in group)-union.area)
        err=max(defect,overlap)/hull.area
        merged=replace(group[0],polygon=hull)
        merr=float(np.max(abs(raw([merged])-raw(group))))
        area_error=max(area_error,err);moment_error=max(moment_error,merr)
        if err>2e-12 or merr>2e-11: raise RuntimeError('material union not convex/conservative')
        result.append(merged)
    return result,area_error,moment_error


def evaluate(sites,w,domain,donors):
    cells=power_cells(sites,w,domain)
    mass=np.array([sum(f.rho*c.intersection(f.polygon).area for f in donors) for c in cells])
    return mass,cells


def capacity(sites,w,domain,target,donors,rtol=5e-13):
    if np.any(target<=0): raise ValueError('positive targets required')
    available=sum(f.rho*f.polygon.area for f in donors)
    if abs(available/target.sum()-1)>2e-12: raise RuntimeError('target mass incompatible')
    w=w.copy();w-=w[-1];N=len(w);eps=1e-6*max(np.ptp(domain,axis=0))**2
    for iteration in range(30):
        m,cells=evaluate(sites,w,domain,donors);error=np.max(abs(m/target-1))
        if error<rtol: return w,cells,iteration,float(error)
        if min(m)<=0: raise RuntimeError('empty capacity cell')
        # Centered differences also handle a density interface crossing a face.
        # This small reference solver favors a direct integral over speed.
        J=np.empty((N-1,N-1))
        for j in range(N-1):
            dw=np.zeros(N);dw[j]=eps
            J[:,j]=(evaluate(sites,w+dw,domain,donors)[0][:-1]-evaluate(sites,w-dw,domain,donors)[0][:-1])/(2*eps)
        step=np.r_[np.linalg.solve(J,target[:-1]-m[:-1]),0.]
        for k in range(30):
            trial=w+step*2.**(-k);mm,_=evaluate(sites,trial,domain,donors)
            if min(mm)>0 and np.max(abs(mm/target-1))<error: w=trial;break
        else: raise RuntimeError('overlap capacity line search failed')
    raise RuntimeError('overlap capacity failed')


def initial(n,contrast):
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack([xx.ravel(),yy.ravel()])+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    normal=np.array([1.,.23]);bound=.53
    b=np.array([.2,0.]);B=np.array([[0.,1.1],[0.,0.]])
    donors=[Piece(Polygon(clip_halfplane(domain,s*normal,s*bound)),rho,.6/(GAMMA-1),b.copy(),B.copy()) for s,rho in [(1,1.),(-1,float(contrast))]]
    w=np.zeros(n*n);cells=power_cells(x,w,domain);state,_=remap(donors,cells)
    return x,domain,w,cells,state,donors


def run(n,contrast,dt,end):
    x,domain,w,cells,state,reference0=initial(n,contrast)
    target=np.array([raw(f)[0] for f in state]);iv=invariant(raw(reference0))
    topology={(i,j) for i,j,l in edges(x,w,cells)};history=[];events=[]
    steps=round(end/dt)
    if abs(steps*dt-end)>1e-12: raise ValueError('end not multiple of dt')
    for step in range(steps):
        flat=[f for group in state for f in group];field=flat[0]
        if any(not np.array_equal(f.b,field.b) or not np.array_equal(f.B,field.B) for f in flat):
            raise RuntimeError('common affine velocity gate only')
        material,merge_area,merge_moment=compress_materials(flat)
        # Common physical velocity is unambiguous even at the density interface.
        x+=dt*np.array([group[0].b+group[0].B@p for group,p in zip(state,x)])
        domain=domain@(np.eye(2)+dt*field.B).T+dt*field.b
        donors=[move(f,dt) for f in material]
        flight=float(np.max(abs(invariant(raw(donors))-invariant(raw(material)))))
        w,cells,it,cap=capacity(x,w,domain,target,donors)
        state,audit=remap(donors,cells);diag=[diagnose(f) for f in state]
        total=sum((raw(f) for f in state));t=(step+1)*dt
        # Reference never enters state or capacity updates.
        exact=[move(f,t) for f in reference0];ref,_=remap(exact,cells)
        moment_error=max(float(np.max(abs(raw(a)-raw(b)))) for a,b in zip(state,ref))
        support_error=max(d.polygon.symmetric_difference(next(e.polygon for e in exact if e.rho==d.rho)).area for d in donors)
        topo={(i,j) for i,j,l in edges(x,w,cells)}
        if topo!=topology: events.append(t)
        topology=topo
        history.append(dict(t=t,P_error=max(abs(d['P']-.6) for d in diag),
            mass_error=float(np.max(abs(np.array([d['mass'] for d in diag])/target-1))),
            invariant_error=float(np.max(abs(invariant(total)-iv))),flight_error=flight,
            min_R=min(d['min_R'] for d in diag),max_R=max(d['max_R'] for d in diag),
            reference_error=moment_error,support_error=support_error,
            merge_area=merge_area,merge_moment=merge_moment,capacity_error=cap,
            newton_iterations=it,mixed_cells=sum(len(group)>1 for group in state),**audit))
    summary=dict(n=n,contrast=contrast,dt=dt,end=end,steps=steps,events=events)
    for key in history[0]:
        if key!='t': summary[key]=(min if key=='min_R' else max)(r[key] for r in history)
    checks=['P_error','mass_error','invariant_error','flight_error','reference_error','support_error','merge_area','merge_moment','closure','remap_error']
    summary['pass_gate']=all(summary[k]<2e-11 for k in checks) and summary['min_R']>=-1e-25 and summary['max_R']<1e-25
    tag=f'contact_n{n}_r{contrast:g}_dt{dt:g}_end{end:g}'
    Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
    Path(tag+'.checkpoint.json').write_text(json.dumps(dict(t=end,sites=x.tolist(),weights=w.tolist(),domain=domain.tolist(),target=target.tolist(),state=[[dict(vertices=list(f.polygon.exterior.coords),rho=f.rho,u=f.u,b=f.b.tolist(),B=f.B.tolist()) for f in group] for group in state]),indent=2))
    return summary

if __name__=='__main__':
    for n,r,dt,end in [(3,2,.005,.02),(4,2,.005,.02),(4,2,.0025,.02),(4,100,.005,.02),(4,100,.0025,.02),(4,100,.005,.1)]:
        result=run(n,r,dt,end);print(json.dumps(result),flush=True)
        if not result['pass_gate']: raise RuntimeError('contact transport gate failed')
