"""Explicit point/normal/ds nodes for true straight and circular faces.

Straight faces retain the audited quadratic-speed cuts; prior cuts use the
support representation. Wall arcs are split at actual prior intersections
and accept an independent Gauss order for a same-state integration audit.
"""
from fractions import Fraction as F
import numpy as np
from curved_disk_cell import ld,sub,dot,gauss,UnresolvedCurvedTopology
from quadratic_upwind_face_nodes import polynomial_cuts
from quadratic_velocity_reconstruction import velocity_samples,velocity_samples_and_mean_batch
from kinetic_factor_rhs import ALPHA
from arc_velocity_roots import arc_speed_cuts


def line_prior_cuts(line,owners,priors):
    delta=sub(line.B,line.A);values=[0.,1.]
    for k in owners:
        if k<0 or priors[k]['flat']:continue
        meta=priors[k];normal=tuple(F(float(v)) for v in meta['normal']);den=dot(normal,delta)
        if den==0:continue
        t=(F(float(meta['offset']))-dot(normal,line.A))/den
        if 0<t<1:
            value=float(t)
            if not 0<value<1:raise UnresolvedCurvedTopology('prior intersection parameter not representable')
            values.append(value)
    return np.unique(values)


def arc_prior_cuts(arc,meta):
    cuts=[arc.lo,arc.hi]
    if meta['flat']:return cuts
    normal=np.asarray(meta['normal']);nf=tuple(F(float(v)) for v in normal);offset=F(float(meta['offset']));radius=F(float(arc.radius))
    discriminant=radius*radius*dot(nf,nf)-offset*offset
    if discriminant<0:return cuts
    norm=np.sqrt(np.sum(normal.astype(np.longdouble)**2));ratio=ld(offset)/(np.longdouble(arc.radius)*norm)
    if discriminant>0 and abs(ratio)>=1:raise UnresolvedCurvedTopology('arc prior intersection unresolved')
    angle=float(np.arctan2(normal[1],normal[0]));opening=0. if discriminant==0 and offset>=0 else (np.pi if discriminant==0 else float(np.arccos(ratio)))
    for theta in (angle-opening,angle+opening):
        theta=theta%(2*np.pi)
        if arc.lo<theta<arc.hi:cuts.append(theta)
    return sorted(set(cuts))


def build_geometry(cells,ids,positions,mass,factors,states,gen,z,priors,line_order=8,arc_order=8,wall_only=False):
    lookup={int(g):i for i,g in enumerate(ids)};g,w=gauss(line_order);g=(g+1)/2;w=w/2
    points=[];normals=[];scales=[];owners=[];distances=[];offsets=[0];keys=[];collapsed=0;arc_audits=[]
    def append(p,n,ds,i,j,distance,key):
        p=np.asarray(p,dtype=float);n=np.asarray(n,dtype=float);ds=np.asarray(ds,dtype=float)
        if not np.isfinite(p).all() or not np.isfinite(n).all() or not np.isfinite(ds).all() or np.any(ds<=0):raise UnresolvedCurvedTopology('invalid explicit face nodes')
        points.append(p);normals.append(n);scales.append(ds);owners.append(np.tile([i,j],(len(p),1)));distances.append(np.full(len(p),distance));offsets.append(offsets[-1]+len(p));keys.append(key)
    for i,cell in enumerate(cells):
        if not wall_only:
            for line in cell.lines:
                if line.neighbor<0:raise ValueError('unclipped square wall in physical mesh')
                j=lookup[line.neighbor]
                if ids[i]>ids[j]:continue
                A=np.array([ld(v) for v in line.A]);delta=np.array([ld(v) for v in sub(line.B,line.A)]);_,normal,measure=line.quadrature(1);n=normal[0];length=measure[0]
                d=float(np.linalg.norm(positions[j]-positions[i]));dx=positions[j]-positions[i]
                if d<=0:raise ValueError('coincident fluid generators')
                if np.array_equal(np.asarray(line.A,dtype=float),np.asarray(line.B,dtype=float)):collapsed+=1
                def relative(t):
                    p=A+np.longdouble(t)*delta
                    vl=velocity_samples(states[i],factors[i],mass[i],p,ALPHA)[0][0]
                    vr=velocity_samples(states[j],factors[j],mass[j],p,ALPHA)[0][0]
                    wn=(dx@gen[j]-(gen[j]-gen[i])@(p-positions[i]))/d+(z[i]-z[j])/(2*d)
                    return np.r_[vl@n-wn,vr@n-wn]
                left,middle,right=relative(0.),relative(.5),relative(1.);c2=2*(right+left-2*middle);c1=right-left-c2
                for t in (.25,.75):
                    if not np.allclose(relative(t),left+t*c1+t*t*c2,atol=2e-12,rtol=2e-12):raise ValueError('straight normal speed is not quadratic')
                cuts=np.unique(np.r_[polynomial_cuts(left,middle,right),line_prior_cuts(line,(i,j),priors)])
                tt=np.concatenate([lo+(hi-lo)*g for lo,hi in zip(cuts[:-1],cuts[1:])]);ww=np.concatenate([(hi-lo)*w for lo,hi in zip(cuts[:-1],cuts[1:])])
                p=A+tt.astype(np.longdouble)[:,None]*delta
                append(p,np.broadcast_to(n,p.shape),length*ww,i,j,d,('line',int(ids[i]),int(ids[j])))
        for k,arc in enumerate(cell.arcs):
            root_audit=dict(id=int(ids[i]),arc=k);speed_cuts=arc_speed_cuts(arc,states[i],factors[i],mass[i],root_audit);cuts=sorted(set(arc_prior_cuts(arc,priors[i])+speed_cuts));root_audit['combined_cuts']=cuts;arc_audits.append(root_audit);gg,ww=gauss(arc_order)
            theta=np.concatenate([np.longdouble(lo)+(gg.astype(np.longdouble)+1)*np.longdouble(hi-lo)/2 for lo,hi in zip(cuts[:-1],cuts[1:])]);measure=np.concatenate([np.longdouble(arc.radius)*np.longdouble(hi-lo)*ww.astype(np.longdouble)/2 for lo,hi in zip(cuts[:-1],cuts[1:])])
            n=np.column_stack([np.cos(theta),np.sin(theta)]);p=np.longdouble(arc.radius)*n
            append(p,n,measure,i,-1,0.,('arc',int(ids[i]),k,arc.lo,arc.hi))
    arrays=dict(points=np.concatenate(points),normals=np.concatenate(normals),owners=np.concatenate(owners).astype(np.int64),record_scale=np.concatenate(scales),record_distance=np.concatenate(distances),offsets=np.asarray(offsets,dtype=np.int64),keys=keys)
    arrays['geometry_stats']=dict(nodes=len(arrays['points']),faces=len(keys),wall_nodes=int((arrays['owners'][:,1]<0).sum()),microfaces_with_collapsed_binary64_endpoints=collapsed,line_order=line_order,arc_order=arc_order,wall_only=wall_only)
    arrays['arc_root_audit']=arc_audits
    return arrays


def prepare_states(geometry,states,factors,mass):
    """Same audited per-cell batch evaluator on explicit per-node geometry."""
    points=geometry['points'];normal=geometry['normals'];owners=geometry['owners'];count=len(points);N=len(states)
    flat=owners.ravel();valid=flat>=0;receivers=flat[valid];positions=np.flatnonzero(valid);order=np.argsort(receivers,kind='stable');position_order=positions[order];sizes=np.bincount(receivers,minlength=N);offsets=np.r_[0,np.cumsum(sizes)]
    rho=np.full((count,2),np.nan);samples=np.full((count,2,4,2),np.nan);coordinates=np.full_like(samples,np.nan);means=np.full((count,2,2),np.nan);branches=np.zeros(3,dtype=np.int64)
    for i in np.flatnonzero(sizes):
        entries=position_order[offsets[i]:offsets[i+1]];rows=entries//2;sides=entries%2;p=points[rows];inward=(2*sides-1)[:,None]*normal[rows]
        values,branch=states[i]['density'].value_trace_batch(p,inward,return_branch=True)
        v,xi,mean=velocity_samples_and_mean_batch(states[i],factors[i],mass[i],p,ALPHA)
        rho[rows,sides]=values;samples[rows,sides]=v;coordinates[rows,sides]=xi;means[rows,sides]=mean;branches+=np.bincount(branch+1,minlength=3)
    cache=dict(geometry,rho=rho,samples=samples,coordinates=coordinates,means=means)
    cache['stats']=dict(nodes=count,cell_batches=int((sizes>0).sum()),point_sides=int(sizes.sum()),maximum_cell_batch=int(sizes.max(initial=0)),branch_counts_low_flat_high=branches.tolist(),cache_scope='this actual outer; explicit point/normal/ds arrays only')
    return cache
