"""Batch the audited straight-face cuts/fill; true arcs stay in the oracle."""
import ctypes as C
import hashlib
import os
from pathlib import Path

import numpy as np
from curved_disk_cell import ld, sub, gauss, UnresolvedCurvedTopology
from quadratic_velocity_reconstruction import velocity_samples
from kinetic_factor_rhs import ALPHA

CAPACITY=24
_library=None


def library():
    global _library
    if _library is None:
        _library=C.CDLL(str(Path(os.environ['LAGFORCE_NATIVE_STRAIGHT_FACE_LIB']).resolve()))
        _library.lag_line_cuts.argtypes=[C.c_int64]+[C.c_void_p]*10+[C.c_int64]+[C.c_void_p]*3
        _library.lag_line_cuts.restype=C.c_int
        _library.lag_line_fill.argtypes=[C.c_int64]*3+[C.c_void_p]*12
        _library.lag_line_fill.restype=C.c_int
    return _library


def ptr(value):return value.ctypes.data_as(C.c_void_p)


def prepare(cells,ids,positions,mass,factors,states,gen,priors,stage_cache,prior_cut):
    if stage_cache is not None and hasattr(stage_cache,'native_straight_constants'):
        result=stage_cache.native_straight_constants
        stage_cache.line_hits+=len(result['keys'])
        return result
    lookup={int(g):i for i,g in enumerate(ids)}
    rows=[];keys=[];prior=[];prior_offsets=[0];collapsed=[]
    for i,cell in enumerate(cells):
        for line in cell.lines:
            if line.neighbor<0:raise ValueError('unclipped square wall in physical mesh')
            j=lookup[line.neighbor]
            if ids[i]>ids[j]:continue
            cached=None if stage_cache is None else stage_cache.lines.get((i,j))
            if cached is None:
                A=np.array([ld(v) for v in line.A]);delta=np.array([ld(v) for v in sub(line.B,line.A)])
                _,normal,measure=line.quadrature(1);n=normal[0];length=measure[0]
                dx=positions[j]-positions[i];distance=float(np.linalg.norm(dx))
                if distance<=0:raise ValueError('coincident fluid generators')
                samples={}
                for t in (0.,.5,1.,.25,.75):
                    point=A+np.longdouble(t)*delta
                    vl=velocity_samples(states[i],factors[i],mass[i],point,ALPHA)[0][0]
                    vr=velocity_samples(states[j],factors[j],mass[j],point,ALPHA)[0][0]
                    samples[t]=(vl@n,vr@n)
                cached=(A,delta,n,length,distance,dx,samples,prior_cut(line,(i,j),priors))
                if stage_cache is not None:stage_cache.lines[i,j]=cached
            elif stage_cache is not None:stage_cache.line_hits+=1
            A,delta,n,length,distance,dx,samples,cuts=cached
            rows.append((A,delta,n,length,distance,dx@gen[j],gen[j]-gen[i],positions[i],
                         np.array([np.r_[samples[t][0],samples[t][1]] for t in (0.,.5,1.,.25,.75)])))
            keys.append((i,j));prior.extend(cuts);prior_offsets.append(len(prior))
            collapsed.append(np.array_equal(np.asarray(line.A,dtype=float),np.asarray(line.B,dtype=float)))
    count=len(rows)
    if not count:raise ValueError('Native straight-face batch is empty')
    result=dict(keys=keys,index={key:k for k,key in enumerate(keys)},collapsed=collapsed,
        A=np.ascontiguousarray([r[0] for r in rows],dtype=np.longdouble),
        delta=np.ascontiguousarray([r[1] for r in rows],dtype=np.longdouble),
        normal=np.ascontiguousarray([r[2] for r in rows],dtype=np.longdouble),
        length=np.ascontiguousarray([r[3] for r in rows],dtype=np.longdouble),
        distance=np.ascontiguousarray([r[4] for r in rows],dtype=float),
        base_dot=np.ascontiguousarray([r[5] for r in rows],dtype=float),
        dg=np.ascontiguousarray([r[6] for r in rows],dtype=float),
        position=np.ascontiguousarray([r[7] for r in rows],dtype=float),
        samples=np.ascontiguousarray([r[8] for r in rows],dtype=np.longdouble),
        prior_offsets=np.asarray(prior_offsets,dtype=np.int64),prior=np.asarray(prior,dtype=float))
    assert result['samples'].shape==(count,5,8)
    if stage_cache is not None:stage_cache.native_straight_constants=result
    return result


def evaluate(constants,z,order):
    n=len(constants['keys']);owners=np.asarray(constants['keys'],dtype=np.int64)
    dz=np.ascontiguousarray(z[owners[:,0]]-z[owners[:,1]],dtype=float)
    cuts=np.empty((n,CAPACITY));sizes=np.empty(n,dtype=np.int64);status=np.empty(n,dtype=np.int64)
    keys=('A','delta','samples','base_dot','dg','position','distance')
    inputs=[constants[k] for k in keys]+[dz,constants['prior_offsets'],constants['prior']]
    code=library().lag_line_cuts(n,*(ptr(v) for v in inputs),CAPACITY,ptr(cuts),ptr(sizes),ptr(status))
    if code or np.any(status):
        bad=np.flatnonzero(status)
        raise UnresolvedCurvedTopology(dict(status='native_straight_speed_cut_failed',code=code,
            indices=bad[:64].tolist(),face_status=status[bad[:64]].tolist()))
    offsets=np.r_[0,np.cumsum((sizes-1)*order)].astype(np.int64)
    count=int(offsets[-1]);points=np.empty((count,2));normals=np.empty((count,2));measure=np.empty(count)
    g,w=gauss(order);g=np.ascontiguousarray((g+1)/2);w=np.ascontiguousarray(w/2)
    inputs=[constants[k] for k in ('A','delta','normal','length')]+[cuts,sizes,offsets,g,w,points,normals,measure]
    code=library().lag_line_fill(n,order,CAPACITY,*(ptr(v) for v in inputs))
    if code or not np.isfinite(points).all() or not np.isfinite(normals).all() or not np.isfinite(measure).all() or np.any(measure<=0):
        raise UnresolvedCurvedTopology('Invalid native straight face quadrature')
    return dict(constants=constants,points=points,normals=normals,measure=measure,offsets=offsets,cuts=cuts,sizes=sizes)


def provenance():
    path=Path(os.environ['LAGFORCE_NATIVE_STRAIGHT_FACE_LIB']).resolve()
    return dict(library=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                arc_algorithm_changed=False,quadrature_order_changed=False,
                ordered_panels=True,coordinate_arithmetic='long double through final binary64 export')
