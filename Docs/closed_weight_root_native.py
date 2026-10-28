# Native primitives with original closed-root diagnostics/Newton/sparse solve.
# Original source SHA256: a7a0a5c737d21009137583a190f6423a6ccd18df4e7234d83e1c8e826b64ff89
"""Closed fixed-mass kinetic root with an explicit single weight gauge.

L=-d(residual)/dz is a positive weighted graph Laplacian. Solve its reduced
system, then return to mean(z)=0. No frozen fluid, ridge, lsmr fallback,
residual correction or boundary leakage is introduced.
"""
import warnings
import native_transport_root
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import MatrixRankWarning,spsolve
from vector_kinetic_root import record_arrays

def evaluate(data,count,z,with_jacobian=True):
    i=data['i'];j=data['j'];inside=j>=0
    prepared=data.get('_native_prepared') or native_transport_root.Prepared(data)
    native=native_transport_root.evaluate(prepared,count,z)
    wl,wr,wn=(native[k] for k in ('wl','wr','wn'))
    fm,primitive,residual=(native[k] for k in ('fm','primitive','residual'))
    eps=np.finfo(float).eps;operations=128+4*len(i)+2*count
    gamma=operations*eps/(1-operations*eps)
    wall=~inside
    # A cold tangential sample has a tiny computed normal speed while its
    # Cartesian dot product subtracts large terms. Bound the reflection and
    # dot arithmetic using those actual absolute Cartesian terms.
    projection=np.asarray(data.get('wall_projection_absolute_sum',primitive))
    normal_error=np.asarray(data.get('normal_norm_defect',np.zeros(len(i))))
    wall_bound=128*eps*projection[wall]+2*normal_error[wall]*primitive[wall]
    compatibility_bound=gamma*primitive.sum()+wall_bound.sum()
    compatibility=float(residual.sum())
    if abs(compatibility)>compatibility_bound:raise RuntimeError(dict(status='incompatible_closed_mass_residual',residual=compatibility,bound=float(compatibility_bound)))
    if np.any(abs(fm[wall])>wall_bound):raise RuntimeError(dict(status='nonzero_closed_wall_massflux',maximum=float(abs(fm[wall]).max()),bound_max=float(wall_bound.max(initial=0.))))
    diag=dict(compatibility_residual=compatibility,compatibility_bound=float(compatibility_bound),wall_mass_flux_max=float(abs(fm[wall]).max(initial=0.)),wall_mass_flux_bound_max=float(wall_bound.max(initial=0.)),primitive_mass_flux_absolute_sum=float(primitive.sum()))
    L=None
    if with_jacobian:
        rho=native['rho']
        weight=native['weight'][inside]
        L=coo_matrix((np.r_[weight,weight,-weight,-weight],(prepared.jacobian_rows,prepared.jacobian_columns)),shape=(count,count)).tocsc();L.eliminate_zeros()
        ones=L@np.ones(count);absolute=np.asarray(abs(L).sum(axis=1)).ravel();degree=np.asarray(L.getnnz(axis=1)).ravel()
        null_bound=(128+2*degree)*eps*absolute
        if np.any(abs(ones)>null_bound):raise RuntimeError('closed kinetic Jacobian lost constant-weight nullspace')
        components,labels=connected_components(L,directed=False)
        if components!=1:raise RuntimeError(dict(status='disconnected_closed_kinetic_jacobian',components=int(components),component_sizes=np.bincount(labels).tolist(),fallback=False))
        diag.update(nullspace_residual_max=float(abs(ones).max(initial=0.)),nullspace_bound_max=float(null_bound.max(initial=0.)),components=int(components),jacobian_nnz=int(L.nnz))
    return residual,L,(wl,wr,wn),diag

def solve_arrays(data,count,initial=None,tolerance=2e-13):
    if count<1:raise ValueError('nonempty fluid registry required')
    i=data['i'];j=data['j'];inside=j>=0
    if np.any((i<0)|(i>=count)|(j<-1)|(j>=count)) or np.any(i[inside]==j[inside]):raise ValueError('invalid closed node owners')
    if (~inside).any() and 'wall_projection_absolute_sum' not in data:raise ValueError('Cartesian projection envelope required for closed wall nodes')
    if np.any(data['distance'][inside]<=0) or np.any(data['base'][~inside]!=0):raise ValueError('positive internal distances and stationary wall required')
    if np.any(data['scale']<=0) or np.any(data['rl']<=0) or np.any(data['rr']<=0):raise ValueError('positive node measure/density required')
    if any(not np.isfinite(data[k]).all() for k in ('scale','base','distance','vl','vr','rl','rr')):raise ValueError('finite closed kinetic inputs required')
    if np.any(data.get('epsilon',0.)!=0):raise ValueError('unsmoothed kinetic closure required')
    data=dict(data);data['_native_prepared']=native_transport_root.Prepared(data)
    z=np.zeros(count) if initial is None else np.asarray(initial,dtype=float).copy()
    if z.shape!=(count,) or not np.isfinite(z).all():raise ValueError('invalid initial gauge rates')
    z-=z.mean();history=[]
    for iteration in range(50):
        residual,L,values,diag=evaluate(data,count,z);norm=float(abs(residual).max());history.append(dict(iteration=iteration,norm=norm,**diag))
        if norm<tolerance:
            return z,values,norm,dict(iterations=iteration,history=history,gauge='mean weight-rate zero; all fluid active',gauge_sum=float(z.sum()),fallback=False)
        if count==1:raise RuntimeError('one-cell closed residual exceeds roundoff')
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            try:step=spsolve(L[:-1,:-1],residual[:-1])
            except MatrixRankWarning as exc:raise RuntimeError('unresolved closed reduced gauge solve') from exc
        if not np.isfinite(step).all():raise RuntimeError('nonfinite closed Newton step')
        dz=np.r_[step,0.];dz-=dz.mean()
        linear_residual=L@dz-residual
        linear_bound=128*np.finfo(float).eps*(abs(L)@abs(dz)+abs(residual))+diag['compatibility_bound']
        if np.any(abs(linear_residual)>linear_bound):raise RuntimeError(dict(status='closed_linear_solve_residual',maximum=float(abs(linear_residual).max()),bound=float(linear_bound.max())))
        for damping in range(30):
            trial=z+dz*2.**(-damping);trial-=trial.mean()
            r,_,_,_=evaluate(data,count,trial,False)
            if abs(r).max()<norm:z=trial;break
        else:raise RuntimeError('closed massroot line search failed')
    raise RuntimeError('closed massroot did not converge')

def solve(records,count,initial=None,tolerance=2e-13,return_diagnostics=False,array_values=False):
    data=record_arrays(records)
    normals=np.asarray([r[9] for r in records]);vl=np.asarray([r[5] for r in records]);vr=np.asarray([r[6] for r in records])
    defect=abs(np.sum(normals*normals,axis=1)-1)
    if np.any(defect>64*np.finfo(float).eps):raise ValueError('unit face normals required')
    projectionL=(abs(vl)@abs(normals)[:,:,None])[:,:,0];projectionR=(abs(vr)@abs(normals)[:,:,None])[:,:,0]
    data['wall_projection_absolute_sum']=data['scale']*.25*(data['rl']*projectionL.sum(1)+data['rr']*projectionR.sum(1))
    data['normal_norm_defect']=defect
    z,(wl,wr,wn),norm,diag=solve_arrays(data,count,initial,tolerance)
    result=(z,(wl,wr,wn) if array_values else [(a,b,float(c)) for a,b,c in zip(wl,wr,wn)],norm)
    return result+(diag,) if return_diagnostics else result
