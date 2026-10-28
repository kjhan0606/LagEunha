"""Positive exponential moment fit relative to a local discontinuous prior.

Research candidate. Uses current cell moments and mesh neighbors only.
No analytic interface, phase values, moment clipping or automatic fit fallback.
"""
import numpy as np
from density_quadrature import rule,validate_order


def clip_halfplane(poly,normal,offset):
    """Closed high half-plane n.x>=offset, preserving convex orientation."""
    poly=np.asarray(poly,dtype=float)
    if len(poly)==0:return np.empty((0,2))
    out=[]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        sa=float(a@normal-offset);sb=float(b@normal-offset)
        if sa>=0:out.append(a)
        if (sa<0 and sb>0) or (sa>0 and sb<0):
            out.append(a+(sa/(sa-sb))*(b-a))
    return np.asarray(out).reshape(-1,2)


def polygon_area(poly):
    if len(poly)<3:return 0.
    # Center before taking cross products to avoid cancellation far from 0.
    p=np.asarray(poly)-np.asarray(poly)[0]
    return float(.5*np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)))


def features(points,origin,scale):
    z=(np.atleast_2d(points)-origin)/scale
    return np.column_stack([z,z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2])


def local_contact_priors(polys,raw,edges,vertex_errors=None,mass_errors=None):
    from geometry_uncertainty import mean_density_interval,subtract,divide_positive
    raw=np.asarray(raw);N=len(polys)
    areas=np.array([polygon_area(p) for p in polys])
    if np.any(areas<=0) or np.any(raw[:,0]<=0):raise ValueError('positive areas and mass required')
    centers=[]
    for poly in polys:
        points,w=rule(np.asarray(poly),6);centers.append(w@points/w.sum())
    centers=np.asarray(centers);means=raw[:,0]/areas
    if vertex_errors is not None and len(vertex_errors)!=N:raise ValueError('vertex error count mismatch')
    if mass_errors is None:mass_errors=np.zeros(N)
    bounds=[mean_density_interval(poly,areas[i],raw[i,0],None if vertex_errors is None else vertex_errors[i],mass_errors[i]) for i,poly in enumerate(polys)]
    rho_bounds=np.asarray([b[1] for b in bounds])
    neighbors=[set() for _ in range(N)]
    for i,j,*_ in edges:
        if j>=0:neighbors[i].add(j);neighbors[j].add(i)
    result=[];eps=64*np.finfo(float).eps
    for i,poly in enumerate(polys):
        stencil={i}|neighbors[i]
        for j in neighbors[i]:stencil.update(neighbors[j])
        indices=np.array(sorted(stencil));low=float(means[indices].min());high=float(means[indices].max())
        low_interval=(float(rho_bounds[indices,0].min()),float(rho_bounds[indices,1].min()))
        high_interval=(float(rho_bounds[indices,0].max()),float(rho_bounds[indices,1].max()))
        metadata=dict(low=low,high=high,stencil_cells=len(indices),area_bounds=list(bounds[i][0]),
            rho_bounds=rho_bounds[i].tolist(),low_bounds=list(low_interval),high_bounds=list(high_interval),
            fraction_bounds=None,
            coordinate_error_model='supplied' if vertex_errors is not None else 'working_affine_gamma4',
            maximum_vertex_error=float(np.max(bounds[i][2])),mass_error=float(mass_errors[i]))
        def constant(reason):return dict(metadata,flat=True,value=float(means[i]),reason=reason)
        if float(rho_bounds[indices,0].max())<=float(rho_bounds[indices,1].min()):
            result.append(constant('unresolved_local_contrast'));continue
        contrast_bounds=subtract(high_interval,low_interval)
        if contrast_bounds[0]<=0:
            result.append(constant('unresolved_local_contrast'));continue
        fraction_bounds=divide_positive(subtract(tuple(rho_bounds[i]),low_interval),contrast_bounds)
        metadata['fraction_bounds']=[float(v) for v in fraction_bounds]
        fraction=float((means[i]-low)/(high-low))
        if fraction_bounds[0]<=0<=fraction_bounds[1] or fraction_bounds[0]<=1<=fraction_bounds[1]:
            result.append(constant('unresolved_phase_fraction'));continue
        others=indices[indices!=i];scale=np.sqrt(areas[i]);delta=(centers[others]-centers[i])/scale
        distance=np.linalg.norm(delta,axis=1)
        if np.any(distance==0):raise ValueError('coincident stencil centroids')
        grad,_,rank,_=np.linalg.lstsq(delta/distance[:,None],((means[others]-means[i])/high)/distance,rcond=None)
        if rank==2 and np.linalg.norm(grad)>eps:
            normal=grad/np.linalg.norm(grad);direction='neighbor_gradient'
        else:
            shift=(raw[i,1:3]/raw[i,0]-centers[i])/scale
            if np.linalg.norm(shift)<=eps:
                result.append(constant('unidentified_direction'));continue
            normal=shift/np.linalg.norm(shift);direction='mass_centroid'
        lo=float(np.min(np.asarray(poly)@normal));hi=float(np.max(np.asarray(poly)@normal))
        for _ in range(64):
            mid=(lo+hi)/2
            area=polygon_area(clip_halfplane(poly,normal,mid))
            if area>fraction*areas[i]:lo=mid
            else:hi=mid
        result.append(dict(metadata,flat=False,normal=normal.tolist(),offset=(lo+hi)/2,
                           fraction=fraction,direction=direction))
    return result


class ContactPriorDensity:
    def __init__(self,poly,m,q,Q,prior_metadata,order=12):
        self.quadrature_order=validate_order(order);poly=np.asarray(poly,dtype=float)
        if poly.ndim!=2 or poly.shape[1]!=2 or len(poly)<3 or not np.isfinite(poly).all():
            raise ValueError('finite convex polygon required')
        if polygon_area(poly)<=0 or not np.isfinite(m) or m<=0:
            raise ValueError('positive oriented area and mass required')
        self.mass=float(m);self.cut_metadata=dict(prior_metadata)
        if prior_metadata['flat']:
            value=float(prior_metadata['value']);parts=[(poly,value)]
        else:
            n=np.asarray(prior_metadata['normal'],dtype=float);offset=float(prior_metadata['offset'])
            if n.shape!=(2,) or not np.isfinite(n).all() or np.linalg.norm(n)==0 or not np.isfinite(offset):
                raise ValueError('finite nonzero cut normal required')
            parts=[(clip_halfplane(poly,-n,-offset),float(prior_metadata['low'])),
                   (clip_halfplane(poly,n,offset),float(prior_metadata['high']))]
        points=[];weights=[];priors=[];self.parts=[]
        for part,rho in parts:
            if not np.isfinite(rho) or rho<=0:raise ValueError('positive finite prior required')
            area=polygon_area(part)
            if area<0:raise ValueError('negative clipped area')
            if area==0:continue
            pts,w=rule(part,self.quadrature_order)
            if np.any(w<=0):raise ValueError('nonpositive quadrature weight')
            points.append(pts);weights.append(w);priors.append(np.full(len(w),rho));self.parts.append(part)
        self.points=np.concatenate(points);self.weights=np.concatenate(weights)
        self.area=float(self.weights.sum());self.origin=self.weights@self.points/self.area;self.scale=np.sqrt(self.area)
        if abs(self.area/polygon_area(poly)-1)>1e-11:raise ValueError('cut partition lost geometric area')
        prior=np.concatenate(priors);F=features(self.points,self.origin,self.scale)
        ld=np.longdouble;o=self.origin.astype(ld);c=np.asarray(q,dtype=ld)/ld(m)
        second=np.asarray(Q,dtype=ld)/ld(m)-np.outer(o,c)-np.outer(c,o)+np.outer(o,o)
        target=np.asarray(np.r_[(c-o)/ld(self.scale),second[0,0]/ld(self.area),second[0,1]/ld(self.area),second[1,1]/ld(self.area)],dtype=float)
        if not np.isfinite(target).all():raise ValueError('nonfinite target moments')
        theta=np.zeros(5)
        def evaluate(t):
            exponent=F@t+np.log(prior);shift=float(exponent.max())
            weighted=self.weights*np.exp(exponent-shift);Z=float(weighted.sum());prob=weighted/Z
            mean=prob@F;center=F-mean;H=center.T@(prob[:,None]*center)
            return mean-target,H,shift+np.log(Z),prob
        for iteration in range(60):
            residual,H,logZ,prob=evaluate(theta);norm=float(abs(residual).max())
            if norm<3e-13:break
            step=np.linalg.solve(H,residual)
            for damping in range(30):
                trial=theta-step*2.**(-damping)
                if abs(evaluate(trial)[0]).max()<norm:theta=trial;break
            else:raise RuntimeError('contact prior moment fit stalled')
        else:raise RuntimeError('contact prior moment fit did not converge')
        self.theta=theta;self.logZ=logZ;self.residual=norm;self.iterations=iteration+1
        self.condition=float(np.linalg.cond(H));self.weighted_mass=self.mass*prob

    def value(self,points):
        points=np.atleast_2d(points);meta=self.cut_metadata
        prior=np.full(len(points),meta['value']) if meta['flat'] else np.where(points@np.asarray(meta['normal'])>=meta['offset'],meta['high'],meta['low'])
        return self._value_with_prior(points,prior)

    def value_trace(self,points,inward_direction):
        """One-sided cell-interior trace at points exactly on the cut."""
        points=np.atleast_2d(points);meta=self.cut_metadata
        if meta['flat']:return self.value(points)
        normal=np.asarray(meta['normal']);signed=points@normal-meta['offset']
        direction=np.broadcast_to(np.asarray(inward_direction,dtype=float),points.shape)
        if not np.isfinite(direction).all():raise ValueError('finite inward direction required')
        derivative=direction@normal;oncut=signed==0
        if np.any(oncut&(derivative==0)):
            raise ValueError('tangent direction does not define a one-sided cut trace')
        high=(signed>0)|(oncut&(derivative>0))
        prior=np.where(high,meta['high'],meta['low'])
        return self._value_with_prior(points,prior)

    def value_trace_batch(self,points,inward_direction,return_branch=False):
        """Cell batch with original one-point dot kernels and exact cut rules.

        A leading batch axis retains the scalar (1,2)@normal kernel rather
        than changing the cut predicate to a large matrix/vector product.
        Branch codes are0 for flat,1 for high,-1 for low; no tolerance.
        """
        points=np.asarray(points,dtype=float)
        if points.ndim!=2 or points.shape[1]!=2:
            raise ValueError('trace batch requires (nodes,2) points')
        meta=self.cut_metadata
        if meta['flat']:
            prior=np.full(len(points),meta['value'])
            branch=np.zeros(len(points),dtype=np.int8)
        else:
            normal=np.asarray(meta['normal'])
            signed=(points[:,None,:]@normal)[:,0]-meta['offset']
            direction=np.broadcast_to(np.asarray(inward_direction,dtype=float),points.shape)
            if not np.isfinite(direction).all():raise ValueError('finite inward direction required')
            derivative=(direction[:,None,:]@normal)[:,0]
            oncut=signed==0
            if np.any(oncut&(derivative==0)):
                raise ValueError('tangent direction does not define a one-sided cut trace')
            high=(signed>0)|(oncut&(derivative>0))
            prior=np.where(high,meta['high'],meta['low'])
            branch=np.where(high,1,-1).astype(np.int8)
        F=features(points,self.origin,self.scale)
        exponent=(F[:,None,:]@self.theta)[:,0]+np.log(prior)-self.logZ
        value=self.mass*np.exp(exponent)
        if not np.isfinite(value).all() or np.any(value<=0):raise RuntimeError('invalid contact density value')
        return (value,branch) if return_branch else value

    def _value_with_prior(self,points,prior):
        value=self.mass*np.exp(features(points,self.origin,self.scale)@self.theta+np.log(prior)-self.logZ)
        if not np.isfinite(value).all() or np.any(value<=0):raise RuntimeError('invalid contact density value')
        return value
