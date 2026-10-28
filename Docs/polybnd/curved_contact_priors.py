"""Local density priors on actual curved domains, with numerical intervals.

Uses numerical cell masses/areas and the graph only. No analytic phase values,
interface radii, density smoothing or geometry/weight steering enters a prior.
"""
from fractions import Fraction
import numpy as np
from geometry_uncertainty import mean_density_interval,outward,subtract,divide_positive


def geometric_density_bounds(cell,area,mass):
    if not cell.arcs:
        polygon=np.array([[float(x),float(y)] for x,y in cell.vertices])
        ab,rb,error=mean_density_interval(polygon,area,mass)
        return dict(area_bounds=list(ab),rho_bounds=list(rb),maximum_vertex_error=float(np.max(error)),coordinate_error_model='working_affine_gamma4',mass_error=0.)
    points,normals,ds=cell.boundary_quadrature()
    u=np.finfo(float).eps/2;gamma4=4*u/(1-4*u)
    point_error=gamma4*(abs(points)+np.sqrt(area))
    normal_error=gamma4*(abs(normals)+1)
    primitive=ds*np.sum(points*normals,axis=1)
    operations=128+8*len(points);gamma=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    error=(gamma*abs(primitive).sum()+np.sum(ds*np.sum(point_error*abs(normals)+abs(points)*normal_error+point_error*normal_error,axis=1)))/2
    error=float(error)+np.finfo(float).eps*abs(area)
    ab=outward(area-error,area+error)
    if ab[0]<=0:raise ValueError('curved cell area unresolved in working envelope')
    rb=divide_positive(outward(mass,mass),ab)
    return dict(area_bounds=list(map(float,ab)),rho_bounds=list(map(float,rb)),maximum_vertex_error=float(point_error.max()),coordinate_error_model='curved_boundary_working_roundoff_not_upstream_certificate',mass_error=0.)


def support_extrema(cell,normal):
    values=[];normal=np.asarray(normal,dtype=float)
    for line in cell.lines:
        values.extend(float(sum(Fraction(float(n))*x for n,x in zip(normal,p))) for p in (line.A,line.B))
    angle=float(np.arctan2(normal[1],normal[0])%(2*np.pi));norm=float(np.linalg.norm(normal))
    for arc in cell.arcs:
        for theta in (arc.lo,arc.hi):values.append(float(arc.radius*(normal@np.array([np.cos(theta),np.sin(theta)]))))
        for theta,value in ((angle,arc.radius*norm),((angle+np.pi)%(2*np.pi),-arc.radius*norm)):
            if arc.lo<=theta<=arc.hi:values.append(value)
    if not values:raise ValueError('empty support')
    return np.nextafter(min(values),-np.inf),np.nextafter(max(values),np.inf)


def region_area(cell):
    if cell is None:return 0.
    if not cell.arcs:
        p=cell.vertices
        twice=sum((a[0]*b[1]-a[1]*b[0] for a,b in zip(p,p[1:]+p[:1])),Fraction(0))
        if twice<=0:raise ValueError('nonpositive exact clipped area')
        return float(twice/2)
    return cell.area


def prior_for_cell(i,cell,raw,areas,centers,bounds,neighbors,means):
    stencil={i}|set(neighbors[i])
    for j in neighbors[i]:stencil.update(neighbors[j])
    indices=np.array(sorted(stencil));low=float(means[indices].min());high=float(means[indices].max())
    rb=np.array([bounds[j]['rho_bounds'] for j in indices])
    low_interval=(float(rb[:,0].min()),float(rb[:,1].min()));high_interval=(float(rb[:,0].max()),float(rb[:,1].max()))
    metadata=dict(bounds[i],low=low,high=high,stencil_cells=len(indices),low_bounds=list(low_interval),high_bounds=list(high_interval),fraction_bounds=None)
    def constant(reason):return dict(metadata,flat=True,value=float(means[i]),reason=reason)
    if float(rb[:,0].max())<=float(rb[:,1].min()):return constant('unresolved_local_contrast')
    contrast_bounds=subtract(high_interval,low_interval)
    if contrast_bounds[0]<=0:return constant('unresolved_local_contrast')
    fraction_bounds=divide_positive(subtract(tuple(bounds[i]['rho_bounds']),low_interval),contrast_bounds)
    metadata['fraction_bounds']=[float(v) for v in fraction_bounds]
    fraction=float((means[i]-low)/(high-low))
    if fraction_bounds[0]<=0<=fraction_bounds[1] or fraction_bounds[0]<=1<=fraction_bounds[1]:return constant('unresolved_phase_fraction')
    others=indices[indices!=i];scale=np.sqrt(areas[i]);delta=(centers[others]-centers[i])/scale;distance=np.linalg.norm(delta,axis=1)
    if np.any(distance==0):raise ValueError('coincident stencil centroids')
    grad,_,rank,_=np.linalg.lstsq(delta/distance[:,None],((means[others]-means[i])/high)/distance,rcond=None);eps=64*np.finfo(float).eps
    if rank==2 and np.linalg.norm(grad)>eps:normal=grad/np.linalg.norm(grad);direction='neighbor_gradient'
    else:
        shift=(raw[i,1:3]/raw[i,0]-centers[i])/scale
        if np.linalg.norm(shift)<=eps:return constant('unidentified_direction')
        normal=shift/np.linalg.norm(shift);direction='mass_centroid'
    lo,hi=support_extrema(cell,normal)
    for _ in range(64):
        mid=(lo+hi)/2
        area=region_area(cell.clip_halfplane(-normal,-mid))
        if area>fraction*areas[i]:lo=mid
        else:hi=mid
    return dict(metadata,flat=False,normal=normal.tolist(),offset=(lo+hi)/2,fraction=fraction,direction=direction)
