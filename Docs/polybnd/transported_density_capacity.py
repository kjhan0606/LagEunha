"""Capacity solve in a retained, positive polygonal numerical density field.

Independent affine-contact reference: it retains subcell information beyond the
six cell moments. No exact physical interface is accepted by this module.
"""
import numpy as np
import shapely
from shapely.geometry import Polygon,LineString
from shapely.strtree import STRtree
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from active_power_cells import power_cells_active

def convex_line_intervals(vertices,start,end,spatial_enclosure=0.):
    """Clip a line against convex pieces using local long-double arithmetic."""
    v=np.asarray(vertices,dtype=np.longdouble);start=np.asarray(start,dtype=np.longdouble)
    d=np.asarray(end,dtype=np.longdouble)-start
    edge=np.roll(v,-1,axis=1)-v;local=v-v[:,:1,:]
    signed=np.sum(local[:,:,0]*np.roll(local[:,:,1],-1,axis=1)-local[:,:,1]*np.roll(local[:,:,0],-1,axis=1),axis=1)
    sign=np.where(signed>=0,1.,-1.)[:,None]
    offset=start-v;value=(edge[:,:,0]*offset[:,:,1]-edge[:,:,1]*offset[:,:,0])*sign
    if spatial_enclosure:
        # Only a generalized Newton matrix uses this enclosure. Exact overlap
        # integrals, field geometry, masses, and acceptance remain unchanged.
        value+=np.longdouble(spatial_enclosure)*np.sqrt(np.sum(edge*edge,axis=2))
    slope=(edge[:,:,0]*d[1]-edge[:,:,1]*d[0])*sign
    roots=np.divide(-value,slope,out=np.zeros_like(value),where=slope!=0)
    lower=np.maximum(0.,np.max(np.where(slope>0,roots,-np.inf),axis=1))
    upper=np.minimum(1.,np.min(np.where(slope<0,roots,np.inf),axis=1))
    valid=~np.any((slope==0)&(value<0),axis=1)
    upper=np.where(valid,upper,-np.inf)
    return np.asarray(lower,dtype=float),np.asarray(upper,dtype=float)


def interval_average(lower,upper,rho):
    lower=np.maximum(0.,lower);upper=np.minimum(1.,upper);use=upper>lower
    lower=lower[use];upper=upper[use];rho=rho[use]
    points,inverse=np.unique(np.r_[0.,1.,lower,upper],return_inverse=True)
    sums=np.zeros(len(points));counts=np.zeros(len(points));n=len(lower)
    np.add.at(sums,inverse[2:2+n],rho);np.add.at(sums,inverse[2+n:],-rho)
    np.add.at(counts,inverse[2:2+n],1.);np.add.at(counts,inverse[2+n:],-1.)
    sums=np.cumsum(sums)[:-1];counts=np.cumsum(counts)[:-1];width=np.diff(points)
    missing=float(width[counts==0].sum())
    value=np.divide(sums,counts,out=np.zeros_like(sums),where=counts>0)
    return float(value@width),missing


def polygon_raw_many(polygons,origin):
    """Seven area moments, computed locally before translation to global x."""
    count=len(polygons);out=np.zeros((count,7))
    coordinates,indices=shapely.get_coordinates(polygons,return_index=True)
    if len(coordinates)==0:return out
    coordinates=coordinates-np.asarray(origin)
    good=indices[:-1]==indices[1:];a=coordinates[:-1][good];b=coordinates[1:][good];ids=indices[:-1][good]
    cross=a[:,0]*b[:,1]-b[:,0]*a[:,1]
    terms=np.column_stack([cross/2,(a[:,0]+b[:,0])*cross/6,(a[:,1]+b[:,1])*cross/6,
        (a[:,0]**2+a[:,0]*b[:,0]+b[:,0]**2)*cross/12,
        (2*a[:,0]*a[:,1]+a[:,0]*b[:,1]+b[:,0]*a[:,1]+2*b[:,0]*b[:,1])*cross/24,
        (a[:,1]**2+a[:,1]*b[:,1]+b[:,1]**2)*cross/12])
    local=np.zeros((count,6))
    for j in range(6):local[:,j]=np.bincount(ids,weights=terms[:,j],minlength=count)
    local*=np.where(local[:,0]>=0,1.,-1.)[:,None]
    m=local[:,0];q=local[:,1:3];Q=np.zeros((count,2,2))
    Q[:,0,0]=local[:,3];Q[:,0,1]=Q[:,1,0]=local[:,4];Q[:,1,1]=local[:,5]
    c=np.asarray(origin)
    globalQ=Q+q[:,:,None]*c[None,None,:]+c[None,:,None]*q[:,None,:]+m[:,None,None]*np.outer(c,c)
    out[:,0]=m;out[:,1:3]=q+m[:,None]*c;out[:,3:7]=globalQ.reshape(-1,4)
    out[shapely.area(polygons)<=0]=0.
    return out


class PolygonDensity:
    def __init__(self,vertices,densities):
        try:self.vertices=np.asarray(vertices,dtype=float)
        except ValueError:self.vertices=[np.asarray(v,dtype=float) for v in vertices]
        self.densities=np.asarray(densities,dtype=float)
        if len(self.vertices)!=len(self.densities) or not len(self.vertices):raise ValueError('invalid polygon density size')
        if not np.all(np.isfinite(self.densities)) or np.any(self.densities<=0):raise ValueError('strictly positive density required')
        self.polygons=shapely.polygons(self.vertices) if isinstance(self.vertices,np.ndarray) else np.asarray([Polygon(v) for v in self.vertices],dtype=object)
        if not np.all(shapely.is_valid(self.polygons)) or np.any(shapely.area(self.polygons)<=0):
            raise ValueError('positive valid polygons required')
        if np.any(shapely.get_num_interior_rings(self.polygons)) or not np.all(shapely.equals(self.polygons,shapely.convex_hull(self.polygons))):
            raise ValueError('convex polygons without holes required')
        self.tree=STRtree(self.polygons)
        self.total_mass=float(self.densities@shapely.area(self.polygons))
        self.jacobian_enclosures=[]

    def affine_image(self,F,shift):
        F=np.asarray(F,dtype=float);shift=np.asarray(shift,dtype=float);J=float(np.linalg.det(F))
        if F.shape!=(2,2) or shift.shape!=(2,) or not np.isfinite(J) or J<=0:raise ValueError('invalid affine map')
        vertices=self.vertices@F.T+shift if isinstance(self.vertices,np.ndarray) else [v@F.T+shift for v in self.vertices]
        return PolygonDensity(vertices,self.densities/J)

    def integrate(self,cells,moments=False,batch_size=32):
        cells=np.asarray(cells,dtype=object);mass=np.zeros(len(cells));coverage=np.zeros(len(cells))
        raw=np.zeros((len(cells),7)) if moments else None
        for start in range(0,len(cells),batch_size):
            group=cells[start:start+batch_size];pairs=self.tree.query(group,predicate='intersects')
            ids,source=pairs
            if not len(ids):continue
            overlaps=shapely.intersection(group[ids],self.polygons[source]);areas=shapely.area(overlaps)
            mass[start:start+len(group)]=np.bincount(ids,weights=areas*self.densities[source],minlength=len(group))
            coverage[start:start+len(group)]=np.bincount(ids,weights=areas,minlength=len(group))
            if moments:
                # Each target's local origin avoids cancellation in tiny overlaps.
                for j,cell in enumerate(group):
                    use=ids==j
                    if not np.any(use):continue
                    center=np.asarray(cell.centroid.coords)[0]
                    raw[start+j]=self.densities[source[use]]@polygon_raw_many(overlaps[use],center)
        return mass,coverage,raw

    def line_integral(self,start,end):
        """Integrate rho ds, averaging coincident one-sided polygon traces.

        Intervals are partitioned before summation, so shared triangle edges
        cannot double the Jacobian. Off source interfaces this is the unique
        density trace; at jumps the average is a generalized Newton choice.
        """
        start=np.asarray(start);end=np.asarray(end);delta=end-start;length=float(np.linalg.norm(delta))
        if length==0:return 0.,0.
        ids=self.tree.query(LineString([start,end]))
        def intervals(enclosure=0.):
            if isinstance(self.vertices,np.ndarray):return convex_line_intervals(self.vertices[ids],start,end,enclosure)
            pairs=[convex_line_intervals(self.vertices[i][None,:,:],start,end,enclosure) for i in ids]
            return np.asarray([v[0][0] for v in pairs]),np.asarray([v[1][0] for v in pairs])
        lower,upper=intervals();value,missing=interval_average(lower,upper,self.densities[ids])
        if missing>2e-10:
            # Independently rounded shared vertices can leave a microscopic
            # strip. Near a tangent face its parameter length is amplified.
            # A fixed coordinate-roundoff enclosure supplies a quasi-Newton
            # trace only; it never fills the density or changes overlap mass.
            enclosure=32*np.finfo(float).eps*max(1.,float(np.max(abs(np.r_[start,end]))))
            box=shapely.box(min(start[0],end[0])-enclosure,min(start[1],end[1])-enclosure,
                max(start[0],end[0])+enclosure,max(start[1],end[1])+enclosure)
            ids=self.tree.query(box)
            lower,upper=intervals(enclosure);newvalue,newmissing=interval_average(lower,upper,self.densities[ids])
            self.jacobian_enclosures.append(dict(spatial_enclosure=enclosure,original_uncovered_fraction=missing,
                enclosed_uncovered_fraction=newmissing,length=length,trace_integral_change=length*(newvalue-value)))
            value,missing=newvalue,newmissing
        if missing>2e-10:raise RuntimeError('uncovered capacity face '+str(dict(fraction=missing,start=start.tolist(),end=end.tolist(),length=length)))
        return float(length*value),missing


def capacity_edges(sites,weights,cells):
    sites=np.asarray(sites);weights=np.asarray(weights)
    for i,cell in enumerate(cells):
        p=np.asarray(cell.exterior.coords[:-1]);delta=sites-sites[i]
        bound=np.sum(delta*delta,axis=1)+weights[i]-weights
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            if np.linalg.norm(b-a)<1e-14:continue
            residual=abs(2*delta@((a+b)/2-sites[i])-bound);residual[i]=np.inf
            j=int(np.argmin(residual))
            if residual[j]<1e-10 and i<j:yield i,j,a,b


def mass_jacobian(sites,weights,cells,field):
    rows=[];columns=[];values=[];missing=0.
    for i,j,a,b in capacity_edges(sites,weights,cells):
        integral,gap=field.line_integral(a,b);missing=max(missing,gap)
        c=integral/(2*np.linalg.norm(sites[j]-sites[i]))
        rows.extend((i,j,i,j));columns.extend((i,j,j,i));values.extend((c,c,-c,-c))
    return coo_matrix((values,(rows,columns)),shape=(len(sites),len(sites))).tocsc(),missing


def solve_capacity(sites,weights,domain,target,field,rtol=5e-11,max_iterations=40,callback=None):
    sites=np.asarray(sites,dtype=float);weights=np.asarray(weights,dtype=float).copy();target=np.asarray(target,dtype=float)
    if np.any(target<=0) or target.shape!=(len(sites),):raise ValueError('positive capacity targets required')
    compatibility=abs(field.total_mass/target.sum()-1.)
    if compatibility>1e-11:raise ValueError('source density and fixed total mass incompatible: '+str(compatibility))
    weights-=weights[-1];history=[]
    for iteration in range(max_iterations):
        cells=power_cells_active(sites,weights,domain)
        if any(p.is_empty for p in cells):raise RuntimeError('empty initial capacity cell')
        mass,coverage,_=field.integrate(cells);error=float(np.max(abs(mass-target)/target))
        cover=float(np.max(abs(coverage/np.asarray([p.area for p in cells])-1.)))
        row=dict(iteration=iteration,mass_relative=error,coverage_relative=cover)
        history.append(row)
        if callback is not None:callback(weights,cells,history)
        if error<rtol:
            if cover>1e-9:raise RuntimeError('transported density partition failed')
            return weights,cells,history
        if np.any(mass<=0):raise RuntimeError('empty initial capacity mass')
        J,gap=mass_jacobian(sites,weights,cells,field);row['uncovered_face_fraction']=gap
        row['jacobian_roundoff_enclosures']=len(field.jacobian_enclosures)
        step=np.r_[spsolve(J[:-1,:-1],target[:-1]-mass[:-1]),0.]
        if not np.all(np.isfinite(step)):raise RuntimeError('singular capacity graph')
        for backtrack in range(30):
            trial=weights+step*2.**(-backtrack);newcells=power_cells_active(sites,trial,domain)
            if any(p.is_empty for p in newcells):continue
            newmass,_,_=field.integrate(newcells)
            newerror=float(np.max(abs(newmass-target)/target))
            if np.min(newmass)>0 and newerror<error:
                weights=trial;row['backtracks']=backtrack;break
        else:raise RuntimeError('transported density capacity line search stalled')
    raise RuntimeError('transported density capacity iterations exhausted')
