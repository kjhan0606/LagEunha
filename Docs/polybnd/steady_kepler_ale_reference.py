"""Exact-field boundary derivatives on moving Laguerre cells.

Returns raw18 plus integral rho*Phi. This is an independent steady-field
reference, not an evolution operator or a fixed-mass assumption.
"""
import numpy as np


def interface_cuts(a,b,radii):
    a=np.asarray(a);d=np.asarray(b)-a
    A=float(d@d)
    if A<=0:raise ValueError('nonzero edge required')
    B=2*float(a@d);cuts=[0.,1.]
    for radius in radii:
        C=float(a@a)-radius*radius;disc=B*B-4*A*C
        if disc<0:continue
        q=-.5*(B+np.copysign(np.sqrt(disc),B))
        roots=[-B/(2*A)] if q==0 else [q/A,C/q]
        cuts.extend(float(t) for t in roots if 0<t<1)
    return np.unique(cuts)


def field_moments(points,gravity_coefficient,potential,pressure,gamma,rin,rout,floor):
    """Pointwise Eulerian density of raw18+potential, planar circular field."""
    x=np.asarray(points);r=np.linalg.norm(x,axis=1)
    rho=np.where((r>=rin)&(r<=rout),1.,floor)
    v=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack([-x[:,1],x[:,0]])
    xx=x[:,:,None]*x[:,None,:]
    xv=x[:,:,None]*v[:,None,:];vv=v[:,:,None]*v[:,None,:]
    energy=pressure/(gamma-1)+.5*rho*np.sum(v*v,axis=1)
    return np.column_stack([rho,rho[:,None]*x,rho[:,None]*xx.reshape(-1,4),
        rho[:,None]*v,rho[:,None]*xv.reshape(-1,4),rho[:,None]*vv.reshape(-1,4),
        energy,rho*potential(x)]),v


def boundary_derivative(positions,generator_rates,weight_rates,edges,field,order=16):
    """Integrate a(x)*w_n on each physical face once, sharing both signs.

    field(points) -> (raw18+potential density, velocity).
    field.radii contains circular discontinuities requiring exact splitting.
    Exterior face normal speed follows the existing frozen-exterior rule 0.
    """
    x=np.asarray(positions);u=np.asarray(generator_rates);z=np.asarray(weight_rates)
    result=np.zeros((len(x),19));g,w=np.polynomial.legendre.leggauss(order)
    for i,j,a,b,n,length in edges:
        cuts=interface_cuts(a,b,field.radii)
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            t=lo+(hi-lo)*(g+1)/2;weights=length*(hi-lo)*w/2
            points=a+t[:,None]*(b-a)
            if j>=0:
                delta=x[j]-x[i];distance=np.linalg.norm(delta)
                wn=(delta@u[j]-(points-x[i])@(u[j]-u[i])+.5*(z[i]-z[j]))/distance
            else:wn=np.zeros(len(points))
            density,_=field(points)
            contribution=(weights*wn)@density
            result[i]+=contribution
            if j>=0:result[j]-=contribution
    return result


def exact_fixed_mass_weights(positions,generator_rates,edges,field,frozen,order=32):
    """Linear exact-field mass constraint, frozen z=0 (Dirichlet).

    Each active component must touch a frozen cell; otherwise reject the
    unanchored gauge explicitly rather than silently regularizing it.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse.linalg import spsolve
    x=np.asarray(positions);u=np.asarray(generator_rates)
    frozen=np.asarray(frozen,dtype=bool);N=len(x)
    b=np.zeros(N);ii=[];jj=[];vv=[];anchored=np.zeros(N,dtype=bool)
    g,w=np.polynomial.legendre.leggauss(order)
    for i,j,a,bb,n,length in edges:
        distance=np.linalg.norm(x[j]-x[i]) if j>=0 else None
        coeff=0.;flux=0.
        cuts=interface_cuts(a,bb,field.radii)
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            t=lo+(hi-lo)*(g+1)/2;weights=length*(hi-lo)*w/2
            points=a+t[:,None]*(bb-a);density,velocity=field(points)
            wn0=((x[j]-x[i])@u[j]-(points-x[i])@(u[j]-u[i]))/distance if j>=0 else np.zeros(len(points))
            flux+=float((weights*density[:,0])@(velocity@n-wn0))
            if j>=0:coeff+=float(weights@density[:,0])/(2*distance)
        b[i]+=flux
        if j>=0:
            b[j]-=flux
            ii.extend((i,j,i,j));jj.extend((i,j,j,i));vv.extend((coeff,coeff,-coeff,-coeff))
            if frozen[j] and not frozen[i]:anchored[i]=True
            if frozen[i] and not frozen[j]:anchored[j]=True
    matrix=coo_matrix((vv,(ii,jj)),shape=(N,N)).tocsr()
    active=np.flatnonzero(~frozen);A=matrix[active][:,active]
    count,labels=connected_components(A,directed=False)
    for k in range(count):
        if not np.any(anchored[active[labels==k]]):
            raise ValueError(f'exact mass constraint has unanchored active component {k}')
    z=np.zeros(N);z[active]=spsolve(A,b[active])
    if not np.isfinite(z).all():raise ValueError('nonfinite exact-field weight rate')
    return z,dict(active_components=int(count),linear_residual_max=float(abs((matrix@z-b)[active]).max()),
                  zero_weight_flux_max=float(abs(b[active]).max()))


def polygon_integral(poly,integrand,radii,order=16):
    """Adaptive vertical slices, explicitly split at circular interfaces.

    Gauss order controls the inner y integral; outer x is adaptive quad_vec.
    No smooth-density fit is involved.
    """
    from scipy.integrate import quad_vec
    p=np.asarray(poly);aa=p;dd=np.roll(p,-1,axis=0)-p
    knots=list(p[:,0]);lo=float(p[:,0].min());hi=float(p[:,0].max())
    for radius in radii:
        knots.extend(t for t in (-radius,radius) if lo<t<hi)
        for a,d in zip(aa,dd):
            for t in interface_cuts(a,a+d,(radius,))[1:-1]:knots.append(float(a[0]+t*d[0]))
    knots=np.unique(knots);g,w=np.polynomial.legendre.leggauss(order)
    dimension=integrand(p[:1]).shape[1]
    def vertical(xx):
        good=abs(dd[:,0])>1e-16
        t=np.zeros(len(dd));t[good]=(xx-aa[good,0])/dd[good,0]
        good&=(t>=-1e-12)&(t<=1+1e-12)
        y=aa[good,1]+t[good]*dd[good,1]
        if len(y)<2:return np.zeros(dimension)
        bottom=float(y.min());top=float(y.max());cuts=[bottom,top]
        for radius in radii:
            if abs(xx)<radius:
                height=np.sqrt(radius*radius-xx*xx)
                cuts.extend(v for v in (-height,height) if bottom<v<top)
        out=np.zeros(dimension)
        for yl,yh in zip(sorted(cuts)[:-1],sorted(cuts)[1:]):
            yy=(yl+yh)/2+(yh-yl)*g/2
            points=np.column_stack([np.full(len(yy),xx),yy])
            out+=(yh-yl)/2*(w@integrand(points))
        return out
    return quad_vec(vertical,lo,hi,points=knots[1:-1],epsabs=1e-13,epsrel=1e-11,limit=1000)


def exact_face_rates(positions,generator_rates,weight_rates,edges,field,pressure,order=32):
    """Independent ALE advection and constant-pressure face contributions.

    Pi's pressure source is zero for this spatially constant pressure field.
    D's pressure face term must still be paired with +PV I in the volume.
    """
    x=np.asarray(positions);u=np.asarray(generator_rates);z=np.asarray(weight_rates)
    adv=np.zeros((len(x),19));press=np.zeros_like(adv)
    g,w=np.polynomial.legendre.leggauss(order)
    for i,j,a,b,n,length in edges:
        cuts=interface_cuts(a,b,field.radii)
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            t=lo+(hi-lo)*(g+1)/2;weights=length*(hi-lo)*w/2
            points=a+t[:,None]*(b-a);density,v=field(points)
            wn=((x[j]-x[i])@u[j]-(points-x[i])@(u[j]-u[i])+.5*(z[i]-z[j]))/np.linalg.norm(x[j]-x[i]) if j>=0 else np.zeros(len(points))
            transport=-(weights*(v@n-wn))@density
            pterm=np.zeros(19);pterm[7:9]=-weights.sum()*pressure*n
            pterm[9:13]=-(np.outer(weights@points,pressure*n)).ravel()
            pterm[17]=-pressure*(weights@(v@n))
            adv[i]+=transport;press[i]+=pterm
            if j>=0:adv[j]-=transport;press[j]-=pterm
    return adv,press
