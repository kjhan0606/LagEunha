"""Numerical Psi diagnostic helpers; no dynamical source or fit changes."""
import numpy as np
from density_quadrature import rule
from contact_prior_density import polygon_area,clip_halfplane

def taylor_coefficients(center,potential,gravity_coefficient,rc=.5):
    center=np.asarray(center);r=np.linalg.norm(center,axis=1);u=r/rc
    g=gravity_coefficient(r);dg_over_r=np.empty(len(r))
    low=u<=.5;mid=(u>.5)&(u<1);high=u>=1
    dg_over_r[low]=(-384/5+96*u[low])/rc**5
    z=u[mid];rr=r[mid]
    dg_over_r[mid]=1/(5*rr**5)+(-48+384/5*z-32*z*z)/(rc**4*rr)
    dg_over_r[high]=-3/r[high]**5
    H=g[:,None,None]*np.eye(2)+dg_over_r[:,None,None]*center[:,:,None]*center[:,None,:]
    return dict(center=center,Phi=potential(center),gradient=g[:,None]*center,H=H)

def p2(points,coeff,k):
    ld=np.longdouble;delta=np.asarray(points,dtype=ld)-coeff['center'][k].astype(ld)
    return ld(coeff['Phi'][k])+delta@coeff['gradient'][k].astype(ld)+.5*np.einsum('ni,ij,nj->n',delta,coeff['H'][k].astype(ld),delta)

def raw_p2(raw,coeff):
    ld=np.longdouble;r=np.asarray(raw,dtype=ld);c=coeff['center'].astype(ld)
    m=r[:,0];q=r[:,1:3];Q=r[:,3:7].reshape(-1,2,2)
    centeredQ=Q-c[:,:,None]*q[:,None,:]-q[:,:,None]*c[:,None,:]+m[:,None,None]*c[:,:,None]*c[:,None,:]
    return coeff['Phi'].astype(ld)*m+np.sum(coeff['gradient'].astype(ld)*(q-m[:,None]*c),axis=1)+.5*np.sum(coeff['H'].astype(ld)*centeredQ,axis=(1,2))

def raw_volume(raw,coeff):
    ld=np.longdouble;r=np.asarray(raw,dtype=ld);c=coeff['center'].astype(ld)
    p=r[:,7:9];D=r[:,9:13].reshape(-1,2,2)
    return np.sum(coeff['gradient'].astype(ld)*p,axis=1)+np.sum(coeff['H'].astype(ld)*(D-c[:,:,None]*p[:,None,:]),axis=(1,2))

def integrate_remainder(densities,polygons,coeff,potential,orders=(24,48)):
    output=[]
    for order in orders:
        values=[]
        for k,d in enumerate(densities):
            meta=d.cut_metadata
            if meta['flat']:priors=[meta['value']]
            else:
                normal=np.asarray(meta['normal']);offset=meta['offset']
                candidates=[(clip_halfplane(polygons[k],-normal,-offset),meta['low']),
                    (clip_halfplane(polygons[k],normal,offset),meta['high'])]
                priors=[rho for poly,rho in candidates if polygon_area(poly)>0]
            assert len(priors)==len(d.parts)
            rem=np.longdouble(0);polyval=np.longdouble(0);abssum=np.longdouble(0)
            for part,rho in zip(d.parts,priors):
                pts,w=rule(part,order);dm=(w*d._value_with_prior(pts,np.full(len(w),rho))).astype(np.longdouble)
                polynomial=p2(pts,coeff,k);remainder=potential(pts.astype(np.longdouble))-polynomial
                rem+=dm@remainder;polyval+=dm@polynomial;abssum+=dm@abs(remainder)
            values.append((rem,polyval,abssum))
        output.append(np.asarray(values,dtype=np.longdouble))
    return np.asarray(output)
