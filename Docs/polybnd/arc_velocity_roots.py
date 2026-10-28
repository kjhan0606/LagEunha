"""Circle normal-speed zeros from the actual quadratic sample velocity.

Tangent-half-angle coordinate t in [-1,1] makes speed a degree-six
numerator over the strictly positive (1+s*s*t*t)**3. Derivative roots
partition monotone intervals; all sign crossings are bracketed there.
"""
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial import polynomial as P
from kinetic_factor_rhs import ALPHA
from quadratic_velocity_reconstruction import velocity_samples_and_mean_batch

from fractions import Fraction
LD=np.longdouble

def exact_real_roots(coeff):
    """Sturm isolation on the exact binary coefficients for uncertain signs."""
    def trim(p):
        while len(p)>1 and p[-1]==0:p.pop()
        return p
    def ev(p,x):
        v=Fraction(0)
        for a in reversed(p):v=v*x+a
        return v
    def rem(a,b):
        a=list(a)
        while len(a)>=len(b) and any(a):
            q=a[-1]/b[-1];offset=len(a)-len(b)
            for k in range(len(b)):a[offset+k]-=q*b[k]
            trim(a)
        return trim(a)
    p=trim([Fraction(*LD(x).as_integer_ratio()) for x in coeff])
    if len(p)<2:return []
    def sequence(p):
        seq=[p,[k*p[k] for k in range(1,len(p))]]
        while True:
            r=[-x for x in rem(seq[-2],seq[-1])]
            if not any(r):return seq
            scale=abs(r[-1]);seq.append([x/scale for x in r])
    seq=sequence(p)
    # Remove the exact polynomial gcd before endpoint variation tests.
    # At a repeated root every member of the original Sturm chain vanishes;
    # a square-free chain gives the correct one-sided root count there.
    gcd=seq[-1]
    if len(gcd)>1:
        remainder=list(p);quotient=[Fraction(0)]*(len(p)-len(gcd)+1)
        while len(remainder)>=len(gcd) and any(remainder):
            q=remainder[-1]/gcd[-1];offset=len(remainder)-len(gcd);quotient[offset]=q
            for k in range(len(gcd)):remainder[offset+k]-=q*gcd[k]
            trim(remainder)
        if any(remainder):raise ValueError('nonzero exact square-free division remainder')
        p=trim(quotient);seq=sequence(p)
    def variation(x):
        vals=[ev(p,x) for p in seq];signs=[1 if v>0 else -1 for v in vals if v]
        return sum(a!=b for a,b in zip(signs[:-1],signs[1:]))
    lo=Fraction(-1);hi=Fraction(1);pending=[(lo,hi,variation(lo),variation(hi))];roots=[]
    while pending:
        a,b,va,vb=pending.pop();count=va-vb
        if not count:continue
        if count==1 and b-a<Fraction(1,2**70):
            r=float((a+b)/2)
            if -1<r<1:roots.append(r)
            continue
        if b-a<Fraction(1,2**100):raise ValueError('arc polynomial roots not separable at supported precision')
        mid=(a+b)/2;vm=variation(mid);pending.extend([(a,mid,va,vm),(mid,b,vm,vb)])
    roots=sorted(roots)
    if len(roots)!=len(set(roots)):raise ValueError('distinct arc roots collapse in binary64')
    return roots


def real_crossings(coeff,audit=None):
    c=np.asarray(coeff,dtype=LD)
    while len(c)>1 and c[-1]==0:c=c[:-1]
    if not np.isfinite(c).all():raise ValueError('nonfinite arc speed polynomial')
    if len(c)<=1:return []
    critical=real_crossings(P.polyder(c),audit)
    # Near a stationary zero, floating sign tests cannot certify whether
    # two crossings or a tangent root are present. Resolve the exact stored
    # coefficient polynomial with rational Sturm arithmetic there.
    for t in critical:
        value=P.polyval(LD(t),c);scale=P.polyval(abs(LD(t)),abs(c))
        if abs(value)<=128*np.finfo(LD).eps*scale:
            if audit is not None:audit['exact_sturm_resolutions']=audit.get('exact_sturm_resolutions',0)+1
            return exact_real_roots(c)
    cuts=[-1.]+critical+[1.];roots=[]
    for a,b in zip(cuts[:-1],cuts[1:]):
        fa=P.polyval(LD(a),c);fb=P.polyval(LD(b),c)
        if fa==0 and -1<a<1:roots.append(a)
        if fb==0 and -1<b<1:roots.append(b)
        if np.signbit(fa)!=np.signbit(fb) and fa!=0 and fb!=0:
            try:roots.append(brentq(lambda t:float(P.polyval(LD(t),c)),a,b,xtol=5e-324,rtol=4*np.finfo(float).eps))
            except RuntimeError:
                if audit is not None:audit['exact_sturm_resolutions']=audit.get('exact_sturm_resolutions',0)+1
                return exact_real_roots(c)
    return sorted(set(roots))

def arc_speed_cuts(arc,state,L,mass,audit=None):
    # Short charts also cover full-circle synthetic fixtures without a pole.
    edges=np.linspace(arc.lo,arc.hi,max(1,int(np.ceil((arc.hi-arc.lo)/(np.pi/2))))+1)
    cuts=list(edges)
    if audit is not None:audit.update(charts=[],exact_sturm_resolutions=0)
    for lo,hi in zip(edges[:-1],edges[1:]):
        mid=(LD(lo)+LD(hi))/2;s=np.tan((LD(hi)-LD(lo))/4)
        co=np.cos(mid);si=np.sin(mid);D=np.array([1,0,s*s],dtype=LD);D2=P.polymul(D,D);D3=P.polymul(D2,D)
        normal=[np.array([co,-2*s*si,-s*s*co]),np.array([si,2*s*co,-s*s*si])]
        delta=[LD(arc.radius)*normal[k]-LD(state['c'][k])*D for k in range(2)]
        affine=[]
        for k in range(2):
            a=LD(state['v'][k])*D2
            for j in range(2):a=P.polyadd(a,LD(state['B'][k,j])*P.polymul(delta[j],D))
            affine.append(a)
        meta=state.get('velocity_quadratic')
        if meta is None:
            coord=[np.zeros(5,dtype=LD),np.zeros(5,dtype=LD)];residual=np.asarray(ALPHA,dtype=LD)
        else:
            Z=[(LD(arc.radius)*normal[k]-LD(meta['center'][k])*D)/LD(meta['scale']) for k in range(2)]
            phi=[P.polymul(Z[0],Z[0]),P.polymul(Z[0],Z[1]),P.polymul(Z[1],Z[1])]
            psi=[]
            for k in range(3):
                a=phi[k]
                for j,b in enumerate([D2,P.polymul(Z[0],D),P.polymul(Z[1],D)]):a=P.polysub(a,LD(meta['projection'][j,k])*b)
                psi.append(a)
            coord=[]
            for k in range(2):
                a=np.zeros(5,dtype=LD)
                for j in range(3):a=P.polyadd(a,LD(meta['standard_coefficient'][k,j])*psi[j])
                coord.append(a)
            residual=np.asarray(ALPHA,dtype=LD)@np.asarray(meta['residual_map'],dtype=LD).T
        numerators=[];chart_roots=[];chart_critical=[]
        for sample in range(len(ALPHA)):
            velocity=[]
            for k in range(2):
                a=affine[k]
                for j in range(2):
                    coordinate=P.polyadd(coord[j],residual[sample,j]*D2)
                    a=P.polyadd(a,LD(L[k,j])/np.sqrt(LD(mass))*coordinate)
                velocity.append(a)
            numerator=P.polyadd(P.polymul(velocity[0],normal[0]),P.polymul(velocity[1],normal[1]));numerators.append(numerator)
            # Stationary points are retained even without a sign crossing.
            # This also resolves integration panels around tangent roots.
            roots=real_crossings(numerator,audit);critical=real_crossings(P.polyder(numerator),audit)
            chart_roots.append(roots);chart_critical.append(critical)
            for root in sorted(set(roots+critical)):
                angle=float(mid+2*np.arctan(s*LD(root)))
                if not lo<angle<hi:raise ValueError('arc speed root not representable inside chart')
                cuts.append(angle)
        if audit is not None:audit['charts'].append(dict(lo=float(lo),hi=float(hi),speed_roots_t=chart_roots,critical_t=chart_critical))
        # Verify the algebra against the actual evaluator, with a Cartesian
        # cancellation envelope rather than a relative normal-speed bound.
        t=np.linspace(-1,1,17,dtype=LD);theta=mid+2*np.arctan(s*t);n=np.column_stack([np.cos(theta),np.sin(theta)]).astype(float);points=(LD(arc.radius)*np.column_stack([np.cos(theta),np.sin(theta)])).astype(float)
        v,_,_=velocity_samples_and_mean_batch(state,L,mass,points,ALPHA)
        actual=np.sum(v*n[:,None,:],axis=2)
        predicted=np.column_stack([P.polyval(t,c)/P.polyval(t,D3) for c in numerators])
        # Propagate the actual reconstruction primitives as well: a small
        # final Cartesian component can itself be a cancellation of v+B dx.
        cart=abs(state['v'])+(abs(points)+abs(state['c']))@abs(state['B']).T
        if meta is None:
            cart=cart[:,None,:]+(abs(ALPHA)@abs(L).T/np.sqrt(mass))[None,:,:]
        else:
            zz=abs((points-meta['center'])/meta['scale']);amp=(abs(points)+abs(meta['center']))/meta['scale']
            phiabs=np.column_stack([zz[:,0]**2,zz[:,0]*zz[:,1],zz[:,1]**2])
            basisabs=phiabs+np.column_stack([np.ones(len(points)),zz])@abs(meta['projection'])
            sensitivity=np.column_stack([2*zz[:,0]*amp[:,0],zz[:,1]*amp[:,0]+zz[:,0]*amp[:,1],2*zz[:,1]*amp[:,1]])+amp@abs(meta['projection'][1:])
            coordabs=(basisabs+sensitivity)@abs(meta['standard_coefficient']).T
            coordabs=coordabs[:,None,:]+(abs(ALPHA)@abs(meta['residual_map']).T)[None,:,:]
            cart=cart[:,None,:]+coordabs@abs(L).T/np.sqrt(mass)
        primitive=np.sum(cart*abs(n[:,None,:]),axis=2)
        bound=2048*np.finfo(float).eps*(primitive+np.sum(cart,axis=2))
        if np.any(abs(actual-predicted)>bound):raise ValueError(('arc speed polynomial disagrees with actual evaluator',float(np.max(abs(actual-predicted))),float(np.max(bound)),float(lo),float(hi)))
    if audit is not None:audit['cuts']=sorted(set(cuts))
    return sorted(set(cuts))
