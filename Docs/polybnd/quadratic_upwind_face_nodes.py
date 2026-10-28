"""Positive Gauss quadrature split at all quadratic normal-speed zeros."""
import numpy as np

def polynomial_cuts(left,middle,right):
    left=np.asarray(left);middle=np.asarray(middle);right=np.asarray(right)
    c2=2*(right+left-2*middle);c1=right-left-c2;cuts=[0.,1.]
    for a,b,c in zip(c2.ravel(),c1.ravel(),left.ravel()):
        scale=max(abs(a),abs(b),abs(c))
        if scale==0:continue
        a,b,c=a/scale,b/scale,c/scale
        if abs(a)<=32*np.finfo(float).eps:
            roots=[] if b==0 else [-c/b]
        else:
            disc=b*b-4*a*c
            if disc<0:continue
            q=-.5*(b+np.copysign(np.sqrt(disc),b))
            roots=[-b/(2*a)] if q==0 else [q/a,c/q]
        cuts.extend(float(r) for r in roots if 0<r<1)
    return np.unique(cuts)

class QuadraticUpwindNodes:
    velocity_degree=2
    def __init__(self,x,mass,factors,z,order,boundary_velocity,extra_cuts=None):
        self.x=x;self.mass=mass;self.factors=factors;self.z=np.array(z,copy=True)
        g,w=np.polynomial.legendre.leggauss(order);self.g=(g+1)/2;self.w=w/2
        self.boundary_velocity=boundary_velocity;self.split_faces=0;self.nodes=0
        self.extra_cuts=extra_cuts
    def __call__(self,edge_id,edge,states,gen):
        from kinetic_factor_rhs import ALPHA
        from quadratic_velocity_reconstruction import velocity_samples
        i,j,a,b,n,length=edge
        def relative(point):
            vl=velocity_samples(states[i],self.factors[i],self.mass[i],point,ALPHA)[0][0]
            if j>=0:
                vr=velocity_samples(states[j],self.factors[j],self.mass[j],point,ALPHA)[0][0]
                d=np.linalg.norm(self.x[j]-self.x[i])
                wn=((self.x[j]-self.x[i])@gen[j]-(gen[j]-gen[i])@(point-self.x[i]))/d+(self.z[i]-self.z[j])/(2*d)
                return np.r_[vl@n-wn,vr@n-wn]
            wn=self.boundary_velocity(point)@n
            return np.r_[vl@n-wn,-(vl@n-wn)]
        left,middle,right=relative(a),relative((a+b)/2),relative(b)
        c2=2*(right+left-2*middle);c1=right-left-c2
        for t in (.25,.75):
            if not np.allclose(relative(a+t*(b-a)),left+t*c1+t*t*c2,atol=2e-12,rtol=2e-12):
                raise ValueError('normal sample speed is not quadratic')
        cuts=polynomial_cuts(left,middle,right)
        if self.extra_cuts is not None:
            cuts=np.unique(np.r_[cuts,self.extra_cuts(a,b)])
            if not np.isfinite(cuts).all() or cuts[0]!=0 or cuts[-1]!=1:
                raise ValueError('invalid extra face cuts')
        self.split_faces+=int(len(cuts)>2)
        nodes=[(lo+(hi-lo)*g,(hi-lo)*w) for lo,hi in zip(cuts[:-1],cuts[1:]) for g,w in zip(self.g,self.w)]
        self.nodes+=len(nodes);return nodes
