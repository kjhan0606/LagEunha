"""Gauss nodes split at affine sample upwind zeros for a supplied weight rate.
The caller must re-solve the mass constraint and iterate until cuts and rates
agree. No velocity smoothing, energy repair, or prescribed site steering.
"""
import numpy as np

class AffineUpwindNodes:
    def __init__(self,x,mass,factors,z,order,boundary_velocity):
        self.x=x;self.mass=mass;self.factors=factors;self.z=np.array(z,copy=True)
        g,w=np.polynomial.legendre.leggauss(order);self.g=(g+1)/2;self.w=w/2
        self.boundary_velocity=boundary_velocity;self.split_faces=0;self.nodes=0

    def __call__(self,edge_id,edge,states,gen):
        from kinetic_factor_rhs import ALPHA
        i,j,a,b,n,length=edge
        def relative(point):
            s=states[i];vl=s['v']+s['B']@(point-s['c'])
            vl=vl+ALPHA@self.factors[i].T/np.sqrt(self.mass[i])
            if j>=0:
                d=np.linalg.norm(self.x[j]-self.x[i]);s=states[j]
                vr=s['v']+s['B']@(point-s['c'])+ALPHA@self.factors[j].T/np.sqrt(self.mass[j])
                wn=((self.x[j]-self.x[i])@gen[j]-(gen[j]-gen[i])@(point-self.x[i]))/d+(self.z[i]-self.z[j])/(2*d)
                return np.r_[vl@n-wn,vr@n-wn]
            # Boundary velocity must be affine along the face for exact cuts.
            wn=self.boundary_velocity(point)@n
            return np.r_[vl@n-wn,-(vl@n-wn)]
        left,right=relative(a),relative(b);delta=right-left
        # Verify affine contract rather than silently using endpoint roots.
        if not np.allclose(relative((a+b)/2),.5*(left+right),atol=1e-12,rtol=1e-12):
            raise ValueError('nonaffine normal sample speed')
        cuts=[0.,1.]
        for l,d in zip(left,delta):
            if d!=0:
                root=-l/d
                if 0<root<1:cuts.append(float(root))
        cuts=np.unique(cuts);self.split_faces+=int(len(cuts)>2)
        nodes=[(lo+(hi-lo)*g,(hi-lo)*w) for lo,hi in zip(cuts[:-1],cuts[1:]) for g,w in zip(self.g,self.w)]
        self.nodes+=len(nodes)
        return nodes
