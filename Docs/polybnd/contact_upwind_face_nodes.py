"""Shared face quadrature split at both inferred cuts and velocity roots."""
import numpy as np
from quadratic_upwind_face_nodes import QuadraticUpwindNodes

class ContactUpwindNodes(QuadraticUpwindNodes):
    def __init__(self,*args,priors,**kwargs):
        super().__init__(*args,**kwargs);self.priors=priors

    def __call__(self,edge_id,edge,states,gen):
        i,j,*_=edge
        def cuts(a,b):
            values=[0.,1.]
            for k in (i,j):
                if k<0 or self.priors[k]['flat']:continue
                meta=self.priors[k];n=np.asarray(meta['normal']);den=float((b-a)@n)
                if den!=0:
                    t=(meta['offset']-a@n)/den
                    if 0<t<1:values.append(float(t))
            return np.unique(values)
        self.extra_cuts=cuts
        return super().__call__(edge_id,edge,states,gen)
