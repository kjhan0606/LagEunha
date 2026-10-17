"""Static 2D weight-rate compatibility diagnostic, not a simulation.
At stationary generators, prescribed normal face speed s_ij requires
z_i-z_j=2 |x_j-x_i| s_ij. Edge differences must sum to zero on cycles.
"""
import json
import numpy as np

def probe():
    x=np.array([[0.,0.],[1.,0.],[.2,.8]])
    P=np.array([.6,.8,.7]);rho=np.array([1.,3.,2.]);Z=np.sqrt((5/3)*rho*P)
    edge=[(0,1),(1,2),(2,0)];A=np.zeros((3,3));h=[];speed=[]
    for k,(i,j) in enumerate(edge):
        s=(P[i]-P[j])/(Z[i]+Z[j]);speed.append(s)
        h.append(2*np.linalg.norm(x[j]-x[i])*s);A[k,i]=1;A[k,j]=-1
    h=np.array(h);z=np.r_[np.linalg.lstsq(A[:,:-1],h,rcond=None)[0],0.]
    # Least squares is a diagnostic only; it is NOT used to change face speeds.
    residual=A@z-h
    assert abs(sum(h))>1e-3
    assert np.max(abs(residual))>1e-3
    return dict(sites=x.tolist(),pressure=P.tolist(),density=rho.tolist(),
        oriented_edges=edge,acoustic_contact_speeds=speed,
        required_weight_differences=h.tolist(),cycle_sum=float(sum(h)),
        least_squares_weight_rates=z.tolist(),max_compatibility_residual=float(max(abs(residual))))

if __name__=='__main__': print(json.dumps(probe(),indent=2))
