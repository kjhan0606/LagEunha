"""Check exact edge torque against independent two-point Gauss integration."""

import numpy as np


def main():
    poly=np.array([[0.,0.],[2.,0.],[.3,1.]],dtype=float)
    centroid=np.mean(poly,axis=0);g=np.array([2.,-1.]);p0=3.
    midpoint=0.;correction=0.;gauss=0.
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        e=b-a;normal=np.array([e[1],-e[0]]);c=.5*(a+b)
        cross=lambda x:x[0]*normal[1]-x[1]*normal[0]
        midpoint-=float((p0+g@c)*cross(c-centroid))
        correction-=float((g@e)*cross(e))/12.
        for t in (.5-.5/np.sqrt(3.),.5+.5/np.sqrt(3.)):
            point=a+t*e;gauss-=.5*float((p0+g@point)*cross(point-centroid))
    if abs(gauss)>1e-14 or abs(midpoint+correction-gauss)>1e-14 or abs(midpoint)<1e-3:
        raise ValueError(f'Incorrect torque integration: {midpoint=} {correction=} {gauss=}')
    print(f'PASS: exact linear torque zero; midpoint proxy={midpoint:.8g}')


if __name__=='__main__':main()
