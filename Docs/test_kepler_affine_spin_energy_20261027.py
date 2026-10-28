"""Analytic square checks for polygon tensor and affine spin projection."""

import numpy as np

from audit_kepler_affine_spin_energy_20261027 import energies, tensor_from_faces


def main():
    a=np.array([[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]],dtype=float)
    b=np.roll(a,-1,axis=0)
    area,centroid,qxx,qxy,qyy=tensor_from_faces(
        a,b,np.zeros(4,dtype=int),np.array([0,4]),np.array([4]),np.array([1.]))
    np.testing.assert_allclose(area,[1.],rtol=0,atol=1e-14)
    np.testing.assert_allclose(centroid,[[0.,0.]],rtol=0,atol=1e-14)
    np.testing.assert_allclose([qxx[0],qxy[0],qyy[0]],[1/12,0,1/12],rtol=0,atol=1e-14)
    def ratio(axx,axy,ayx,ayy):
        spin,rot,total,res=energies(*(np.array([x]) for x in (axx,axy,ayx,ayy)),qxx,qxy,qyy)
        return spin[0],rot[0]/total[0],res[0]
    spin,part,res=ratio(0,-2,2,0)
    np.testing.assert_allclose([spin,part,res],[1/3,1,0],rtol=0,atol=1e-14)
    np.testing.assert_allclose(ratio(0,0,1,0)[1],.5,rtol=0,atol=1e-14)
    np.testing.assert_allclose(ratio(0,-1,-.5,0)[1],.1,rtol=0,atol=1e-14)
    print('PASS: square inertia, rigid rotation, simple shear and local Kepler ratio')


if __name__=='__main__':main()
