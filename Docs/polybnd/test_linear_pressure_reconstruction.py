import unittest
import numpy as np
from adaptive_moment_batch import setup
from compatible_moment_rhs import rhs
import moment_ale_evolve as model
class LinearPressureTests(unittest.TestCase):
    def test_force_from_cell_average_pressure(self):
        m,y=setup('wave',3,0.);a=y[:-4].reshape(9,20)
        polygons,_=model.geometry(a[:,:2],a[:,2],0.);g=np.array([.2,-.1])
        for i,poly in enumerate(polygons):
            raw=np.r_[m[i],a[i,3:]];s=model.decode(poly,raw)
            a[i,19]=(.8+g@s['density'].origin)*s['density'].area/(model.GAMMA-1)
        dy,d=rhs(0.,y,m,linear_pressure=True)
        dr=dy[:-4].reshape(9,20)[:,3:]
        np.testing.assert_allclose(dr[:,6:8]/m[:,None],np.tile(-g,(9,1)),atol=2e-14,rtol=0)
        for i in range(9):
            qdotp=dr[i,8:12].reshape(2,2)
            q=y[:-4].reshape(9,20)[i,3:5]
            np.testing.assert_allclose(qdotp,-np.outer(q,g),atol=3e-15,rtol=0)
if __name__=='__main__':unittest.main()
