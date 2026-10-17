import unittest
from unittest.mock import patch
import numpy as np
from shapely.geometry import box
from compare_physical_overlap import density_envelope,velocity_envelope,pressure_envelope,intersections


class PhysicalOverlapTests(unittest.TestCase):
    def test_partition_strip_area_is_the_region_with_different_cell_labels(self):
        a=[box(0.,0.,.5,1.),box(.5,0.,1.,1.)];b=[box(0.,0.,.51,1.),box(.51,0.,1.,1.)];common,strip=intersections(a,b)
        self.assertAlmostEqual(strip,.01,places=14);self.assertAlmostEqual(sum(p.area for p in common)+strip,1.,places=14)

    def test_density_envelope_contains_exponential_quadratic_values(self):
        class Density:
            origin=np.array([.5,.5]);scale=1.;mass=2.;logZ=.2;theta=np.array([.2,-.1,.3,-.2,-.1])
        density=Density();poly=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]]);points=np.random.default_rng(8).random((500,2));z=points-.5;F=np.column_stack([z,z[:,0]**2,z[:,0]*z[:,1],z[:,1]**2]);rho=2*np.exp(F@density.theta-.2)
        self.assertLessEqual(float(rho.max()),density_envelope(density,poly))

    def test_velocity_and_pressure_envelopes_contain_actual_reconstruction(self):
        from test_quadratic_velocity_reconstruction import fixture
        from quadratic_velocity_reconstruction import prepare,correction
        from quadratic_pressure_reconstruction import prepare as pressure_prepare,value
        mass,state,L,cells,edges,states=fixture()
        with patch.dict('os.environ',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VECTOR_QUADRATIC_BUDGET='1',LAGFORCE_PRESSURE_DEGREE='2'):
            prepare(states,cells,edges,mass,L);polys,_=pressure_prepare(states,cells)
        for s,p,meta in zip(states,cells,polys):
            points=np.vstack([p,.5*(p+np.roll(p,-1,axis=0)),p.mean(0)]);v=s['v']+(points-s['c'])@s['B'].T+correction(s,points)
            self.assertLessEqual(float(np.linalg.norm(v,axis=1).max()),velocity_envelope(s,p));self.assertLessEqual(float(abs(value(meta,points)).max()),pressure_envelope(meta,p))
