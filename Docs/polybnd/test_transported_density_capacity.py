import unittest
import numpy as np
from shapely.geometry import Polygon,LineString
from moment_aligned_density_audit import polygon_moments
from active_power_cells import power_cells_active
from transported_density_capacity import PolygonDensity,polygon_raw_many,mass_jacobian,solve_capacity,convex_line_intervals


class TransportedDensityCapacityTests(unittest.TestCase):
    def setUp(self):
        self.domain=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        left=np.array([[0.,0.],[.43,0.],[.43,1.],[0.,1.]])
        right=np.array([[.43,0.],[1.,0.],[1.,1.],[.43,1.]])
        self.field=PolygonDensity([left,right],[1.,100.])
        self.sites=np.array([[.22,.23],[.74,.19],[.27,.79],[.76,.75]])

    def test_local_polygon_moments_and_zero_area(self):
        cells=[Polygon(self.domain),Polygon(self.domain[::-1]),LineString([[0.,1.],[1.,.3]])]
        got=polygon_raw_many(np.asarray(cells,dtype=object),np.array([.47,.58]))
        np.testing.assert_allclose(got[:2],np.tile(polygon_moments(self.domain),(2,1)),rtol=1e-14,atol=1e-14)
        np.testing.assert_array_equal(got[2],np.zeros(7))
        mass,coverage,raw=self.field.integrate([cells[0]],moments=True)
        expected=sum((r*polygon_moments(v) for r,v in zip(self.field.densities,self.field.vertices)))
        np.testing.assert_allclose(raw[0],expected,rtol=1e-14,atol=1e-14)
        self.assertAlmostEqual(mass[0],expected[0],places=12);self.assertAlmostEqual(coverage[0],1.)

    def test_line_integrals_partition_shared_edges(self):
        value,gap=self.field.line_integral([0.,.5],[1.,.5])
        self.assertAlmostEqual(value,.43+57.,places=12);self.assertEqual(gap,0.)
        value,gap=self.field.line_integral([.43,0.],[.43,1.])
        self.assertAlmostEqual(value,50.5,places=12);self.assertEqual(gap,0.)
        tri=PolygonDensity([self.domain[[0,1,2]],self.domain[[0,2,3]]],[7.,7.])
        value,gap=tri.line_integral([0.,0.],[1.,1.])
        self.assertAlmostEqual(value,7*np.sqrt(2),places=12)

    def test_jacobian_matches_finite_differences_and_gauge(self):
        weights=np.array([.007,-.003,.001,0.]);cells=power_cells_active(self.sites,weights,self.domain)
        J,gap=mass_jacobian(self.sites,weights,cells,self.field);J=J.toarray()
        numerical=np.zeros_like(J);epsilon=1e-7
        for i in range(4):
            dw=np.zeros(4);dw[i]=epsilon
            plus=self.field.integrate(power_cells_active(self.sites,weights+dw,self.domain))[0]
            minus=self.field.integrate(power_cells_active(self.sites,weights-dw,self.domain))[0]
            numerical[:,i]=(plus-minus)/(2*epsilon)
        np.testing.assert_allclose(J,numerical,rtol=3e-8,atol=3e-8)
        np.testing.assert_allclose(J@np.ones(4),0.,atol=2e-14)
        self.assertGreater(np.linalg.eigvalsh(J)[1],0.);self.assertLess(gap,1e-12)

    def test_fixed_mass_affine_transport_uses_only_numerical_field(self):
        target=self.field.integrate(power_cells_active(self.sites,np.zeros(4),self.domain))[0]
        F=np.array([[1.,.22],[0.,1.]]);shift=np.array([.04,0.])
        moved=self.field.affine_image(F,shift);sites=self.sites@F.T+shift;domain=self.domain@F.T+shift
        weights,cells,history=solve_capacity(sites,np.zeros(4),domain,target,moved,rtol=5e-12)
        mass,coverage,_=moved.integrate(cells)
        np.testing.assert_allclose(mass,target,rtol=5e-12,atol=0.)
        self.assertGreater(len(history),1);self.assertEqual(weights[-1],0.)
        self.assertAlmostEqual(moved.total_mass,self.field.total_mass,places=12)
        with self.assertRaisesRegex(ValueError,'incompatible'):
            solve_capacity(sites,weights,domain,target*1.01,moved)

    def test_invalid_density_and_nonconvex_support_rejected(self):
        with self.assertRaises(ValueError):PolygonDensity([self.domain],[0.])
        with self.assertRaises(ValueError):PolygonDensity([[[0.,0.],[1.,0.],[.2,.2],[0.,1.]]],[1.])

    def test_numerical_velocity_is_inferred_and_restricted(self):
        from curved_contact_rhs import affine_state_from_position_moments
        from curved_transport_reference import numerical_affine_flow
        raw=polygon_moments(self.domain)[None,:]
        skeleton=np.zeros(24);skeleton[:2]=.5
        mass,state,_=affine_state_from_position_moments(skeleton,raw,raw)
        b,B,P,diag=numerical_affine_flow(state,mass,np.ones(1))
        np.testing.assert_allclose(b,[.2,0.],atol=1e-14)
        np.testing.assert_allclose(B,[[0.,1.1],[0.,0.]],atol=1e-14)
        self.assertAlmostEqual(P,.6,places=14)
        state[11]+=.01
        with self.assertRaises(ValueError):numerical_affine_flow(state,mass,np.ones(1))

    def test_compressive_affine_image_preserves_density_mass(self):
        F=np.array([[1.2,.1],[0.,.8]]);shift=np.array([2.,-3.])
        moved=self.field.affine_image(F,shift)
        cells=[Polygon(self.domain@F.T+shift)]
        mass,coverage,raw=moved.integrate(cells,moments=True)
        source=self.field.integrate([Polygon(self.domain)],moments=True)[2][0]
        expected_q=F@source[1:3]+shift*source[0]
        self.assertAlmostEqual(mass[0],source[0],places=11)
        np.testing.assert_allclose(raw[0,1:3],expected_q,rtol=1e-13,atol=1e-12)

    def test_thin_triangle_intervals_remain_local(self):
        for origin in (np.array([0.,0.]),np.array([.7,.8]),np.array([10.,-7.])):
            h=1e-3;thickness=1e-10
            vertices=origin+np.array([[[0.,0.],[h,0.],[h,thickness]],[[0.,0.],[h,thickness],[0.,thickness]]])
            start=origin+np.array([0.,thickness/2]);end=origin+np.array([h,thickness/2])
            lower,upper=convex_line_intervals(vertices,start,end)
            self.assertAlmostEqual(float(np.sum(upper-lower)),1.,places=12)
            field=PolygonDensity(vertices,[3.,3.]);value,gap=field.line_integral(start,end)
            self.assertLess(gap,1e-12);self.assertAlmostEqual(value,3*np.linalg.norm(end-start),places=14)

    def test_jacobian_enclosure_does_not_change_field_or_accept_real_gap(self):
        gap=2*np.finfo(float).eps
        left=np.array([[0.,0.],[.5-gap,0.],[.5-gap,1.],[0.,1.]])
        right=np.array([[.5+gap,0.],[1.,0.],[1.,1.],[.5+gap,1.]])
        field=PolygonDensity([left,right],[1.,100.]);mass=field.total_mass
        value,uncovered=field.line_integral([.5,0.],[.5,1.])
        self.assertAlmostEqual(value,50.5);self.assertLess(uncovered,2e-10)
        self.assertEqual(field.total_mass,mass);self.assertEqual(len(field.jacobian_enclosures),1)
        self.assertLess(field.jacobian_enclosures[0]['spatial_enclosure'],1e-13)
        right[:,0]=np.where(right[:,0]<.6,.50001,right[:,0])
        broken=PolygonDensity([left,right],[1.,100.])
        with self.assertRaisesRegex(RuntimeError,'uncovered capacity face'):
            broken.line_integral([.500001,0.],[.500001,1.])


if __name__=='__main__':unittest.main()
