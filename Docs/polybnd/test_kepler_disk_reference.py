import unittest
import numpy as np
from kepler_disk_reference import polar_rule,physical_raw,fields,frozen,CENTER,PRESSURE,GAMMA
from kepler_spline_gravity import acceleration
from periodic2d_hllc import quadrature


class KeplerDiskReferenceTests(unittest.TestCase):
    def test_evolved_near_origin_angular_layer(self):
        p=np.array([[2.7610171141540185,3.1673297731510233],[2.714978371692607,3.0365450452742033],[2.753853239319096,2.9694010801637436],[2.9999999953031566,2.999999975558807],[2.8497712232354653,3.1987376605574207]])
        u,v=quadrature(p)
        for order in (16,24):
            x,w=polar_rule(p,order,radii=(.25,.5,2.,2.5))
            np.testing.assert_allclose(np.r_[w.sum(),w@x,(x.T@(w[:,None]*x)).ravel()],np.r_[v.sum(),v@u,(u.T@(v[:,None]*u)).ravel()],atol=3e-12,rtol=0)

    def test_roundoff_vertex_loop_does_not_remove_polar_wedge(self):
        p=np.array([[2.9315588590434056,2.7301563470057153],[3.069052879290041,2.711494525885677],[3.1534767548590246,2.7969350755665667],[3.,3.],[3.,3.],[2.9999999999999996,3.],[2.9999999999999996,2.9999999999999996],[2.9027219642028057,2.7699118775096414]])
        x,w=polar_rule(p,24);u,v=quadrature(p)
        np.testing.assert_allclose(np.r_[w.sum(),w@x],np.r_[v.sum(),v@u],atol=2e-14,rtol=0)
        raw,_,d=physical_raw(p,24)
        self.assertAlmostEqual((raw[-1]-.5*(raw[13]+raw[16]))*(GAMMA-1)/d['volume'],PRESSURE,places=14)

    def test_genuinely_nonconvex_polygon_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'nonconvex'):
            polar_rule(CENTER+np.array([[0.,0.],[1.,0.],[.2,.2],[1.,1.],[0.,1.]]))

    def test_exact_repeated_vertices_preserve_polygon_integrals(self):
        p=np.array([[3.,3.],[3.4,3.],[3.3,3.5]])
        a,_,_=physical_raw(p,16);b,_,_=physical_raw(p[[0,0,1,1,2,2,0]],16)
        np.testing.assert_allclose(a,b,atol=2e-13,rtol=0)

    def test_polar_quadrature_matches_exact_polynomial_integrals(self):
        for poly in (np.array([[0.,0.],[6.,0.],[6.,6.],[0.,6.]]),np.array([[3.3,3.1],[3.8,3.2],[3.5,3.7]])):
            x,w=polar_rule(poly,20);u,v=quadrature(poly)
            self.assertGreater(w.min(),0.)
            np.testing.assert_allclose(np.r_[w.sum(),w@x,(x.T@(w[:,None]*x)).ravel()],np.r_[v.sum(),v@u,(u.T@(v[:,None]*u)).ravel()],atol=2e-11,rtol=0)

    def test_full_box_disk_mass_and_symmetry(self):
        poly=np.array([[0.,0.],[6.,0.],[6.,6.],[0.,6.]])
        raw,L,d=physical_raw(poly,20);expected=.01*36+.99*np.pi*(4-.25)
        self.assertAlmostEqual(raw[0],expected,places=11)
        np.testing.assert_allclose(raw[1:3],raw[0]*CENTER,atol=2e-11,rtol=0)
        np.testing.assert_allclose(raw[7:9],0.,atol=2e-12,rtol=0)
        self.assertAlmostEqual((raw[-1]-.5*(raw[13]+raw[16]))*(GAMMA-1)/d['volume'],PRESSURE,places=14)

    def test_interface_cell_reference_converges_independently(self):
        for lo,hi in ((3.3,3.7),(4.8,5.2)):
            poly=np.array([[lo,2.85],[hi,2.85],[hi,3.15],[lo,3.15]])
            a,_,_=physical_raw(poly,16);b,_,_=physical_raw(poly,24)
            np.testing.assert_allclose(a,b,atol=3e-11,rtol=0)

    def test_circular_velocity_balances_canonical_force(self):
        x=CENTER+np.column_stack([np.array([.01,.2,.5,1.,2.,2.6]),np.zeros(6)])
        rho,v,P=fields(x);g=acceleration(x,center=CENTER,rc=.5);r=np.linalg.norm(x-CENTER,axis=1)
        np.testing.assert_allclose(np.sum(v*v,axis=1)/r,-g[:,0],atol=2e-14,rtol=0)
        np.testing.assert_array_equal(frozen(0.,x),[False,False,False,False,False,True]);np.testing.assert_array_equal(P,np.full(6,PRESSURE))

if __name__=='__main__':unittest.main()
