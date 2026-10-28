import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from moment_state_3d import pack
from power_cells_3d import face_quadrature
from kinetic_moment_rhs_3d import rhs
from split_face_mass_3d import solve


class SharedFaceRefinementTests(unittest.TestCase):
    def fixture(self):
        sites=np.array([[.25,.5,.5],[.75,.5,.5]])
        vertices=np.array([[.5,0,0],[.5,1,0],[.5,1,1],[.5,0,1]])
        n=np.array([1.,0,0]);points,weights=face_quadrature(vertices,n,4)
        face=dict(i=0,j=1,vertices=vertices,normal=n,points=points,weights=weights)
        geo=SimpleNamespace(build=lambda *a,**k:([None,None],[face],{}))
        raw=[];states=[]
        for m,c,vx in zip([1.,2.],sites,[.2,-.1]):
            v=np.array([vx,0.,0]);S=.1*np.eye(3);R=.01*np.eye(3)
            raw.append(pack(m,m*c,S+m*np.outer(c,c),m*v,m*np.outer(c,v),m*np.outer(v,v)+R,10.))
            density=SimpleNamespace(value=lambda p,m=m:np.full(len(p),m),volume=1.,origin=c)
            states.append(dict(m=m,c=c,v=v,S=S,B=np.zeros((3,3)),R=R,P=1.,density=density))
        return geo,sites,np.asarray(raw),np.tile(.1*np.eye(3),(2,1,1)),states

    def test_refinement_resolves_new_root_and_preserves_inputs(self):
        geo,x,raw,L,states=self.fixture();before=raw.copy();calls=[]
        def root(*args,**kwargs):
            calls.append(kwargs['order'])
            if kwargs['order']==4:return np.zeros(2),dict(iterations=0,residual=0.)
            return solve(*args,**kwargs)
        decode=lambda p,r:states[int(r[0])-1]
        with patch('kinetic_moment_rhs_3d.decode',side_effect=decode),patch('split_face_mass_3d.decode',side_effect=decode),patch('split_face_mass_3d.solve',side_effect=root):
            dr,gen,z,d=rhs(geo,x,np.zeros(2),raw,L,face_integration='split')
        self.assertEqual(calls,[4,6]);self.assertEqual(d['effective_face_order'],6)
        self.assertLess(d['mass_flux'],2e-13);self.assertGreater(abs(z[0]),0.)
        np.testing.assert_array_equal(raw,before)
        np.testing.assert_array_equal(gen,np.array([s['v'] for s in states]))
        self.assertEqual(d['face_quadrature_retries'][0]['order'],4)

    def test_exhaustion_fails_without_relaxation(self):
        geo,x,raw,L,states=self.fixture();decode=lambda p,r:states[int(r[0])-1]
        with patch('kinetic_moment_rhs_3d.decode',side_effect=decode),patch('split_face_mass_3d.solve',return_value=(np.zeros(2),dict(iterations=0,residual=0.))):
            with self.assertRaisesRegex(RuntimeError,'at order 4'):
                rhs(geo,x,np.zeros(2),raw,L,face_integration='split',max_face_order=4)


if __name__=='__main__':unittest.main()
