import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from quadratic_pressure_reconstruction import prepare,value
from periodic2d_hllc import quadrature

class TensorPressureTests(unittest.TestCase):
    def setUp(self):
        import moment_ale_evolve as model
        self.addCleanup(setattr,model,'domain',model.domain)
        self.addCleanup(setattr,model,'boundary_velocity',model.boundary_velocity)
        self.enterContext(patch.dict('os.environ',LAGFORCE_PRESSURE_FIT='tensor_free_mean',LAGFORCE_PRESSURE_REFLECT_Y=''))

    def test_step_is_constant_transversely_and_keeps_cell_means(self):
        n=8;cells=[];states=[]
        for j in range(n):
            for i in range(n):
                cells.append(np.array([[i,j],[i+1,j],[i+1,j+1],[i,j+1]],dtype=float)/n)
                states.append(dict(P=1. if i<n/2 else .1,density=SimpleNamespace(area=1/n**2)))
        for degree in (2,4):
            with patch.dict('os.environ',LAGFORCE_PRESSURE_REFLECT_Y='0,1',LAGFORCE_PRESSURE_DEGREE=str(degree)):
                pp,d=prepare(states,cells)
            for k,(poly,meta) in enumerate(zip(cells,pp)):
                j,i=divmod(k,n);points=np.column_stack([np.full(5,(i+.4)/n),(j+np.linspace(0,1,5))/n])
                values=value(meta,points);np.testing.assert_allclose(values,values[0],rtol=0,atol=3e-12)
                x,w=quadrature(poly);self.assertAlmostEqual(float(w@value(meta,x)/w.sum()),states[k]['P'],places=12)
                self.assertGreater(value(meta,x).min(),0.)

    def test_global_quartic_and_thermal_energy_mean(self):
        from test_polynomial_pressure_reconstruction import data,pressure_and_gradient
        mass,state,L,cells,states=data()
        with patch.dict('os.environ',LAGFORCE_PRESSURE_DEGREE='4'):pp,d=prepare(states,cells)
        for poly,s,meta in zip(cells,states,pp):
            x,w=quadrature(poly);exact,_=pressure_and_gradient(x)
            np.testing.assert_allclose(value(meta,x),exact,rtol=0,atol=3e-8)
            self.assertAlmostEqual(float(w@value(meta,x)),s['P']*s['density'].area,places=12)

    def test_positivity_limiter_still_preserves_mean(self):
        from test_polynomial_pressure_reconstruction import data
        from polynomial_pressure_reconstruction import lower_bound
        _,_,_,cells,states=data()
        for i,s in enumerate(states):s['P']=1e-5 if i==len(states)//2 else .5+i%3
        with patch.dict('os.environ',LAGFORCE_PRESSURE_DEGREE='4'):pp,d=prepare(states,cells)
        for poly,s,meta in zip(cells,states,pp):
            self.assertGreaterEqual(lower_bound(meta,poly),.01*s['P']-2e-11)
            x,w=quadrature(poly);self.assertAlmostEqual(float(w@value(meta,x)),s['P']*s['density'].area,places=12)

    def test_pressure_fit_config_mismatch_is_rejected(self):
        import json
        from simple_wave_evolution import run_with_executor
        cfg=dict(mass_root_backend='legacy_dense',velocity_reconstruction='affine',pressure_reconstruction='quadratic_cell_mean')
        with patch.dict('os.environ',LAGFORCE_QUADRATIC_PRESSURE='1',LAGFORCE_PRESSURE_DEGREE='2',LAGFORCE_QUADRATIC_VELOCITY='0',LAGFORCE_VECTOR_KINETIC_ROOT='0'),patch('simple_wave_evolution.Path.read_bytes',return_value=json.dumps(cfg).encode()):
            with self.assertRaisesRegex(ValueError,'pressure fit configuration'):run_with_executor()

if __name__=='__main__':unittest.main()
