import unittest
from unittest.mock import patch
import numpy as np
import test_dense_path_3d as fixture
from moment_evolution_3d import step


class IntegrationModeTests(unittest.TestCase):
    def test_both_stages_receive_mode_and_record_it(self):
        r=fixture.DensePath3DTests().record()
        diag=dict(r['diag_mid'],mass_flux=0.,ledger=np.zeros(7))
        response=(r['dr0'],r['dx0'],r['dw0'],diag)
        with patch('moment_evolution_3d.rhs',return_value=response) as mocked:
            out=step(None,0.,r['sites'],r['weights'],r['raw'],r['factors'],.001,
                lambda t:None,None,face_order=4,pressure_reconstruction='quadratic',
                return_record=True,face_integration='split')
        self.assertEqual(mocked.call_count,2)
        self.assertTrue(all(c.kwargs['face_integration']=='split' for c in mocked.call_args_list))
        record=out[-1]['dense_record']
        self.assertEqual(record['face_integration'],'split');self.assertEqual(record['face_order'],4)
        np.testing.assert_array_equal(out[2][:,0],r['raw'][:,0])
