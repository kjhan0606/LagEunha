"""Exercise rejection paths of the native padding-reuse acceptance gate."""
from copy import deepcopy
import unittest
from unittest.mock import patch

import numpy as np

import audit_kepler_padding_reuse_pair_20261028 as audit


class AcceptanceGateTests(unittest.TestCase):
    def fixture(self):
        state = {field: np.ones(2, dtype=np.float64) for field in audit.FIELDS}
        state['x'][0] = 0.
        trace = {'steps': [('2.66', '23070', '1e-4')],
                 'tessellation': [('147456', '336000', '884790')],
                 'force': [('1', '147456', '884790')]}
        metrics = {'time': 2.66, 'wall_seconds': 10.,
                   'all_native_physics_checks': True,
                   'every_ledger_positive_and_injection_free': True,
                   'distinct_GPU_per_rank': True,
                   'parsed_GPU_assignments': [(0, 0, 2, 'A100')]}
        return state, trace, metrics

    def compare(self, a, b, exact=True):
        with patch.object(audit, 'load', side_effect=[a, b]):
            return audit.compare('control', 'candidate', exact)

    def test_identical_state_and_stage_counts_pass(self):
        a = self.fixture()
        self.assertTrue(self.compare(a, deepcopy(a))['passed'])

    def test_signed_zero_is_not_bitwise_equal(self):
        a = self.fixture()
        b = deepcopy(a)
        b[0]['x'][0] = -0.
        self.assertTrue(np.array_equal(a[0]['x'], b[0]['x']))
        result = self.compare(a, b)
        self.assertFalse(result['passed'])
        self.assertFalse(result['checks']['x_exact'])

    def test_stage_face_change_rejected_despite_equal_endpoint(self):
        a = self.fixture()
        b = deepcopy(a)
        b[1]['tessellation'][0] = ('147456', '336000', '884788')
        result = self.compare(a, b)
        self.assertFalse(result['passed'])
        self.assertFalse(result['checks']['stage_tessellation_counts_identical'])

    def test_different_rank_counts_allow_different_local_topology(self):
        a = self.fixture()
        b = deepcopy(a)
        b[1]['tessellation'][0] = ('73728', '177300', '442395')
        b[1]['force'][0] = ('1', '73728', '442395')
        self.assertTrue(self.compare(a, b, exact=False)['passed'])

    def test_invalid_device_mapping_rejected_despite_equal_endpoint(self):
        a = self.fixture()
        b = deepcopy(a)
        b[2]['distinct_GPU_per_rank'] = False
        self.assertFalse(self.compare(a, b)['passed'])

    def test_single_node_device_assignments(self):
        good = ('[GPU] rank 0 -> device 0/2: A100\n'
                '[GPU] rank 1 -> device 1/2: A100\n')
        self.assertTrue(audit.gpu_assignments(good, 2)[1])
        self.assertTrue(audit.gpu_assignments(good + good, 2)[1])
        self.assertFalse(audit.gpu_assignments(good.replace('device 1/2', 'device 0/2'), 2)[1])
        self.assertFalse(audit.gpu_assignments(good.splitlines()[0], 2)[1])
        self.assertFalse(audit.gpu_assignments('', 1)[1])


if __name__ == '__main__':
    unittest.main()
