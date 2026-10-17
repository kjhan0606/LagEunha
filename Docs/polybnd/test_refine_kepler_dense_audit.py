import unittest
from refine_kepler_dense_audit import interleave_seed


class RefineKeplerTests(unittest.TestCase):
    def test_all_seed_nodes_reused_at_exact_positions(self):
        seed=[{'value':i} for i in range(17)]
        for intervals in (32,64):
            result=interleave_seed(seed,intervals)
            self.assertEqual(len(result),17)
            for i,row in enumerate(seed):self.assertIs(result[str(i*intervals//16)],row)

    def test_invalid_seed_or_lattice(self):
        for count,n in [(16,32),(17,16),(17,33)]:
            with self.assertRaises(ValueError):interleave_seed([{}]*count,n)
