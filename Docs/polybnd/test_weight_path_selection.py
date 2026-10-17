import unittest
from weight_path_followup import strict_predecessor_gate
class SelectionTests(unittest.TestCase):
    def test_strict_configuration_is_required(self):
        rows=[dict(kind=k,duration=h,rtol=1e-10,status='completed',sampled_path_gate=True) for k in ('shear','wave') for h in (.01,.05)]
        self.assertTrue(strict_predecessor_gate(rows+[dict(kind='shear',duration=.05,rtol=1e-8,status='completed',sampled_path_gate=False)]))
        self.assertFalse(strict_predecessor_gate(rows[:-1]))
        self.assertFalse(strict_predecessor_gate(rows+[rows[0]]))
        rows[-1]['sampled_path_gate']=False
        self.assertFalse(strict_predecessor_gate(rows))
if __name__=='__main__':unittest.main()
