import json,tempfile,unittest
from pathlib import Path
from factor_followup_batch import gate
class FollowupTests(unittest.TestCase):
    def test_gate(self):
        rows=[dict(kind=k,amplitude=a,dt=dt,status='completed',min_factor_R=1e-8,min_pressure=.5,final={'ledger_error':1e-15}) for k,a in [('wave',.001),('wave',.1),('shear',0.),('contact',0.)] for dt in (.001,.0005)]
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'results.jsonl'
            def write():p.write_text('\n'.join(json.dumps(r) for r in rows))
            write();self.assertTrue(gate(Path(d))[0])
            rows[0]['min_factor_R']=-1e-5;write();self.assertFalse(gate(Path(d))[0])
            rows[0]['min_factor_R']=0;rows[0]['kind']='unknown';write();self.assertFalse(gate(Path(d))[0])
if __name__=='__main__':unittest.main()
