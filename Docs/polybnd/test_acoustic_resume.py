import json,os,tempfile,unittest
from pathlib import Path
import numpy as np
from check_3d_acoustic_evolution import load_resume

class AcousticResumeTests(unittest.TestCase):
    def test_original_mass_and_configuration_are_required(self):
        with tempfile.TemporaryDirectory() as directory:
            old=os.getcwd();os.chdir(directory)
            try:
                raw=np.zeros((1,35));raw[0,0]=1.
                ic=dict(n=1,dt=.001,end=.02,pressure_model='quadratic',mass=[1.],raw=raw.tolist())
                cp=dict(ic,t=.005,sites=[[.5,.5,.5]],weights=[0.],factors=np.zeros((1,3,3)).tolist(),ledger=[0.]*7)
                Path('initial_state.json').write_text(json.dumps(ic));Path('checkpoint.json').write_text(json.dumps(cp))
                saved,initial=load_resume(1,.001,.02,'quadratic');self.assertEqual(saved['t'],.005)
                with self.assertRaisesRegex(ValueError,'configuration'):load_resume(1,.002,.02,'quadratic')
                cp['mass']=[2.];Path('checkpoint.json').write_text(json.dumps(cp))
                with self.assertRaisesRegex(ValueError,'masses'):load_resume(1,.001,.02,'quadratic')
                cp['mass']=[1.];cp['factors'][0][0][0]=float('nan');Path('checkpoint.json').write_text(json.dumps(cp))
                with self.assertRaisesRegex(ValueError,'nonfinite'):load_resume(1,.001,.02,'quadratic')
            finally:os.chdir(old)

if __name__=='__main__':unittest.main()
