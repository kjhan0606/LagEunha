import json,os,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import continue_enriched_kepler as driver


class EnrichedResumeTests(unittest.TestCase):
    def test_saved_trial_survives_interrupted_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);old=Path.cwd();os.chdir(root)
            try:
                y=np.zeros(24);y[:2]=3.;y[19]=1.
                source=dict(t=0.,state=y.tolist(),mass=[1.],factors=[np.eye(2).tolist()],potential_moments=[-1.])
                (root/'seed.json').write_text(json.dumps(source));(root/'cfg.json').write_text('{}');(root/'SHA256SUMS').write_text('fixture')
                diagnostic=dict(stage_mass_flux=0.,generator_velocity_defect=0.,hydro_record={})
                args=(str(root/'seed.json'),str(root/'cfg.json'),.01,.01)
                def audit(*args):
                    driver.atomic('report.json',dict(integrated_potential_boundary=0.,maximum_generator_defect=0.))
                    driver.atomic('nodes.json',[dict(flux=[0.]) for _ in range(17)])
                with patch.object(driver,'__file__',str(root/'driver.py')),patch.object(driver,'configure'):
                    with patch.object(driver,'step',return_value=(y,[np.eye(2)],np.array([-1.]),diagnostic)),patch.object(driver,'audit_path',side_effect=RuntimeError('interrupted')):
                        with self.assertRaisesRegex(RuntimeError,'interrupted'):driver.run(*args)
                    self.assertTrue((root/'step_00001/trial.json').exists());self.assertFalse((root/'checkpoint.json').exists())
                    with patch.object(driver,'step',side_effect=AssertionError('trial must not repeat')),patch.object(driver,'audit_path',side_effect=audit):driver.run(*args)
                    cp=json.loads((root/'checkpoint.json').read_text());self.assertEqual(cp['accepted_steps'],1);self.assertEqual(cp['t'],.01)
                    with patch.object(driver,'step',side_effect=AssertionError('completed prefix must not repeat')):driver.run(*args)
            finally:os.chdir(old)


if __name__=='__main__':unittest.main()
