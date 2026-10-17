import hashlib,json,tempfile,unittest
from pathlib import Path
from extend_validated_wave import prepare


class ExtensionTests(unittest.TestCase):
    def fixture(self,root):
        source=root/'source';source.mkdir();cfg=dict(end=.01,validation_end=10.,n=32)
        payload=json.dumps(cfg).encode();identity=hashlib.sha256(payload).hexdigest()
        initial=dict(t=0.,mass=[1.],state=[2.],factors=[[3.]],config_sha256=identity)
        cp=dict(initial,t=.01,steps=12,rejected=3,next_dt=.001,absolute=[.000001],peak={'rate':.00001})
        report=dict(t=.01,config=cfg,completed=True,numerical_passed=True,physical_accuracy_passed=True,path_budget_duration=10.)
        (source/'simple_wave_config.json').write_bytes(payload)
        for name,value in [('initial_state.json',initial),('checkpoint.json',cp),('simple_wave_report.json',report)]:
            (source/name).write_text(json.dumps(value))
        (source/'execution.txt').write_text('exit=0\n');return source,initial,cp

    def test_preserves_every_state_and_budget_field(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source,initial,cp=self.fixture(root);before={p.name:p.read_bytes() for p in source.iterdir()}
            prepare(source,root/'next',10.)
            for name,old in [('initial_state.json',initial),('checkpoint.json',cp)]:
                new=json.loads((root/'next'/name).read_text());new['config_sha256']=old['config_sha256'];self.assertEqual(new,old)
            self.assertEqual(before,{p.name:p.read_bytes() for p in source.iterdir()})
            with self.assertRaises(FileExistsError):prepare(source,root/'next',10.)

    def test_failure_pending_and_budget_changes_rejected_before_writes(self):
        for kind in ('physical','numerical','duration','pending','horizon','identity','receipt'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);source,_,_=self.fixture(root);p=source/'simple_wave_report.json';r=json.loads(p.read_text());end=10.
                if kind in ('physical','numerical'):r[{'physical':'physical_accuracy_passed','numerical':'numerical_passed'}[kind]]=False
                if kind=='duration':r['path_budget_duration']=.01
                if kind=='pending':(source/'pending.json').write_text('{}')
                if kind=='horizon':end=11.
                if kind=='identity':(source/'simple_wave_config.json').write_text('{}')
                if kind=='receipt':(source/'execution.txt').write_text('exit=124\n')
                p.write_text(json.dumps(r))
                with self.assertRaises(ValueError):prepare(source,root/'next',end)
                self.assertFalse((root/'next').exists())


if __name__=='__main__':unittest.main()
