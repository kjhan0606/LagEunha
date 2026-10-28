import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from audit_kepler_vertex_halo_pair_20261029 import EXPECTED, launch_checks, path_checks


def log(reuse=True, alias=True):
    text=('[VERTEX_HALO] force=1 actual_rectangle=[-1.25,7.25]x[-1.25,7.25] W0_only=1\n'
          '[VERTEX_HALO] final_volume=1 actual_rectangle=[-1.25,7.25]x[-1.25,7.25] W0_only=1\n'
          '[VERTEX_HALO] recertified=180 incomplete=0\n')
    if reuse: text+='[TESS_PADDING_REUSE] rank=0 np=147456 npad=180000 reuse=1 refresh_validation=1\n'
    if alias: text+='[LOCAL_PERIODIC_SYNC] rank=0 device alias refresh npad=180000\n'
    return text


class ComboGateTests(unittest.TestCase):
    def test_executed_combination_passes(self):
        checks,_=path_checks('vertex20_reuse_local_sync',log())
        self.assertTrue(all(checks.values()))

    def test_noop_reuse_or_alias_rejected(self):
        for text in [log(reuse=False),log(alias=False),log()+'[HALO_SYNC] host refresh\n']:
            checks,_=path_checks('vertex20_reuse_local_sync',text)
            self.assertFalse(all(checks.values()))

    def test_wrong_domain_or_incomplete_certificate_rejected(self):
        for text in [log().replace('7.25','7.5'),log().replace('incomplete=0','incomplete=1')]:
            checks,_=path_checks('vertex20_reuse_local_sync',text)
            self.assertFalse(all(checks.values()))

    def contexts(self,root):
        keys=['GFS_GPU_VERTEX_HALO_CERT','GFS_KEPLER_HALO_CELLS','GFS_REUSE_TESS_PADDING','GFS_GPU_LOCAL_PERIODIC_SYNC']
        for name,flags in EXPECTED.items():
            folder=root/name;folder.mkdir()
            value=dict(zip(keys,map(str,flags)))
            value.update(SLURM_JOB_ID='12345',SLURM_NODELIST='syn103',CUDA_VISIBLE_DEVICES='0')
            (folder/'launch_context.json').write_text(json.dumps(value))

    def test_allocation_and_seven_case_completeness(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);self.contexts(root)
            self.assertTrue(all(launch_checks(root).values()))
            p=root/'vertex20_reuse_local_sync/launch_context.json'
            value=json.loads(p.read_text());value['SLURM_JOB_ID']='54321';p.write_text(json.dumps(value))
            self.assertFalse(all(launch_checks(root).values()))
            p.unlink()
            self.assertFalse(all(launch_checks(root).values()))


if __name__=='__main__':unittest.main()
