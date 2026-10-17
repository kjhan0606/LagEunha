import tempfile,unittest
from pathlib import Path
import numpy as np
from dense_record_3d import save,load
import test_dense_path_3d as fixture
from moment_evolution_3d import dense_state

class DenseRecordTests(unittest.TestCase):
    def test_lossless_path_and_endpoint_roundtrip(self):
        record=fixture.DensePath3DTests().record();endpoint=dict(raw=record['raw'],factors=record['factors'],ledger=np.zeros(7))
        record.update(face_integration='split',face_order=4,pressure_reconstruction='quadratic')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'step.npz';save(path,record,endpoint);restored,end=load(path)
            for key in ('face_integration','face_order','pressure_reconstruction'):self.assertEqual(record[key],restored[key])
            for s in (0.,.37,1.):
                for a,b in zip(dense_state(record,s),dense_state(restored,s)):np.testing.assert_array_equal(a,b)
            for k in endpoint:np.testing.assert_array_equal(end[k],endpoint[k])
            with self.assertRaises(FileExistsError):save(path,record,endpoint)

if __name__=='__main__':unittest.main()
