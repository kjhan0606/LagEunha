import ast,os,unittest
from unittest.mock import patch
from pathlib import Path
import numpy as np
import p2_geometry_backend as backend
class GeometryBackendTests(unittest.TestCase):
    def test_signed_geometry_matches_reference(self):
        tree=ast.parse(Path(__file__).with_name('stagP2_perf_reference.py').read_text())
        nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('zone_geom','all_geom')]
        ns={'np':np};exec(compile(ast.Module(body=nodes,type_ignores=[]),'reference','exec'),ns)
        X=np.random.default_rng(417885).normal(size=(30,2));keys=[(i,i+1,i+2) for i in range(28)]+[(2,1,0),(0,0,0)]
        expected=ns['all_geom'](keys,X)
        for mode in ('numpy','openmp'):
            with patch.dict(os.environ, P2_GEOMETRY=mode):
                actual=backend.all_geom(keys,X)
            for a,b in zip(actual,expected):np.testing.assert_allclose(a,b,rtol=0,atol=1e-15)
        self.assertEqual(backend.threads_used(),int(os.environ['OMP_NUM_THREADS']))
if __name__=='__main__':unittest.main()
