import unittest
import numpy as np
from smooth_path_quadrature import node_indices,upgrade_nine_samples,assess_values


class SmoothPathQuadratureTests(unittest.TestCase):
    def config(self):return dict(path_mass_budget=1.,path_quadrature_budget=1.,path_rate_budget=1.,generator_path_budget=1.)

    def test_nested_indices_and_lossless_legacy_upgrade(self):
        old={str(j):dict(flux=[float(j)]) for j in range(9)};new=upgrade_nine_samples(old)
        self.assertEqual([new[str(j)] for j in node_indices(9)],list(old.values()))
        for a,b in zip((9,17,33,65),(17,33,65,129)):self.assertTrue(set(node_indices(a))<=set(node_indices(b)))
        with self.assertRaises(ValueError):upgrade_nine_samples({'16':{}})

    def test_quartic_signed_and_absolute_quadrature_order(self):
        deltas=[]
        for count in (9,17,33):
            values=[dict(flux=[(j/128)**4],diag=dict(generator_defect=0.,min_pressure=.6,min_factor_R=0.)) for j in node_indices(count)]
            r=assess_values(values,[1.],.1,.1,self.config(),True)
            deltas.append(r['quad_delta'][0]);self.assertAlmostEqual(r['absolute_quad_delta'][0],r['quad_delta'][0])
        self.assertAlmostEqual(deltas[0]/deltas[1],16.,places=8);self.assertAlmostEqual(deltas[1]/deltas[2],16.,places=7)

    def test_absolute_integral_quadrature_is_checked_separately(self):
        values=[dict(flux=[(-1.)**j*(j+.3)],diag=dict(generator_defect=0.,min_pressure=.6,min_factor_R=0.)) for j in range(9)]
        a=assess_values(values,[1.],.01,.05,self.config(),False);b=assess_values(values,[1.],.01,.05,self.config(),True)
        self.assertGreaterEqual(b['quadrature_ratio'],a['quadrature_ratio'])
        self.assertGreaterEqual(b['ratio'],a['ratio'])



class ExtendedLatticeTests(unittest.TestCase):
    def test_finer_nodes_preserve_old_coordinates(self):
        from smooth_path_quadrature import node_indices
        for count in [9,17,33,65,129]:
            self.assertEqual(node_indices(count,1024),[8*j for j in node_indices(count)])
        self.assertEqual(len(node_indices(1025,1024)),1025)
        with self.assertRaises(ValueError):node_indices(257,128)

    def test_saved_samples_reindexed_without_recomputation(self):
        import tempfile,json
        from pathlib import Path
        from adaptive_coupled_mass import atomic
        from smooth_path_quadrature import ensure_sample_lattice
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);trial=dict(record=dict(example='saved'),sample_lattice=128)
            atomic(p/'pending.json',trial);atomic(p/'samples.json',{'0':dict(flux=[1.]),'64':dict(flux=[2.]),'128':dict(flux=[3.])})
            trial,samples=ensure_sample_lattice(trial,p/'pending.json',p/'samples.json',1024)
            self.assertEqual(samples,{'0':dict(flux=[1.]),'512':dict(flux=[2.]),'1024':dict(flux=[3.])})
            again,cached=ensure_sample_lattice(trial,p/'pending.json',p/'samples.json',1024)
            self.assertEqual(samples,cached)

    def test_interrupted_reindex_recovers_from_original_archive(self):
        import tempfile,hashlib
        from pathlib import Path
        from adaptive_coupled_mass import atomic
        from smooth_path_quadrature import ensure_sample_lattice
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);archive=p/'old.json';atomic(archive,{'64':dict(flux=[2.])})
            trial=dict(record=dict(example='saved'),sample_lattice=128,sample_lattice_upgrade=dict(source=128,target=1024,archive=archive.name,sha256=hashlib.sha256(archive.read_bytes()).hexdigest()))
            atomic(p/'pending.json',trial);atomic(p/'samples.json',{'512':dict(flux=[999.])})
            recovered,samples=ensure_sample_lattice(trial,p/'pending.json',p/'samples.json',1024)
            self.assertEqual(samples,{'512':dict(flux=[2.])});self.assertNotIn('sample_lattice_upgrade',recovered)

if __name__=='__main__':unittest.main()
