import hashlib,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from adaptive_moment_batch import D0
from curved_contact_rhs import affine_state_from_position_moments
from moment_aligned_density_audit import polygon_moments
from saved_curved_contact import initialize_saved,advected_boundary,reference_raw
import predictive_evolution_batch as batch
import moment_ale_evolve as model


class SavedCurvedContactTests(unittest.TestCase):
    def fixture(self,directory):
        theta=np.arange(128)*2*np.pi/128
        boundary=.5+.22*np.column_stack([np.cos(theta),np.sin(theta)])
        geo=polygon_moments(D0)[None,:];raw=geo+99*polygon_moments(boundary)[None,:]
        skeleton=np.zeros(24);skeleton[:2]=.5
        mass,state,factors=affine_state_from_position_moments(skeleton,raw,geo)
        source=dict(n=1,t=0.,mass=mass.tolist(),state=state.tolist(),factors=[v.tolist() for v in factors])
        reference=dict(raw=raw.tolist(),geometry=geo.tolist(),reference_boundary=boundary.tolist())
        a=Path(directory)/'initial.json';b=Path(directory)/'reference.json'
        a.write_text(json.dumps(source));b.write_text(json.dumps(reference))
        cfg=dict(kind='saved_curved_contact',n=1,initial_state=str(a),reference=str(b),
            initial_state_sha256=hashlib.sha256(a.read_bytes()).hexdigest(),
            reference_sha256=hashlib.sha256(b.read_bytes()).hexdigest())
        return cfg,source,reference

    def test_dispatch_preserves_physical_mass_without_shear_ic(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg,source,reference=self.fixture(directory)
            with patch.object(batch,'setup',side_effect=AssertionError('must not regenerate another IC')):
                mass,state=batch.initialize_case(cfg)
            np.testing.assert_array_equal(mass,source['mass'])
            np.testing.assert_array_equal(state,source['state'])
            self.assertAlmostEqual(mass.sum(),1+99*64*.22**2*np.sin(2*np.pi/128),places=10)
            self.assertGreater(abs(mass.sum()-50.5),30.)

    def test_hash_and_mass_corruption_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg,source,reference=self.fixture(directory)
            with self.assertRaisesRegex(ValueError,'hash mismatch'):
                initialize_saved(dict(cfg,initial_state_sha256='0'*64))
            reference['raw'][0][0]*=2
            Path(cfg['reference']).write_text(json.dumps(reference))
            cfg['reference_sha256']=hashlib.sha256(Path(cfg['reference']).read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError,'mass mismatch'):initialize_saved(cfg)

    def test_exact_advected_polygon_moments(self):
        with tempfile.TemporaryDirectory() as directory:
            cfg,source,reference=self.fixture(directory);initialize_saved(cfg)
            for t in (0.,.00011,.02):
                raw,geo,_=reference_raw([model.domain(t)],reference['reference_boundary'],t)
                expected=polygon_moments(model.domain(t))+99*polygon_moments(advected_boundary(reference['reference_boundary'],t))
                np.testing.assert_allclose(raw[0],expected,rtol=2e-13,atol=2e-13)
                self.assertAlmostEqual(raw[0,0],source['mass'][0],places=11)

    def test_path_workers_install_only_identical_boundary_law(self):
        from adaptive_moment_batch import setup
        setup('shear_contact',1,100.)
        domains=[model.domain(t) for t in (0.,.3)];vel=model.boundary_velocity(np.array([.2,.7]))
        with patch.object(batch,'setup',side_effect=AssertionError('unnecessary interior IC')):
            batch.initialize_path_worker(64,100.,True)
        for t,domain in zip((0.,.3),domains):np.testing.assert_array_equal(model.domain(t),domain)
        np.testing.assert_array_equal(model.boundary_velocity(np.array([.2,.7])),vel)


if __name__=='__main__':unittest.main()
