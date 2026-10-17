"""Small deterministic gates, including failure paths of the actual runner."""
import ast
from pathlib import Path
import unittest
import numpy as np
from stag_checks import signed_area, check_orientation, overlap_remap, atomic_substeps


class Gates(unittest.TestCase):
    def setUp(self):
        self.x = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
        self.old = self.x[[[0, 1, 2], [0, 2, 3]]]
        self.new = self.x[[[0, 1, 3], [1, 2, 3]]]

    def test_constant_pressure_density_jump_and_cycle(self):
        gamma = 5/3; p = 2.78e-5
        mass = np.array([.5, .005]); internal = np.full(2, p*.5/(gamma-1))
        for _ in range(20):
            mass, internal, residual = overlap_remap(self.old, self.new, mass, internal)
            np.testing.assert_allclose((gamma-1)*internal/.5, p, rtol=2e-15)
            np.testing.assert_allclose(mass.sum(), .505, rtol=2e-15)
            self.old, self.new = self.new, self.old

    def test_nonuniform_pressure_uses_area_weights(self):
        mass, internal, _ = overlap_remap(self.old, self.new, [.5, .005], [1., 3.])
        np.testing.assert_allclose(internal, [2., 2.])

    def test_missing_receiver_rejected(self):
        with self.assertRaisesRegex(RuntimeError, 'closure'):
            overlap_remap(self.old, self.new[:1], [.5, .005], [1., 1.])

    def test_duplicate_receiver_rejected(self):
        with self.assertRaisesRegex(RuntimeError, 'closure'):
            overlap_remap(self.old, [self.new[0], self.new[0]], [.5, .005], [1., 1.])

    def test_zero_coverage_rejected(self):
        with self.assertRaisesRegex(RuntimeError, 'closure'):
            overlap_remap(self.old, self.new+10, [.5, .005], [1., 1.])

    def test_orientation_of_actual_runner(self):
        # Import functions without executing the CLI simulation at module level.
        tree = ast.parse(Path(__file__).with_name('stagP_checked.py').read_text())
        selected = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name in ('zone_geom', 'all_geom')]
        ns = dict(np=np, check_orientation=check_orientation)
        exec(compile(ast.Module(body=selected, type_ignores=[]), '<geometry>', 'exec'), ns)
        key = tuple((i, (0, 0)) for i in (0, 1, 2))
        ref = {key: signed_area(self.x[:3])}
        ns['all_geom']({key: None}, self.x, ref)
        crossed = self.x.copy(); crossed[2, 1] = -1
        with self.assertRaisesRegex(RuntimeError, 'orientation'):
            ns['all_geom']({key: None}, crossed, ref)

    def test_partial_attempt_ledger_discarded(self):
        def advance(x, dt):
            if dt == 1 or (dt == .5 and x > 0):
                raise RuntimeError('injected late failure')
            return x+dt, dt
        x, ledger, rejected = atomic_substeps(0., 1., advance)
        self.assertEqual(x, 1.)
        self.assertEqual(sum(ledger), 1.)
        self.assertEqual(len(ledger), 4)
        self.assertEqual(len(rejected), 2)

    def test_terminal_rejection(self):
        def advance(x, dt):
            raise RuntimeError('injected permanent failure')
        with self.assertRaisesRegex(RuntimeError, 'macrostep failed'):
            atomic_substeps(0., 1., advance)

    def test_rayleigh_taper(self):
        for outer, expected in [(2.9, -.525980552), (3.6, -.144695453)]:
            r = (2.5+outer)/2; f = .5; fp = -.5*np.pi/(outer-2.5)
            kappa2 = (f*f+2*r*f*fp)/r**3
            self.assertAlmostEqual(kappa2, expected, places=8)


if __name__ == '__main__':
    unittest.main(verbosity=2)
