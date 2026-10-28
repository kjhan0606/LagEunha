"""One original-seed midpoint provider replay, after contact fine completion.

Uses original decode/evaluate AST prefixes without rewriting arithmetic. No
update, new trajectory, backend shadow, warmup, or tuning replay is executed.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np

base = Path('/gpfs/kjhan/LagForce')
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--preflight-only', action='store_true', help='Resolve imports, inputs, hashes and function definitions; never call provider')
args = parser.parse_args()
if not args.preflight_only:
    fine = json.loads((base/'Codex_contact_prior_evolution_dt5e5_20261018/summary.json').read_text())
    assert fine['status'] == 'completed' and fine['completed_time'] == .001
modules = base/'Codex_provider_profile_modules_20261018'
sys.path[:0] = [str(modules), str(base/'Codex_compiled_midpoint_modules_20261018'),
               str(base/'Codex_quadratic_thermal_audit_modules_20261018'),
               str(base/'Codex_quadratic_split_modules_20261018'),
               str(base/'Codex_coupled_relaxation_evolution_20261018'), str(base)]
import provider_profile_clock as clock
import moment_ale_evolve as model
import fixedmass_hllc
import quadratic_velocity_reconstruction as quad
import frozen_weight_root as frozen_root
import vector_kinetic_root as vector_root
import euler_star_pressure
from kinetic_factor_rhs import rhs, ALPHA
from contact_prior_density import ContactPriorDensity, local_contact_priors
from contact_upwind_face_nodes import ContactUpwindNodes
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from kinetic_face_capture import capture
from moment_gravity_source import source
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read

flags = dict(LAGFORCE_QUADRATIC_VELOCITY='1', LAGFORCE_QUADRATIC_BUDGET='projected',
             LAGFORCE_VELOCITY_FIT='radial_fixed_mean', LAGFORCE_VELOCITY_STENCIL='area',
             LAGFORCE_QUADRATIC_PRESSURE='0', LAGFORCE_KINETIC_SMOOTHING='0',
             LAGFORCE_PRESSURE_RIEMANN='exact_euler_star', LAGFORCE_VECTOR_KINETIC_ROOT='0',
             LAGFORCE_SHARED_FACE_VELOCITY='1', LAGFORCE_CAPTURE_PRESSURE_RATES='1')
os.environ.update(flags)
# Original runner did not set this flag. A non-default inherited flag requires
# provenance resolution before this single replay can be used.
assert os.environ.get('LAGFORCE_VECTOR_MOMENT_FLUX', '0') == '0'
assert not os.environ.get('LAGFORCE_ANGULAR_FACE_CELLS', '')
model.GAMMA = float(np.float32(1.666667))
fixedmass_hllc.GAMMA = model.GAMMA
saved_root = base/'Codex_compiled_midpoint_once_20261018'
original_manifest = json.loads((saved_root/'manifest.json').read_text())
profile_manifest = json.loads((modules/'profile_manifest.json').read_text())
runner = base/'Codex_compiled_midpoint_modules_20261018/run_compiled_midpoint_once.py'
assert hashlib.sha256(runner.read_bytes()).hexdigest() == original_manifest['runner_sha256']
for name, meta in profile_manifest['modules'].items():
    assert hashlib.sha256(Path(meta['original_path']).read_bytes()).hexdigest() == meta['original_sha256']
    assert hashlib.sha256((modules/(name+'.py')).read_bytes()).hexdigest() == meta['instrumented_sha256']


class RemoveTimers(ast.NodeTransformer):
    def visit_ImportFrom(self, node):
        return None if node.module == 'provider_profile_clock' else node

    def visit_With(self, node):
        if len(node.items) == 1 and isinstance(node.items[0].context_expr, ast.Call):
            call = node.items[0].context_expr
            if isinstance(call.func, ast.Name) and call.func.id == '_profile_section':
                return [self.visit(n) for n in node.body]
        return self.generic_visit(node)


for name in profile_manifest['instrumented_blocks']:
    original = ast.parse(Path(profile_manifest['modules'][name]['original_path']).read_text())
    instrumented = ast.parse((modules/(name+'.py')).read_text())
    stripped = RemoveTimers().visit(instrumented)
    assert ast.dump(original, include_attributes=False) == ast.dump(stripped, include_attributes=False)

mid = np.load(saved_root/'midpoint/inputs.npz')
stage0 = np.load(saved_root/'stage0/inputs.npz')
expected = np.load(saved_root/'midpoint/reference.npz')
mass, ids, frozen = mid['mass'], mid['ids'], mid['frozen'].astype(bool)
n = len(ids)
y, L, seed = mid['state'].copy(), mid['factors'].copy(), stage0['weight_rates'].copy()
ic, t0 = read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
assert t0 == 0.
spec = importlib.util.spec_from_file_location('clip', base/'Codex_weighted_compression_geometry_20261017/audit.py')
clip = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clip)
builder = EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3, ic['y']-3]), ids, clip.clip)
model.geometry = builder
model.boundary_velocity = lambda point: np.zeros(2)

# Wrapping bound aliases preserves original function objects and caches.
EmbeddedLaguerreGeometry.__call__ = clock.wrap(EmbeddedLaguerreGeometry.__call__, 'geometry.build')
ContactPriorDensity.__init__ = clock.wrap(ContactPriorDensity.__init__, 'density.fit')
ContactPriorDensity.value_trace = clock.wrap(ContactPriorDensity.value_trace, 'density.face_trace')
ContactUpwindNodes.__call__ = clock.wrap(ContactUpwindNodes.__call__, 'face.quadrature_split')
local_contact_priors = clock.wrap(local_contact_priors, 'density.local_contact_priors')
quad.prepare = clock.wrap(quad.prepare, 'velocity.prepare')
quad.velocity_samples_and_mean = clock.wrap(quad.velocity_samples_and_mean, 'velocity.face_samples')
quad.correction = clock.wrap(quad.correction, 'velocity.correction')
model.decode = clock.wrap(model.decode, 'moment.decode')
euler_star_pressure.star = clock.wrap(euler_star_pressure.star, 'pressure.exact_star')
frozen_root.record_arrays = clock.wrap(frozen_root.record_arrays, 'massroot.record_arrays')
frozen_root.evaluate = clock.wrap(frozen_root.evaluate, 'massroot.residual_and_optional_jacobian')
frozen_root.spsolve = clock.wrap(frozen_root.spsolve, 'massroot.spsolve')
frozen_root.lsmr = clock.wrap(frozen_root.lsmr, 'massroot.lsmr')
frozen_root.solve = clock.wrap(frozen_root.solve, 'massroot.solve')
source = clock.wrap(source, 'gravity.raw_source')
capture = clock.wrap(capture, 'diagnostic.face_capture')
rhs = clock.wrap(rhs, 'rhs.total')
events = []
progress = {}


def emit(row):
    events.append(row)
    with (args.output/'events.jsonl').open('a') as stream:
        stream.write(json.dumps(row)+'\n')
    print(json.dumps(row), flush=True)


def sha_array(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


# Use the original definitions. evaluate ends immediately after its original
# gravity/moment-defect block, before reference/backend shadows and file dumps.
tree = ast.parse(runner.read_text())
decode_fn = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == 'decode')
evaluate_fn = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == 'evaluate')
cut = next(k for k, node in enumerate(evaluate_fn.body)
           if ast.get_source_segment(runner.read_text(), node).startswith('pv=np.array'))
evaluate_fn.body = evaluate_fn.body[:cut]+[ast.parse('return dy,diag,z,gravity,moment_defect,len(records),iteration+1').body[0]]
for index, node in enumerate(evaluate_fn.body):
    text = ast.get_source_segment(runner.read_text(), node) if hasattr(node, 'end_lineno') else ''
    if isinstance(node, ast.For) and text.startswith('for i,s in enumerate(states):'):
        context = ast.Call(func=ast.Attribute(value=ast.Name(id='clock', ctx=ast.Load()), attr='section', ctx=ast.Load()),
                           args=[ast.Constant(value='gravity.quadrature_and_moment_diagnostics')], keywords=[])
        evaluate_fn.body[index] = ast.With(items=[ast.withitem(context_expr=context)], body=[node])
exec(compile(ast.fix_missing_locations(ast.Module(body=[decode_fn, evaluate_fn], type_ignores=[])), str(runner)+'[provider prefix]', 'exec'), globals())
decode = clock.wrap(decode, 'provider.initial_decode')
runtime_sources = {}
for name, module in list(sys.modules.items()):
    path = getattr(module, '__file__', None)
    if path is not None and Path(path).is_relative_to(base) and Path(path).suffix == '.py':
        runtime_sources[name] = dict(path=str(path), sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest())
if args.preflight_only:
    print(json.dumps(dict(status='preflight_only_no_provider_evaluation', cells=n,
        input_time=float(mid['time']), seed_sha256=sha_array(seed), output_created=False,
        resolved_modules={name: str(Path(module.__file__).resolve()) for name, module in
                          [('model', model), ('quad', quad), ('frozen_root', frozen_root),
                           ('vector_root', vector_root), ('clock', clock)]},
        instrumented_ast_matches_original=True,
        runtime_source_sha256=runtime_sources,
        limits='Imports/inputs/function definitions only; no timer or numerical provider result.'), indent=2))
    raise SystemExit(0)
args.output.mkdir(exist_ok=False)
with clock.section('provider.total'):
    dy, diag, z, gravity, moment_defect, node_count, iterations = evaluate(float(mid['time']), y, L, seed, 'midpoint')
hydro = np.column_stack([diag['mass_rates'], dy[:-4].reshape(-1, 20)[:, 3:]])
total = hydro.copy()
total[~frozen] += gravity[~frozen]
old_total = expected['original_hydro'].copy()
old_total[~frozen] += expected['gravity'][~frozen]
old_events = [json.loads(line) for line in (saved_root/'events.jsonl').read_text().splitlines()]
old_iterations = [row for row in old_events if row.get('stage') == 'midpoint' and 'iteration' in row]
comparisons = {}
for name, actual, old in (('weight_rates', z, mid['weight_rates']), ('hydro', hydro, expected['original_hydro']),
                          ('gravity', gravity, expected['gravity']), ('raw_source', total, old_total)):
    tolerance = 1e-14+1e-10*abs(old)
    comparisons[name] = dict(bitwise_equal=actual.tobytes() == old.tobytes(),
        maximum_absolute_difference=float(abs(actual-old).max()),
        maximum_tolerance_fraction=float((abs(actual-old)/tolerance).max()))
assert iterations == len(old_iterations)
assert all(row['maximum_tolerance_fraction'] <= 1 for row in comparisons.values())
report = dict(status='completed', scope=__doc__, flags=flags, cells=n, nodes=node_count,
    time=float(mid['time']), initial_seed='stage0/inputs.npz:weight_rates',
    seed_sha256=sha_array(seed), state_sha256=sha_array(y), factors_sha256=sha_array(L),
    iterations=iterations, original_iterations=len(old_iterations), mass_residual=diag['mass_flux'],
    comparisons=comparisons, timers=clock.report(),
    profile_manifest=profile_manifest,
    runtime_source_sha256=runtime_sources,
    environment=dict(host=platform.node(), python=sys.version, numpy=np.__version__,
                     scipy=sys.modules['scipy'].__version__,
                     OPENBLAS_NUM_THREADS=os.environ.get('OPENBLAS_NUM_THREADS'),
                     OMP_NUM_THREADS=os.environ.get('OMP_NUM_THREADS')),
    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    limits='One instrumented provider replay. Timers include their overhead; exclusive time removes nested measured spans, not uninstrumented callees. Face capture/moment-defect diagnostics included and labeled where measured. Not independent trajectory, speedup measurement, optimization, long-time or physical accuracy proof.')
np.savez(args.output/'replayed.npz', state=y, factors=L, seed=seed, weight_rates=z, hydro=hydro, gravity=gravity, raw_source=total)
(args.output/'report.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2), flush=True)
