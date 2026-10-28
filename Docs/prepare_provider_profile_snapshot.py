"""Prepare instrumented immutable module copies; never evaluates a provider."""
import argparse
import ast
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
base = Path('/gpfs/kjhan/LagForce')
manifest_path = base/'Codex_compiled_midpoint_once_20261018/manifest.json'
manifest = json.loads(manifest_path.read_text())
args.output.mkdir(exist_ok=False)
provenance = {}
for name, meta in manifest['modules'].items():
    path = Path(meta['path'])
    source = path.read_bytes()
    assert hashlib.sha256(source).hexdigest() == meta['sha256']
    (args.output/path.name).write_bytes(source)
    provenance[name] = dict(original_path=str(path), original_sha256=meta['sha256'])
# frozen/vector are dependencies not included in the original runner manifest.
original_dir = base/'Codex_compiled_midpoint_modules_20261018'
search_dirs = [original_dir, base/'Codex_quadratic_thermal_audit_modules_20261018',
               base/'Codex_quadratic_split_modules_20261018',
               base/'Codex_coupled_relaxation_evolution_20261018', base]
for name in ('frozen_weight_root', 'vector_kinetic_root', 'euler_star_pressure', 'vector_moment_flux'):
    path = next(folder/(name+'.py') for folder in search_dirs if (folder/(name+'.py')).is_file())
    source = path.read_bytes()
    (args.output/path.name).write_bytes(source)
    provenance[name] = dict(original_path=str(path), original_sha256=hashlib.sha256(source).hexdigest(),
                            note='dependency source SHA captured now; not in original runner manifest')


def instrument(name, function, select):
    path = args.output/(name+'.py')
    source = path.read_text()
    lines = source.splitlines(keepends=True)
    fn = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == function)
    spans = []
    for node in fn.body:
        label = select(node, ast.get_source_segment(source, node))
        if label:
            spans.append((node.lineno-1, node.end_lineno, node.col_offset, label))
    for lo, hi, indent, label in sorted(spans, reverse=True):
        lines[lo:hi] = [' '*indent+f'with _profile_section({label!r}):\n']+['    '+line for line in lines[lo:hi]]
    result = 'from provider_profile_clock import section as _profile_section\n'+''.join(lines)
    ast.parse(result)
    path.write_text(result)
    return [label for _, _, _, label in spans]


def rhs_selection(node, text):
    if text.startswith('for edge_id,'): return 'rhs.face_record_construction'
    if text.startswith("if backend_name()=='batched'"): return 'rhs.moment_flux_accumulation'
    if text.startswith('for i,s in enumerate(states):'): return 'rhs.compatible_volume_force'


changes = dict(kinetic_factor_rhs=instrument('kinetic_factor_rhs', 'rhs', rhs_selection),
               vector_kinetic_root=instrument('vector_kinetic_root', 'evaluate',
                   lambda node, text: 'massroot.jacobian_build' if text.startswith('if with_jacobian:') else None))
assert len(changes['kinetic_factor_rhs']) == 3 and len(changes['vector_kinetic_root']) == 1
clock = Path(__file__).parent/'polybnd/provider_profile_clock.py'
(args.output/clock.name).write_bytes(clock.read_bytes())
for name in provenance:
    provenance[name]['instrumented_sha256'] = hashlib.sha256((args.output/(name+'.py')).read_bytes()).hexdigest()
report = dict(status='prepared_not_executed', scope=__doc__, original_manifest=str(manifest_path),
              modules=provenance, instrumented_blocks=changes,
              timer_sha256=hashlib.sha256(clock.read_bytes()).hexdigest(),
              prepare_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              limits='Source snapshot preparation only. No physics/model/order/cache/stopping changes. Timer overhead included in execution measurements.')
(args.output/'profile_manifest.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
