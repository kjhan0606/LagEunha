"""Summarize actual completed native stages; never infer a missing endpoint."""
from pathlib import Path
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('--root', type=Path, required=True)
parser.add_argument('--version', default='v2')
parser.add_argument('--batch', default='native_v2')
args = parser.parse_args()
base = args.root/('radial_resolution_20261021_'+args.version)
ic = json.loads((base/'manifest.json').read_text())
result = dict(scope='smooth radial IC comparison in the original native CUDA/RK4 solver',
              long_goal_achieved=False, cases={})
for label in ('uniform', 'graded'):
    path = base/args.batch/label/'manifest.json'
    if not path.exists():
        result['cases'][label] = dict(status='NOT_STARTED')
        continue
    data = json.loads(path.read_text())
    init = ic['cases'][label]
    row = dict(n=init['n'], IC_sha256=init['IC_sha256'], initial_mass=init['total_mass'],
               initial_U=init['total_internal_energy'], initial_K=init['total_kinetic_energy'],
               initial_L=init['total_angular_momentum'], radial_mesh=init['radial_mesh'],
               status=data.get('status', 'INCOMPLETE; inspect live handle'),
               successful_evolution_checkpoints=[], failures=[])
    if 'preflight' in data:
        pre = data['preflight']
        rhs = pre['first_RHS']
        row.update(native_volume_relative_max=pre['native_volume_relative_max'],
                   first_RHS=rhs,
                   positive_Udot_cells_initial_R_lt_half=sum(r['positive_Udot'] for r in rhs if r['Rlo'] < .5),
                   signed_Udot_cells_initial_R_lt_half=sum(r['signed_Udot'] for r in rhs if r['Rlo'] < .5),
                   face_geometry=pre['explicit_tiny_face_audit'])
    if 'mpi1_mpi4' in data:
        row['mpi1_mpi4'] = data['mpi1_mpi4']
    for stage in data['stages']:
        if stage.get('exit_code'):
            row['failures'].append(stage)
        elif 'physical_audit' in stage:
            row['successful_evolution_checkpoints'].append(stage['physical_audit'])
    result['cases'][label] = row
a, b = (result['cases'][label] for label in ('uniform', 'graded'))
if all('positive_Udot_cells_initial_R_lt_half' in v for v in (a, b)):
    numerator, denominator = a['positive_Udot_cells_initial_R_lt_half'], b['positive_Udot_cells_initial_R_lt_half']
    result['uniform_over_graded_first_positive_Udot_ratio'] = numerator/denominator if denominator else None
result['comparison_caveats'] = [
    'Positive initial Udot is not pure irreversible heating and not a long-run stability proof.',
    'Only checkpoints explicitly listed have completed physical times.',
    'Failure snapshots contain invalid RK stages; their normal starts are the valid restart states.',
    'Both paths use the same frozen v4 binary; intrinsic-spin and new remap operators are not present.',
    'The two paths share the Ada GPU; elapsed times do not establish exclusive GPU speed or scaling.'
]
output = base/('summary_'+args.batch+'.json')
output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k != 'cases'}, indent=2))
for name, row in result['cases'].items():
    print(name, row.get('n'), row['status'],
          'completed_times', [x['time'] for x in row.get('successful_evolution_checkpoints', [])])
