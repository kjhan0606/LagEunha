"""Independent saved-output checks only. No physical query or evolution."""
import csv
import hashlib
import json
import pathlib
from decimal import Decimal, getcontext

getcontext().prec = 60
root = pathlib.Path(__file__).resolve().parent
run = root / 'material_full_GPU_mass_v1_h200_mpi4'
dest = root / 'material_full_GPU_mass_audit_v1'
dest.mkdir()
result = json.loads((run / 'run/full_GPU_mass_continuation_result.json').read_text())
assert (run / 'native_exit_code.txt').read_text().strip() == '0'
assert result['first_GPU_mass_continuation_pass'] and not result['replayed_accepted_trajectory']
assert result['same_state_v2_restart_roundtrip']
checked_files = 0
for name in ('source.sha256', 'binary.sha256', 'actual_retained_restart.sha256'):
    for line in (run / name).read_text().splitlines():
        digest, filename = line.split(maxsplit=1)
        assert hashlib.file_digest(pathlib.Path(filename).open('rb'), 'sha256').hexdigest() == digest
        checked_files += 1
step = json.loads((run / 'run/first_GPU_mass_step/step.json').read_text(), parse_float=Decimal)
assert step['finite_step_gates_pass']
for field, bound in [('energy_defect_over_U', '1e-8'), ('thermal_increment_defect_over_U', '1e-8'), ('angular_defect_over_2K', '1e-11')]:
    assert abs(step[field]) < Decimal(bound)
summaries = []
for stage in range(1, 5):
    d = run / f'run/first_GPU_mass_step/stage{stage}'
    rhs = json.loads((d / 'rhs.json').read_text(), parse_float=Decimal)
    gpu = json.loads((d / 'GPU_mass.json').read_text(), parse_float=Decimal)
    assert rhs['RHS_gates_pass'] and rhs['GPU_mass_acceleration_implemented']
    assert max(gpu['true_relative_residual']) < Decimal('1e-16')
    assert max(gpu['outer_corrections']) <= 3
    assert gpu['binary80_to_hi_lo_restore_error'] == 0
    seen = set()
    energy = angular = Decimal(0)
    for rank in range(4):
        with (d / f'cells_r{rank}.csv').open() as f:
            for row in csv.DictReader(f):
                cell = int(row['cell'])
                assert cell not in seen
                seen.add(cell)
                assert Decimal(row['volume']) > 0 and Decimal(row['density']) > 0 and Decimal(row['U']) > 0
                er = Decimal(row['energy_rate_error_over_U'])
                lr = Decimal(row['angular_rate_error_over_2K'])
                energy, angular = max(energy, er), max(angular, lr)
                spin, inertia, rot = map(Decimal, (row['spin'], row['inertia'], row['Erot']))
                assert inertia > 0 and abs(rot - spin*spin/(2*inertia)) <= Decimal('1e-18') * max(rot, Decimal('1e-40'))
    assert len(seen) == 71080 and energy < Decimal('1e-8') and angular < Decimal('1e-11')
    assert abs(energy-rhs['every_cell_energy_rate_error_over_U']) < Decimal('1e-30')
    assert abs(angular-rhs['every_cell_angular_rate_error_over_2K']) < Decimal('1e-30')
    summaries.append({'stage': stage, 'cells': len(seen), 'maximum_E_rate_error_over_U': str(energy),
                      'maximum_L_rate_error_over_2K': str(angular),
                      'mass_solve_seconds': str(rhs['mass_solve_seconds_rank0'])})
output = {'scope': 'saved native full GPU mass continuation, all four stages and 71080 cells; no replay',
          'accepted_time': str(step['time_after']), 'native_hashes_checked': checked_files,
          'source_mass_fixed': True, 'no_independent_spin_reservoir': True,
          'stage_audits': summaries, 'audit_pass': True, 'long_Kepler_validated': False}
(dest / 'contract.json').write_text(json.dumps(output, indent=2)+'\n')
print(json.dumps(output, indent=2))
