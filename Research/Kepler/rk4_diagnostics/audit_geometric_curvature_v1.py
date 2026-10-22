#!/usr/bin/env python3
"""Merge saved engineering checks; never replace finite by zero curvature."""
import json
from pathlib import Path

r = Path(__file__).resolve().parent
base = json.loads((r / 'GEOMETRIC_PRESSURE_LEDGER_20261021.json').read_text())
curvature = {}
for n in (1, 4):
    p = r / f'geometric_curvature_v1_mpi{n}'
    assert (p / 'native_exit_code.txt').read_text().strip() == '0'
    for name in ('final_source_integrity.txt', 'final_input_integrity.txt'):
        assert all(x.endswith(': OK') for x in (p / name).read_text().splitlines())
    states = {s: json.loads((p / f'{s}_result.json').read_text()) for s in ('initial', 'finite')}
    assert all(s['actual_nu_temporal_gate_pass'] for s in states.values())
    a = states['initial']
    assert a['actual_nu_pressure_curvature_chi_max'] <= .002
    assert a['actual_nu_every_cell_mean_chi_max'] <= .001
    for s in states.values():
        assert s['nu_vs_material_curvature_correction_chi_max'] <= .002
        assert s['every_cell_mean_correction_chi_max'] <= .001
    curvature[str(n)] = states

comparison = {}
for state in ('initial', 'finite'):
    fields = {}
    for key in ('actual_nu_pressure_curvature_chi_max',
                'actual_nu_every_cell_mean_chi_max',
                'nu_vs_material_curvature_correction_chi_max'):
        a, b = curvature['1'][state][key], curvature['4'][state][key]
        fields[key] = {'absolute_difference_in_chi_units': abs(a - b),
                       'relative_difference': abs(a - b) / max(abs(a), abs(b), 1e-300)}
    comparison[state] = fields

# A preliminary inline postprocessing assertion used a FORCE relative
# tolerance for a cancellation difference near1e-14. That diagnostic
# assertion failed, while every actual chi/positivity/physical gate passed.
# Preserve the discrepancy explicitly; never change a native physical gate
# or rerun a simulation to hide it. Relative comparison of a near-zero
# correction is not a new physical accuracy criterion.
ledger = {
    'scope': 'compatible pressure/heat saved-state bundle plus actual nu curvature, native CPU/OpenMP/MPI/CUDA',
    'pressure_audit_ledger': 'GEOMETRIC_PRESSURE_LEDGER_20261021.json',
    'pressure_same_state_checks_pass': base['all_same_state_native_checks_pass'],
    'actual_nu_curvature': curvature,
    'MPI_curvature_reduction_comparison': comparison,
    'preliminary_inline_audit_assertion_failed': True,
    'inline_audit_failure_reason': 'force relative1e-11 was incorrectly applied to the near-zero nu-minus-J curvature correction; native physics gates unchanged',
    'initial_physical_curvature_gates_pass': True,
    'finite_control_is_an_accepted_time_step': False,
    'finite_ideal_Kepler_curvature_pass_claimed': False,
    'thermal_restart_for_general_evolved_state_implemented': False,
    'whole_target_68526_transfer_certified': False,
    'production_long_Kepler_validated': False,
    'all_native_exec_sessions_closed': [90245, 42750, 16281, 33998, 60880, 20383],
}
with (r / 'GEOMETRIC_PRESSURE_COMPLETE_LEDGER_20261021.json').open('x') as out:
    json.dump(ledger, out, indent=2)
    out.write('\n')
print(json.dumps({'pressure_and_actual_nu_physical_gates': 'PASS', 'MPI_reduction_comparison': comparison}))
