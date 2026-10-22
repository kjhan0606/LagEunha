"""Saved native evidence, no physics run, state change or tolerance change."""
import csv
import json
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def doc(folder, name='contract.json'):
    return json.loads((ROOT / folder / name).read_text())


def main():
    factories = {str(r): doc(f'shared_trace_factory_v1_mpi{r}') for r in (1, 4)}
    for ranks, entry in factories.items():
        assert (ROOT / f'shared_trace_factory_v1_mpi{ranks}' / 'native_exit_code.txt').read_text().strip() == '0'
        assert entry['factory_initial_gate_pass']
        assert entry['new_scalar_dofs'] == 1853 + 56 + 2 * 560 + 2 * 384 == 3797
        assert entry['positive_points'] == 817152 and entry['invalid_points'] == 0
        assert entry['canonical_trace_jump'] == entry['bubble_trace'] == 0
        assert entry['checkpoint_and_transaction_pass']
    files = [ROOT / f'shared_trace_factory_v1_mpi{r}/shared_coefficients.bin' for r in (1, 4)]
    data = [f.read_bytes() for f in files]
    arrays = [struct.unpack('<%dd' % (len(b) // 8), b) for b in data]
    assert arrays[0] == arrays[1]
    signed_zero = 0
    for j, (a, b) in enumerate(zip(*arrays)):
        if data[0][8*j:8*j+8] != data[1][8*j:8*j+8]:
            assert a == b == 0
            signed_zero += 1
    assert signed_zero == 49
    rhs = doc('shared_trace_rhs_v1_mpi4')
    assert (ROOT / 'shared_trace_rhs_v1_mpi4/native_exit_code.txt').read_text().strip() == '2'
    assert rhs['global_physical_algebra_gate_pass']
    assert rhs['initial_chi_max'] > .002 and rhs['initial_every_cell_mean_chi_max'] <= .001
    assert not rhs['old_curvature_gate_pass'] and not rhs['full_factory_bundle_pass']
    witness = doc('shared_trace_witness_v3')
    mass = doc('shared_trace_witness_v3', 'mass_metric_audit.json')
    pressure = doc('shared_trace_witness_v3', 'pressure_load_audit.json')
    assert witness['old_curvature_failure_reproduced'] and witness['worst_source'] == 5
    assert mass['represented_acceleration_action_relative_difference'] < 2e-17
    assert pressure['pressure_fraction_at_saved_worst_point'] > .9999999
    assert pressure['gravity_remainder_chi_max'] < .002
    with (ROOT / 'shared_trace_witness_v3/curvature_witnesses.csv').open() as f:
        worst = next(csv.DictReader(f))
    assert int(worst['source']) == 5 and int(worst['cell']) == 17894
    result = {
        'scope': 'same24-cell prospective shared state and native sparse RHS; decisive initial curvature failure',
        'factory_initial_native': factories,
        'MPI1_MPI4_coefficient_values_identical': True,
        'MPI1_MPI4_coefficient_bytes_identical': False,
        'signed_zero_byte_differences': signed_zero,
        'new_RHS_finite_control': rhs,
        'saved_state_worst_witness': worst,
        'same_basis_mass_audit': mass,
        'same_load_linear_pressure_decomposition': pressure,
        'old_degree_mesh_order_precision_guards_or_tolerances_changed': False,
        'old_broken_packet_projected_or_replaced': False,
        'new_shared_state_promoted_to_transfer_or_physical_evolution': False,
        'local_hybrid_traction_and_rotation_rate_complete': False,
        'analytic_cache_error_certificate_and_new_target_compiler_complete': False,
        'physical_t6_8_10_and_longest_GIZMO_stability_achieved': False,
        'goal_remains_active': True,
    }
    with (ROOT / 'SHARED_TRACE_FACTORY_LEDGER_20261021.json').open('x') as f:
        json.dump(result, f, indent=2)
    print('Initial factory/native algebra pass; unchanged max-curvature gate FAIL.')
    print('Saved pressure fraction at worst point:', pressure['pressure_fraction_at_saved_worst_point'])
    print('Gravity remainder max chi:', pressure['gravity_remainder_chi_max'])


if __name__ == '__main__':
    main()
