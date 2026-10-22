#!/usr/bin/env python3
"""Saved-output audit only. Does not evaluate physical fields or evolve state."""
import csv
import hashlib
import json
import struct
import sys
from decimal import Decimal as D
from pathlib import Path

run, root, dest = map(Path, sys.argv[1:])
dest.mkdir()
def read_json(p):
    return json.loads(p.read_text(), parse_float=D)
def verify_manifest(p):
    count = 0
    for line in p.read_text().splitlines():
        expected, name = line.split('  ', 1)
        h = hashlib.sha256()
        with Path(name).open('rb') as f:
            for chunk in iter(lambda: f.read(1048576), b''):
                h.update(chunk)
        assert h.hexdigest() == expected, name
        count += 1
    return count
manifests = {p: verify_manifest(run / p) for p in
             ['source.sha256', 'retained_input.sha256', 'binary.sha256', 'output.sha256']}
contract = read_json(run / 'contract.json')
assert contract['bundle_pass'] and contract['actual_steps_executed'] == 3
assert contract['actual_RHS_executed'] == 12 and contract['MPI_ranks'] == 4
assert contract['restart_same_backend_bitmatch'] and not contract['long_Kepler_validated']
saved = (root / 'shared_trace_rhs_v1_mpi4/shared_state.bin').read_bytes()
h = struct.unpack_from('<8Q', saved)
assert h[3:7] == (3797, 384, 1064, 32)
offset = 64 + h[4] * h[5]
old_coeff = [None] * h[3]
for i in range(h[3]):
    ident = struct.unpack_from('<4Q', saved, offset)
    assert ident[0] < h[3] and ident[3] == ident[0] % 4
    old_coeff[ident[0]] = saved[offset + 32:offset + 64]
    offset += 64
assert offset == len(saved) and all(old_coeff)
old_coeff = b''.join(old_coeff)
names = ['accepted_first', 'uninterrupted_second', 'restarted_second']
metrics = {}
bound = {
    'global_energy_rate_error_over_U': D('1e-8'),
    'every_cell_energy_rate_error_over_U': D('1e-8'),
    'global_angular_rate_error_over_2K': D('1e-11'),
    'every_cell_angular_rate_error_over_2K': D('1e-11'),
    'CPU_CUDA_force_relative_difference': D('1e-11'),
    'hybrid_global_reaction_residual_relative': D('1e-11'),
    'zero_trace_bubble_reaction_relative': D('1e-11'),
    'external_traction_physical_relative_difference': D('1e-8'),
    'corrected_volume_identity_relative_error_max': D('1e-10'),
    'stage_EOS_consistency_relative_error_max': D('1e-10'),
}
stages = []
for name in names:
    for stage in range(1, 5):
        p = run / name / f'stage{stage}'
        r = read_json(p / 'result.json')
        assert r['evolved_stage_algebra_pass']
        assert all(r[k] > 0 for k in ['min_J', 'min_nu', 'min_physical_U'])
        assert r['arbitrary_constant_face_force'] == 0
        for k, limit in bound.items():
            assert abs(r[k]) < limit, (name, stage, k)
            metrics[k] = max(metrics.get(k, D(0)), abs(r[k]))
        witness = [(p / f'coeff_rank{rank}.bin').read_bytes() for rank in range(4)]
        assert len(witness[0]) == 3797 * 32 and all(x == witness[0] for x in witness)
        if name == 'accepted_first' and stage == 1:
            assert witness[0] == old_coeff  # actual saved state, not a fresh IC
        else:
            assert witness[0] != old_coeff  # genuinely updated x/v
        if name == 'restarted_second':
            assert witness[0] == (run / 'uninterrupted_second' / f'stage{stage}/coeff_rank0.bin').read_bytes()
        stages.append({'step': name, 'stage': stage,
                       'shared_coefficient_sha256': hashlib.sha256(witness[0]).hexdigest()})
steps = [read_json(run / n / 'step.json') for n in names]
for s in steps:
    assert s['finite_step_gate_pass'] and s['RK4_stages'] == 4
    assert s['mass_before'] == s['mass_after']
    assert abs(s['energy_defect_over_U']) < D('1e-8')
    assert abs(s['angular_defect_over_2K']) < D('1e-11')
    assert abs(s['thermal_increment_defect_over_U']) < D('1e-8')
assert steps[1] == steps[2]
# Compare unchanged physical quantities at the first actual stage to the
# already saved source-cell diagnostics, without another field integration.
def rows(p):
    with p.open() as f:
        return {r['cell']: r for r in csv.DictReader(f)}
old = rows(root / 'geometric_pressure_v2_mpi4/finite/cell_physics_rotation.csv')
new = rows(run / 'accepted_first/stage1/cell_physics_rotation.csv')
assert old.keys() == new.keys() and len(old) == 24
for cell in old:
    for key in ['U', 'K', 'internal_L', 'inertia', 'omega', 'rotation_energy']:
        assert D(old[cell][key]) == D(new[cell][key]), (cell, key)
report = {
    'scope': 'independent audit of actual saved native CUDA/RK4/restart outputs; no physical queries or simulation replay',
    'manifest_files_verified': manifests,
    'actual_stages_audited': len(stages),
    'all4rank_shared_coefficients_bitmatch_every_stage': True,
    'first_stage_matches_original_finite_CONTROL_bitwise': True,
    'later_stages_are_actual_changed_coefficients': True,
    'all4restart_stage_witnesses_bitmatch_uninterrupted': True,
    'restart_finite_step_ledgers_exactly_match': True,
    'all24_initial_cell_U_K_internal_L_I_omega_Krot_match_saved': True,
    'stage_metric_maxima': {k: str(v) for k, v in metrics.items()},
    'stage_witnesses': stages,
    'audit_pass': True,
    'full_disk_or_long_Kepler_validated': False,
}
(dest / 'contract.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'stage_witnesses'}, indent=2))
