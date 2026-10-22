"""Saved-output audit and radial heating summaries; no field evolution."""
import csv
import hashlib
import json
import pathlib
import struct
import sys
from decimal import Decimal, getcontext

getcontext().prec = 60
root = pathlib.Path(__file__).resolve().parent
run = pathlib.Path(sys.argv[1]).resolve()
dest = pathlib.Path(sys.argv[2]).resolve()
assert run.parent == root and dest.parent == root
assert (run / 'native_exit_code.txt').read_text().strip() == '0'
result = json.loads((run / 'run/full_continuous_result.json').read_text(), parse_float=Decimal)
assert result['new_accepted_steps'] > 0 and result['same_accepted_trajectory']
assert result['no_replayed_accepted_interval'] and result['GPU_mass_DD'] and result['GPU_point_cache']
assert result['accepted_time'] > result['time_before']
dest.mkdir()
checked = 0
for name in ('source.sha256', 'binary.sha256', 'accepted_input.sha256'):
    for line in (run / name).read_text().splitlines():
        digest, filename = line.split(maxsplit=1)
        assert hashlib.file_digest(pathlib.Path(filename).open('rb'), 'sha256').hexdigest() == digest
        checked += 1
prefix = result['last_checkpoint_prefix']
step_dir = run / 'run' / pathlib.Path(prefix).parent
step = json.loads((step_dir / 'step.json').read_text(), parse_float=Decimal)
assert step['finite_step_gates_pass'] and step['time_after'] == result['accepted_time']
for field, bound in [('energy_defect_over_U', '1e-8'), ('thermal_increment_defect_over_U', '1e-8'), ('angular_defect_over_2K', '1e-11')]:
    assert abs(step[field]) < Decimal(bound)
for rank in range(4):
    assert (step_dir / f'accepted_r{rank}.bin').stat().st_size > 86000000
for stage in range(1, 5):
    d = step_dir / f'stage{stage}'
    rhs = json.loads((d / 'rhs.json').read_text(), parse_float=Decimal)
    gpu = json.loads((d / 'GPU_mass.json').read_text(), parse_float=Decimal)
    assert rhs['RHS_gates_pass'] and max(gpu['true_relative_residual']) < Decimal('1e-16')
    assert max(gpu['outer_corrections']) <= 3
    seen = set()
    for rank in range(4):
        with (d / f'cells_r{rank}.csv').open() as f:
            for row in csv.DictReader(f):
                cell = int(row['cell'])
                assert cell not in seen
                seen.add(cell)
                assert Decimal(row['U']) > 0 and Decimal(row['volume']) > 0
                assert Decimal(row['energy_rate_error_over_U']) < Decimal('1e-8')
                assert Decimal(row['angular_rate_error_over_2K']) < Decimal('1e-11')
    assert len(seen) == 71080
cache = json.loads((step_dir / 'stage1/cache_equivalence.json').read_text(), parse_float=Decimal)
assert cache['pass'] and cache['maximum_relative_reduced_difference'] < Decimal('1e-11')
initial_U = {}
for rank in range(4):
    with (root / 'material_full_initial_ledger_v1_mpi4/run' / f'full_initial_local_ledger_r{rank}.csv').open() as f:
        for row in csv.DictReader(f):
            initial_U[int(row['cell'])] = Decimal(row['U'])
# Original gamma is stored as the native float converted to double.
gamma = Decimal(struct.unpack('f', struct.pack('f', 1.666667))[0])
P0 = Decimal('0.000001')
edges = list(map(Decimal, ['0', '.5', '1', '2', '4', '8', '12', '100']))
bins = [{'count': 0, 'U': Decimal(0), 'U_initial_same_cells': Decimal(0),
         'max_abs_cell_U_ratio_minus1': Decimal(0), 'max_abs_area_P_ratio_minus1': Decimal(0)} for _ in edges[:-1]]
seen = set()
for rank in range(4):
    with (step_dir / 'endpoint' / f'physical_cells_r{rank}.csv').open() as f:
        for row in csv.DictReader(f):
            cell = int(row['cell'])
            assert cell not in seen and cell in initial_U
            seen.add(cell)
            mass, volume, density, P, U = [Decimal(row[k]) for k in ('mass', 'volume', 'density', 'area_mean_pressure', 'U')]
            assert min(mass, volume, density, P, U) > 0
            assert abs(density-mass/volume) < Decimal('1e-18')*density
            assert abs(P-(gamma-1)*U/volume) < Decimal('1e-18')*P
            x, y = Decimal(row['x_mass']), Decimal(row['y_mass'])
            radius = (x*x+y*y).sqrt()
            i = next(i for i in range(len(edges)-1) if edges[i] <= radius < edges[i+1])
            b = bins[i]; b['count'] += 1; b['U'] += U; b['U_initial_same_cells'] += initial_U[cell]
            b['max_abs_cell_U_ratio_minus1'] = max(b['max_abs_cell_U_ratio_minus1'], abs(U/initial_U[cell]-1))
            b['max_abs_area_P_ratio_minus1'] = max(b['max_abs_area_P_ratio_minus1'], abs(P/P0-1))
assert len(seen) == 71080
out = {'scope': 'actual next full accepted RK4 stage/endpoint/hash audit; saved data only',
       'accepted_time': str(result['accepted_time']), 'hashes_checked': checked,
       'cells_audited_per_stage_and_endpoint': len(seen), 'cache_comparison': str(cache['maximum_relative_reduced_difference']),
       'radial_bin_edges_current_mass_centroid': list(map(str, edges)),
       'radial_heating': [{k: str(v) if isinstance(v, Decimal) else v for k,v in b.items()} for b in bins],
       'global_map_injectivity_proven': False, 'long_Kepler_validated': False, 'audit_pass': True}
(dest / 'contract.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
