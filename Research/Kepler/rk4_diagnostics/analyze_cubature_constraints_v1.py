"""Analyze saved native matrices; no physics query, rule change or evolution."""
import hashlib
import json
import os
import pathlib
import time

os.environ['OPENBLAS_NUM_THREADS'] = '8'
os.environ['MKL_NUM_THREADS'] = '8'
os.environ['OMP_NUM_THREADS'] = '8'
import numpy as np

root = pathlib.Path(__file__).resolve().parent
source = root / 'material_cubature_feasibility_v1_mpi4'
dest = root / 'material_cubature_matrix_analysis_v1'
dest.mkdir()
adoption = json.loads((source / 'export_adoption.json').read_text())
assert adoption['data_adopted_for_saved_matrix_analysis_only']
reports = []
for rank in range(4):
    entries = [json.loads(line) for line in (source / 'run' / f'constraints_rank{rank}.jsonl').read_text().splitlines()]
    for entry in entries:
        start = time.monotonic()
        p = source / 'run' / entry['file']
        with p.open('rb') as f:
            h = np.fromfile(f, dtype='<u8', count=8)
            assert int(h[0]) == 0x4d46435542434e31 and int(h[1]) == 1
            rows, n, size = map(int, h[3:6])
            assert size == np.dtype(np.longdouble).itemsize == 16
            index = np.fromfile(f, dtype='<u4', count=n)
            w = np.fromfile(f, dtype='<f8', count=n)
            A = np.fromfile(f, dtype=np.longdouble, count=rows*n).reshape(rows, n)
            target = np.fromfile(f, dtype=np.longdouble, count=rows)
            absolute = np.fromfile(f, dtype=np.longdouble, count=rows)
            assert f.read() == b''
        assert np.all(w > 0) and np.isfinite(A).all()
        recomputed = np.sum(A*w.astype(np.longdouble), axis=1, dtype=np.longdouble)
        error = np.abs(recomputed-target)
        assert np.all(error <= 64*np.finfo(np.longdouble).eps*absolute)
        # Floating-point singular values only diagnose conditioning. No row or
        # small direction is deleted, and no reduced rule is accepted here.
        weighted = A * np.sqrt(w.astype(np.longdouble))
        scale = np.sqrt(np.sum(weighted*weighted, axis=1, dtype=np.longdouble))
        nonzero = scale > 0
        B = np.asarray(weighted[nonzero]/scale[nonzero, None], dtype=np.float64)
        singular = np.linalg.svd(B, compute_uv=False)
        floor = 64*max(B.shape)*np.finfo(np.float64).eps*singular[0]
        resolved = int(np.count_nonzero(singular > floor))
        report = {'rank': rank, 'triangle': entry['triangle'], 'cell': entry['cell'],
                  'constraints': rows, 'positive_points': n,
                  'float64_resolved_singular_directions_diagnostic_only': resolved,
                  'diagnostic_roundoff_floor': float(floor),
                  'singular_values': list(map(float, singular)),
                  'mathematical_rank_certified': False,
                  'numerically_small_directions_discarded': False,
                  'positive_reduced_rule_constructed': False,
                  'candidate_rule_committed': False,
                  'maximum_native_target_recompute_error_over_abs_sum': float(np.max(error[nonzero]/absolute[nonzero])),
                  'source_sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest(),
                  'matrix_analysis_seconds': time.monotonic()-start}
        reports.append(report)
        (dest / f'triangle_{entry["triangle"]}_rank{rank}.json').write_text(json.dumps(report, indent=2)+'\n')
output = {'scope': 'fixed 24-source saved native constraint matrix conditioning; not a reduced physical step',
          'representatives': len(reports), 'positive_points': sum(x['positive_points'] for x in reports),
          'constraint_count_range': [min(x['constraints'] for x in reports), max(x['constraints'] for x in reports)],
          'diagnostic_resolved_direction_range': [min(x['float64_resolved_singular_directions_diagnostic_only'] for x in reports), max(x['float64_resolved_singular_directions_diagnostic_only'] for x in reports)],
          'maximum_target_recompute_error_over_abs_sum': max(x['maximum_native_target_recompute_error_over_abs_sum'] for x in reports),
          'native_MPI_shutdown_error_preserved': True,
          'positive_reduced_rule_constructed': False, 'candidate_rule_committed': False,
          'mass_operator_changed': False, 'new_physical_evolution': False}
(dest / 'contract.json').write_text(json.dumps(output, indent=2)+'\n')
print(json.dumps(output, indent=2))
