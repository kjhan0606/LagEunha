"""Compare partitioned MPI gravity accumulation with retained serial backends."""
import hashlib
import json
from pathlib import Path
import numpy as np

base = Path('/gpfs/kjhan/LagForce/Codex_gravity_quadrature_backends_20261018')
root = Path('/gpfs/kjhan/LagForce/Codex_gravity_quadrature_mpi_20261018')
rows = []
for ranks in (2, 4):
    for backend in ('cpu1', 'omp4', 'cuda'):
        name = f'mpi{ranks}_{backend}.bin'
        actual = np.fromfile(root / name, dtype='<f8').reshape(1592, 11)
        reference = np.fromfile(base / f'{backend}.bin', dtype='<f8').reshape(1592, 11)
        assert np.isfinite(actual).all()
        assert np.all(actual[:1590, 0] == 0)
        assert np.all(actual[1590:, 0] == 1)
        assert np.all(actual[1590:, 1:] == 12345.)
        bitwise = actual.tobytes() == reference.tobytes()
        assert bitwise, name
        rows.append(dict(file=name, ranks=ranks, backend=backend,
                         bitwise_equal_to_serial_backend=bitwise,
                         maximum_absolute_difference=float(np.max(abs(actual-reference)))))

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

report = dict(valid_cells=1590, invalid_cases=2, results=rows,
              sha256={name: sha(root/name) for name in
                      ('moment_gravity_quadrature.h', 'check_gravity_quadrature_backends.cu', 'probe_mpi')},
              reference_report_sha256=sha(base/'backend_report.json'),
              limits='Independent cell partition and MPI result collection on one node; GPU ranks share one Ada. No halo exchange, full hydro, multi-GPU or long-term stability validation.')
with (root/'mpi_report.json').open('x') as f:
    json.dump(report, f, indent=2)
print(json.dumps(report, indent=2))
