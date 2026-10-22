#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
envelope=${1:?immutable native envelope required}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
sha256sum -c "$envelope/source.sha256"
source "$envelope/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
# Independent layouts share the idle Ada GPU; retain separate native logs.
shift
python3 - "$envelope" "$root" "$@" <<'PY'
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json, subprocess, sys
envelope, root = map(Path, sys.argv[1:3])
labels = sys.argv[3:] or ['uniform', 'graded']
if not labels or any(x not in ('uniform', 'graded') for x in labels):
    raise ValueError('expected uniform and/or graded')
def execute(label):
    with (envelope/f'{label}.out').open('x') as out:
        command = [sys.executable, str(envelope/'run_radial_resolution_v1.py'),
            '--root', str(root), '--case', label, '--version', 'v2']
        previous = root/'radial_resolution_20261021_v2/native/uniform/preflight_mpi1'
        if label == 'uniform' and previous.exists():
            command += ['--reuse-preflight', str(previous)]
        result = subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)
    (envelope/f'native_{label}_exit_code.txt').write_text(str(result.returncode)+'\n')
    print(label, 'exit', result.returncode, flush=True)
    return dict(case=label, exit_code=result.returncode)
with ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(execute, labels))
(envelope/('completion_'+'_'.join(labels)+'.json')).write_text(json.dumps(results, indent=2)+'\n')
raise SystemExit(0 if all(x['exit_code'] == 0 for x in results) else 2)
PY
