"""Verify native thermal64 checkpoint/restart against continuous execution."""
import hashlib
import json
from pathlib import Path
import sys
from audit_kepler_restart_bytes_20261029 import compare

root=Path(sys.argv[1]);continuous=root/'bridge_untruncated';restart=root/'bridge_restart';reference=root/'omp'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def accepted(d):
    return [s[s.index('[WALL_ACCEPT]'):] for s in (d/'log').read_text(errors='replace').splitlines() if '[WALL_ACCEPT]' in s]
paths=[continuous/'khout.000026.dat',restart/'khout.000026.dat']
comparison=compare(*paths);checkpoint_output=compare(reference/'khout.000026.dat',paths[0])
ca=accepted(continuous);ra=accepted(restart)
checks={'full_restart_record_and_footer_parity':comparison['passed'],
        'intermediate_output_does_not_change_final':checkpoint_output['passed'],
        'continuous_has_four_accepted_steps':len(ca)==4,
        'restart_has_two_accepted_steps':len(ra)==2,
        'resumed_accepted_ledgers_exact':ca[2:]==ra,
        'both_exit_zero':all((d/'exit_status.txt').read_text().strip()=='0' for d in (continuous,restart)),
        'both_enable_omp_and_thermal64':all('GFS_THERMAL64=1 GFS_WALL_OMP_VALIDATE=1' in (d/'run.sh').read_text() for d in (continuous,restart)),
        'restart_reads_exact_intermediate':sha(restart/'khout.000024.dat')==sha(continuous/'khout.000024.dat'),
        'unchanged_binary':sha(root/'project/eunha2.exe')=='da70844984a52fec77e2d6dec451e9b5695e761ba93a42d0f0d079c1091edf44'}
result={'passed':all(checks.values()),'checks':checks,'restart_comparison':comparison,
        'intermediate_output_comparison':checkpoint_output,
        'hashes':{str(p):sha(p) for p in (continuous/'khout.000024.dat',continuous/'log',restart/'log',continuous/'run.sh',restart/'run.sh',Path(__file__))},
        'scope':'Executed native thermal64 Q1 short no-impact restart from an untruncated intermediate step, OpenMP guard enabled in both paths. Shared GPU, one rank; no physical wall crossing, performance or long-time stability proof.'}
with (root/'bridge_verified.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit(1)
