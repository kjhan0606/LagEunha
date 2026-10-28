import json,re,subprocess
from pathlib import Path
root=Path('/gpfs/kjhan/LagForce/Codex_mpi_multistep_20261018')
source=Path('/gpfs/kjhan/LagForce/Codex_mpi_coupled_step_20261017')
records={};snapshots={}
for case in ['serial','np4']:
    run=root/case
    assert (run/'exit.code').read_text().strip()=='0'
    log=(run/'log').read_text()
    states=[dict(re.findall(r'(\w+)=([^\s]+)',s)) for s in log.splitlines() if s.startswith('[ACCEPTED_STATE]')]
    assert len(states)==8,(case,len(states))
    assert [int(s['step']) for s in states]==list(range(1,9))
    assert all(s['finite']=='1' and s['retries']=='0' for s in states)
    assert abs(float(states[-1]['elapsed'])-.004)<1e-14
    budgets=[dict(re.findall(r'(\w+)=([^\s]+)',s)) for s in log.splitlines() if s.startswith('[ACCEPTED_RK_BUDGET]')]
    assert len(budgets)==8 and all(s['stage_mask_min']==s['stage_mask_max']=='15' for s in budgets)
    phases=re.findall(r'\[MPI_COUPLED_STAGE\] step=(\d+) retry=(\d+) phase=(\d+) rk=(\d+) epoch=(\d+)',log)
    if case=='np4':
        assert len(phases)==48
        assert [int(s[4]) for s in phases]==list(range(1,49))
    paths=[p for p in run.iterdir() if p.name.startswith('khout.') and not p.is_symlink()]
    assert len(paths)==1,(case,paths)
    snapshots[case]=paths[0]
    records[case]={'accepted':states,'phase_count':len(phases)}
report=root/'particle_comparison.json'
subprocess.run(['python3',str(source/'compare_particles.py'),'--layout',
 '/gpfs/kjhan/LagForce/Codex_DPQ_backend_parity_20261017/layout.json',
 '--output',str(report),str(snapshots['serial']),str(snapshots['np4'])],check=True)
comparison=json.loads(report.read_text())
out={'scope':'Eight experimental coupled HG/RK4/HG steps; OMP2, serial vs MPI4. CPU hydro/HG, CUDA nearest neighbors.',
     'records':records,'particle_comparison':comparison}
(root/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'equal':comparison['all_physical_fields_bitwise_equal'],'time':comparison['t'],'phases':records['np4']['phase_count']}))
assert comparison['all_physical_fields_bitwise_equal']
