"""One saved-data postprocessing bundle after accepted Ada coarse/fine .05.

No polling, RHS, reconstruction fits, time update, retry or extension launch.
Missing native/restart evidence remains pending for the final evaluation.
"""
import argparse,concurrent.futures,hashlib,json,os,subprocess,sys
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--existing-platform-comparison',type=Path)
    p.add_argument('--existing-platform-budget',type=Path)
    args=p.parse_args();assert not args.output.exists()
    assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='1'
    b=Path('/gpfs/kjhan/LagForce');root=b/'Codex_full_disk_continuous_two_dt_to_t1_20261019'
    checkpoints={name:root/name/'checkpoints'/f'{name}_{step:04d}' for name,step in [('coarse',10),('fine',20)]}
    input_evidence={}
    for branch,cp in checkpoints.items():
        result=json.loads((root/branch/'result.json').read_text())
        latest=json.loads((root/branch/f'latest_{branch}.json').read_text())
        gate=json.loads((cp/'physical_gate.json').read_text())
        assert result['status']=='completed_comparison_branch' and result['time']==.05
        assert latest['time']==.05 and latest['step']==(10 if branch=='coarse' else 20)
        assert gate['status']=='passed_checkpoint_physical_monitor' and gate['time']==.05
        for f in (root/branch/'result.json',root/branch/f'latest_{branch}.json',cp/'physical_gate.json',cp/'accepted.npz'):
            input_evidence[str(f)]=sha(f)
    assert bool(args.existing_platform_comparison)==bool(args.existing_platform_budget)
    reused_platform=None
    if args.existing_platform_comparison:
        from compare_full_moment_checkpoint_errors import load
        comparison=json.loads(args.existing_platform_comparison.read_text())
        platform=json.loads(args.existing_platform_budget.read_text())
        native=b/'Codex_full_disk_native_samehost_original_and_coarse_20261019/mpi_coarse'
        accepted=json.loads((native/'accepted_0010.json').read_text())
        assert accepted['accepted_commit'] and accepted['time']==.05
        assert comparison['status']=='completed_checkpoint_actual_cell_error_difference_diagnostics'
        assert comparison['time']==platform['time']==.05
        assert platform['status']=='passed_platform_physical_difference_budget'
        assert platform['kind']=='platform' and platform['regions']
        assert all(row['passed'] for row in platform['regions'].values())
        assert platform['comparison_sha256']==sha(args.existing_platform_comparison)
        assert Path(platform['comparison_path']).resolve()==args.existing_platform_comparison.resolve()
        for key,folder in [('source_A',checkpoints['coarse']/'physics'),('source_B',native/'candidate_audit_0010/physics')]:
            _,_,current=load(folder)
            assert current==comparison[key], 'Existing platform comparison source changed'
        reused_platform=dict(comparison=comparison,budget=platform)
        for f in (args.existing_platform_comparison,args.existing_platform_budget):
            input_evidence[str(f)]=sha(f)
    args.output.mkdir()
    modules=Path(__file__).resolve().parent;tasks=[]
    def cmd(name,*arguments):return [sys.executable,str(modules/name),*map(str,arguments)]
    cavity=b/'Codex_full_disk_checkpoint_error_diagnostics_20261019/saved_density_volume_drift/initial_density_volume.npz'
    tasks.append(('density_temporal_budget',cmd('evaluate_full_moment_saved_checkpoint_bundle.py','--root',root,'--initial-cavity',cavity,'--output',args.output/'density_temporal_budget')))
    tasks.append(('fine_finite_pressure_budget',cmd('diagnose_full_moment_finite_pressure_budget.py','--a',root/'fine/checkpoints/fine_0010/physics','--b',checkpoints['fine']/'physics','--output',args.output/'fine_finite_pressure_budget')))
    for branch,cp in checkpoints.items():
        tasks.append((branch+'_subcell_energy',cmd('diagnose_full_moment_subcell_energy.py','--checkpoint',cp/'physics','--output',args.output/(branch+'_subcell_energy'))))
    tasks.append(('density_image',cmd('plot_full_moment_two_dt_density.py','--a',checkpoints['coarse']/'physics','--b',checkpoints['fine']/'physics','--geometry-a',checkpoints['coarse']/'geometry','--geometry-b',checkpoints['fine']/'geometry','--output',args.output/'density_image')))
    def run(task):
        label,command=task
        with (args.output/(label+'.log')).open('x') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        return label
    manifest=dict(status='started_saved_comparison_postprocessing',input_sha256=input_evidence,
        source_sha256={str(modules/name):sha(modules/name) for name in (Path(__file__).name,'compare_full_moment_checkpoint_errors.py','evaluate_full_moment_saved_checkpoint_bundle.py','diagnose_full_moment_finite_pressure_budget.py','diagnose_full_moment_subcell_energy.py','plot_full_moment_two_dt_density.py','evaluate_full_moment_checkpoint_difference_budget.py')},
        subprocesses_maximum=3,BLAS_threads=1,OMP_threads=1,physical_providers=0,accepted_time_steps=0,automatic_extension_authorized=False)
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:completed=list(pool.map(run,tasks))
    assert all(sha(f)==digest for f,digest in input_evidence.items())
    budget=json.loads((args.output/'density_temporal_budget/report.json').read_text())
    assert budget['all_three_times_available'] and not budget['pending']
    native=b/'Codex_full_disk_native_samehost_original_and_coarse_20261019/mpi_coarse'
    pending=[]
    if reused_platform is not None:
        (args.output/'reused_platform_evidence.json').write_text(json.dumps(reused_platform,indent=2)+'\n')
        completed+=['reused_same_h_native_ada_t0p05','reused_same_h_native_ada_t0p05_platform_budget']
    elif (native/'accepted_0010.json').exists():
        accepted=json.loads((native/'accepted_0010.json').read_text());assert accepted['accepted_commit'] and accepted['time']==.05
        run(('same_h_native_ada_t0p05',cmd('compare_full_moment_checkpoint_errors.py','--a',checkpoints['coarse']/'physics','--b',native/'candidate_audit_0010/physics','--output',args.output/'same_h_native_ada_t0p05')))
        run(('same_h_native_ada_t0p05_platform_budget',cmd('evaluate_full_moment_checkpoint_difference_budget.py','--comparison',args.output/'same_h_native_ada_t0p05/report.json','--kind','platform','--output',args.output/'same_h_native_ada_t0p05_platform_budget')))
        completed+=['same_h_native_ada_t0p05','same_h_native_ada_t0p05_platform_budget']
    else:pending.append('Native accepted t=.05 and same-h platform comparison')
    pending.append('Actual native restart equivalence: inspect Slurm427904 terminal comparison report separately')
    assert all(sha(f)==digest for f,digest in input_evidence.items())
    manifest.update(status='completed_saved_comparison_postprocessing',completed=completed,pending=pending,
        available_density_passed=budget['available_density_passed'],available_temporal_passed=budget['available_temporal_passed'],
        limits='Postprocessing completion is not a physics pass or extension approval. Existing accepted states untouched; contact-growth and missing evidence require the completed-bundle evaluation.')
    (args.output/'report.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ('status','completed','pending','available_density_passed','available_temporal_passed','automatic_extension_authorized')}),flush=True)


if __name__=='__main__':main()
