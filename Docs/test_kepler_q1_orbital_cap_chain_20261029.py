"""Exercise real checkpoint and rendered environment guards without MPI."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from unittest.mock import patch
import run_kepler_q1_orbital_cap_chain_20261029 as chain


def main():
    preflight=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/gravity_cap_preflight_v1')
    contract=json.loads((preflight/'contract.json').read_text())
    old=json.loads((preflight/'restart_eligibility.json').read_text())
    geometry=json.loads((preflight/'geometry_t0825.json').read_text())
    actual=chain.core_screen(Path(contract['initial']),Path(contract['earliest_chain_restart_candidate']),geometry,contract)
    assert not actual['passed']
    for a,b in zip(actual['regions'],old['regions']):assert a['checks']==b['checks']
    staged=Path(sys.argv[1]).resolve()
    staged_cfg=json.loads((staged/'chain.json').read_text())
    starting_step=str(staged_cfg['starting_step'])
    tmp=Path(tempfile.mkdtemp(prefix='cap-chain-test-',dir='/tmp'))
    for name in ['chain.json','params.dat','native_template.sh','contract.json']:
        shutil.copy2(staged/name,tmp/name)
    original_run=subprocess.run;observed=[]
    def intercepted(args,**kwargs):
        if args[0]=='bash' and len(args)==3 and args[-1]==starting_step:
            script=Path(args[1]);s=script.read_text()
            assert s.count('set +e\n')==1
            # Stop before the real MPI command even if a guard is broken.
            dry=script.with_name('dry_guard.sh');dry.write_text(s.replace('set +e\n','exit 0\nset +e\n'))
            case=script.with_name('case.json');cfg=json.loads(case.read_text())
            bad={**cfg,'target':99.};case.write_text(json.dumps(bad))
            rejected=original_run(['bash',str(dry)],cwd=script.parent,capture_output=True,text=True)
            case.write_text(json.dumps(cfg))
            accepted=original_run(['bash',str(dry)],cwd=script.parent,capture_output=True,text=True)
            assert rejected.returncode!=0 and accepted.returncode==0,(rejected.stderr,accepted.stderr)
            observed.append({'wrong_target_rejected':True,'valid_rendered_environment_accepted':True})
            return subprocess.CompletedProcess(args,99)
        return original_run(args,**kwargs)
    old_cwd=Path.cwd()
    try:
        os.chdir(tmp)
        env=dict(os.environ)
        env.update({'GFS_RK4_AUDIT_STEP_MIN':'0','GFS_RK4_AUDIT_STEP_MAX':'1000000'})
        if staged_cfg['execution_mode']=='ada':
            env.pop('SLURM_JOB_ID',None);env.pop('SLURM_CPUS_PER_TASK',None)
            env.update({'OMP_NUM_THREADS':'8','CUDA_VISIBLE_DEVICES':staged_cfg['ada_gpu_uuid']})
        else:env.update({'SLURM_JOB_ID':'TEST_NO_MPI','SLURM_CPUS_PER_TASK':'8'})
        with patch.dict(os.environ,env,clear=True), \
             patch.object(chain.socket,'gethostname',return_value='LagEunha'), \
             patch.object(chain.subprocess,'run',side_effect=intercepted):
            try:chain.run(tmp)
            except RuntimeError as exc:assert 'Native run failed' in str(exc)
            else:raise AssertionError('Fake execution unexpectedly succeeded')
    finally:os.chdir(old_cwd)
    assert observed and json.loads((tmp/'chain_state.json').read_text())['status']=='stopped'
    print('PASS real t=.825 core failure, actual rendered bounded environment, wrong-target rejection and no-MPI stop')


if __name__=='__main__':main()
