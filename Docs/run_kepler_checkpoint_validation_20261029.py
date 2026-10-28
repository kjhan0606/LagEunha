"""One finite auxiliary validation after an accepted t=3 checkpoint.

Runs only inside its existing Slurm allocation, after the native MPI program
has exited and both physics gates passed. Never changes the next simulation.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time


def atomic_json(path, data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(path)


def fields(text):
    return dict(token.split('=',1) for token in text.split() if '=' in token)


def numeric_steps(text, parent):
    result=[]
    for line in text.splitlines():
        value=line.split('|',1)[0].strip()
        if not value: continue
        if value in (parent+'.batch',parent+'.extern'): continue
        if not re.fullmatch(re.escape(parent)+r'\.\d+',value):
            raise ValueError('Unexpected allocation step identity: '+value)
        result.append(value)
    return result


def gate_policy(gate, radial, environment, parent):
    checks=gate.get('checks',{})
    radial_checks=radial.get('checks',{})
    if (len(checks)!=15 or not all(checks.values()) or len(radial_checks)!=9
            or not all(radial_checks.values()) or abs(gate.get('time',-1)-3.)>1e-10
            or abs(radial.get('time',-1)-3.)>1e-10):
        raise ValueError('Exact t=3 native 15+9 gates must pass first')
    if environment.get('SLURM_JOB_ID')!=parent:
        raise ValueError('This is not the intended existing allocation')
    if int(environment.get('SLURM_CPUS_PER_TASK','0'))<2:
        raise ValueError('Need the original two CPU threads')
    if int(environment.get('SLURM_MEM_PER_NODE','0'))<32768:
        raise ValueError('Need at least the original 32 GiB allocation')
    nodes=environment.get('SLURM_JOB_NUM_NODES',environment.get('SLURM_NNODES','0'))
    if nodes!='1': raise ValueError('Only the one-node allocation was qualified')


def remaining_seconds(environment, parent_info, now):
    end=environment.get('SLURM_JOB_END_TIME','')
    if end.isdecimal(): return int(end)-now
    text=parent_info.get('EndTime','')
    if not text or text in ('Unknown','None'): raise ValueError('No allocation end time')
    return datetime.datetime.fromisoformat(text).timestamp()-now


def command(*args):
    result=subprocess.run(args,text=True,capture_output=True,timeout=30)
    if result.returncode: raise RuntimeError(f'{args}: {result.stdout}{result.stderr}')
    return result.stdout


def steps(parent):
    return numeric_steps(command('squeue','--steps','--jobs',parent,'--noheader','--format=%i'),parent)


def cleanup_steps(parent):
    """Only numeric steps of this allocation, after its original MPI exit.

    No batch/extern step, whole allocation or other job can be signalled.
    State rechecks are separated by 60 seconds; there is no monitor loop.
    """
    actions=[]
    for sig in ['TERM','KILL']:
        current=steps(parent)
        if not current: return True,actions
        command('scancel','--signal='+sig,*current)
        actions.append({'signal':sig,'numeric_auxiliary_steps':current})
        time.sleep(60)
    return not steps(parent),actions


def stop_group(process):
    assert process.pid!=os.getpgrp()
    for sig,wait in [(signal.SIGTERM,30),(signal.SIGKILL,10)]:
        try: os.killpg(process.pid,sig)
        except ProcessLookupError: pass
        try: return process.wait(timeout=wait)
        except subprocess.TimeoutExpired: pass
    raise RuntimeError('Auxiliary process group did not terminate')


def has_case_output(root, cases):
    if any((root/name).exists() for name in ['comparison.json','qualifier.json','incomplete.json']):return True
    return any((root/case/name).exists() for case in cases for name in ['log','runner.out','launch_context.json'])


def run(config_path):
    config=json.loads(config_path.read_text())
    stage=Path(config['stage']).resolve();root=Path(config['benchmark']).resolve()
    parent=str(config['parent_job']);queued=str(config['queued_job'])
    if Path.cwd().resolve()!=stage or os.environ.get('SLURM_JOB_ID')!=parent:return 0
    report_path=stage/'post_checkpoint_validation.json'
    try:
        with (stage/'post_checkpoint_validation.lock').open('x') as lock:
            lock.write(str(os.getpid())+'\n')
    except FileExistsError:return 0  # Serialize duplicate calls; never replay a partial dispatch.
    if report_path.exists():return 0
    report={'started_at':datetime.datetime.now().astimezone().isoformat(),
            'parent_job':parent,'original_queued_job':queued,'benchmark':str(root),
            'original_native_evolution_finished':False,'automatic_promotion':False}
    next_before=(stage/'stage.json').read_bytes()
    claimed=False;held=False;cancelled=False;launched=False;clean=True;process=None;interrupted=False
    try:
        gate=json.loads((stage/'audit_gate.json').read_text())
        radial=json.loads((stage/'continuation_gate.json').read_text())
        gate_policy(gate,radial,os.environ,parent)
        report['original_15_plus_9_gates_passed']=True
        if json.loads(next_before).get('next_job')!=config['expected_next_job']:
            raise ValueError('Continuation changed; leave its scheduling untouched')
        parent_info=fields(command('scontrol','show','job',parent,'-o'))
        if parent_info.get('JobState')!='RUNNING':raise ValueError('Parent allocation is not running')
        available=remaining_seconds(os.environ,parent_info,time.time())
        if available<2700:raise ValueError('Need at least 45 minutes remaining in the existing allocation')
        report['parent_seconds_remaining']=available
        report['auxiliary_timeout_seconds']=min(2400,available-240)
        if steps(parent):raise ValueError('Original numeric MPI steps have not all exited')
        report['original_native_evolution_finished']=True
        if has_case_output(root,config['cases']):raise ValueError('Existing benchmark output; do not replay')
        if (root/'inline_execution_claim.json').exists():raise ValueError('Another inline execution already claimed this benchmark')
        for filename,expected in config['frozen_hashes'].items():
            if hashlib.sha256((root/filename).read_bytes()).hexdigest()!=expected:
                raise ValueError('Frozen auxiliary file changed: '+filename)
        info=fields(command('scontrol','show','job',queued,'-o'))
        report['queued_state_before']=info
        if info.get('JobState')!='PENDING' or info.get('RunTime')!='00:00:00':
            raise ValueError('Standalone is no longer unstarted and pending; leave it alone')
        if Path(info.get('WorkDir','')).resolve()!=root or info.get('JobName')!='KepVertexHaloPair':
            raise ValueError('Queued request belongs to a different workload')
        if info.get('Reason','').startswith('JobHeld'):
            raise ValueError('Do not override an existing operator hold')
        owner=re.search(r'\((\d+)\)',info.get('UserId',''))
        if not owner or int(owner[1])!=os.getuid():raise ValueError('Queued job owner mismatch')
        command('scontrol','hold',queued);held=True
        info=fields(command('scontrol','show','job',queued,'-o'))
        if info.get('JobState')!='PENDING' or info.get('RunTime')!='00:00:00':
            raise ValueError('Standalone changed state while holding it')
        claim=root/'inline_execution_claim.json'
        with claim.open('x') as stream:
            json.dump({'parent_job':parent,'queued_job':queued,'status':'claimed_before_cancel',
                       'config_sha256':hashlib.sha256(config_path.read_bytes()).hexdigest()},stream,indent=2)
        claimed=True
        command('scancel','--ctld','--state=PENDING',queued)
        info=fields(command('scontrol','show','job',queued,'-o'))
        report['queued_state_after_cancel']=info
        if info.get('JobState')!='CANCELLED' or info.get('RunTime')!='00:00:00':
            raise ValueError('Unstarted cancellation not confirmed; no inline execution')
        cancelled=True;held=False
        if has_case_output(root,config['cases']):raise ValueError('A standalone output appeared; no duplicate execution')
        verify=subprocess.run(['sha256sum','-c','audit_inputs.sha256'],cwd=root,capture_output=True,timeout=30)
        (root/'inline_preflight.stdout').write_bytes(verify.stdout)
        (root/'inline_preflight.stderr').write_bytes(verify.stderr)
        verify.check_returncode()
        available=remaining_seconds(os.environ,parent_info,time.time())
        if available<2700:raise ValueError('Allocation margin fell below 45 minutes before launch')
        report['auxiliary_timeout_seconds']=min(2400,available-240)
        report['absolute_auxiliary_deadline_epoch']=time.time()+report['auxiliary_timeout_seconds']
        report['status']='RUNNING_AUXILIARY';atomic_json(report_path,report)
        print('[POST_CHECKPOINT_VALIDATION] accepted t=3; starting frozen seven-case validation inside existing allocation',flush=True)
        with (root/'inline_batch.stdout').open('xb') as out,(root/'inline_batch.stderr').open('xb') as err:
            process=subprocess.Popen(['bash',str(root/'job_with_combinations.slurm')],cwd=root,
                                     stdout=out,stderr=err,start_new_session=True)
            launched=True;report['auxiliary_process_group']=process.pid;atomic_json(report_path,report)
            try: code=process.wait(timeout=max(0.,report['absolute_auxiliary_deadline_epoch']-time.time()))
            except subprocess.TimeoutExpired:
                report['deadline_exceeded']=True;code=stop_group(process)
        report['auxiliary_returncode']=code
        clean,actions=cleanup_steps(parent)
        report['auxiliary_steps_reaped']=clean;report['step_cleanup_actions']=actions
        comparison=root/'comparison.json'
        passed=comparison.is_file() and json.loads(comparison.read_text()).get('passed') is True
        report['auxiliary_passed']=code==0 and passed and not report.get('deadline_exceeded',False) and clean
        report['status']='AUXILIARY_PASS' if report['auxiliary_passed'] else 'AUXILIARY_FAILED_OR_INCOMPLETE'
    except BaseException as error:
        interrupted=isinstance(error,(KeyboardInterrupt,InterruptedError,SystemExit))
        report['status']='SKIPPED_OR_DISPATCH_ERROR' if not launched else 'AUXILIARY_EXCEPTION'
        report['error']=repr(error)
        if launched:
            try:
                if process.poll() is None:stop_group(process)
                clean,actions=cleanup_steps(parent)
                report['auxiliary_steps_reaped']=clean;report['step_cleanup_actions']=actions
            except BaseException as cleanup_error:
                clean=False;report['cleanup_error']=repr(cleanup_error)
    finally:
        if held and not cancelled:
            try:
                info=fields(command('scontrol','show','job',queued,'-o'))
                if info.get('JobState')=='PENDING':command('scontrol','release',queued)
            except BaseException as release_error:report['hold_release_error']=repr(release_error)
        report['queued_request_cancelled']=cancelled;report['inline_claim_created']=claimed
        report['continuation_metadata_unchanged']=(stage/'stage.json').read_bytes()==next_before
        report['finished_at']=datetime.datetime.now().astimezone().isoformat()
        atomic_json(report_path,report)
        if claimed:atomic_json(root/'inline_execution_claim.json',report)
    # A failed auxiliary test leaves the original continuation in place.
    # Unverified residual ranks require an operational stop at the saved
    # checkpoint; this is not a failed physics gate or a request to replay it.
    return 32 if (launched and not clean) or interrupted else 0


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path)
    def interrupted_signal(signum,frame):raise InterruptedError('Signal '+str(signum))
    signal.signal(signal.SIGTERM,interrupted_signal)
    signal.signal(signal.SIGINT,interrupted_signal)
    raise SystemExit(run(parser.parse_args().config))
