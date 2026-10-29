"""Fault-injection checks: no premature restart and no accidental full TSV."""
import json
import os
from pathlib import Path
import tempfile
from unittest.mock import patch
import hashlib
from run_kepler_q1_restart_dt_bundle_20261029 import parent_checkpoint, validate_environment


def rejected(call, substring):
    try:call()
    except ValueError as exc:
        assert substring in str(exc),str(exc)
    else:raise AssertionError('Invalid input was accepted')


def main():
    expected={'target':10.1,'cache':0,'cap':None}
    good={'GFS_KEPLER_WALL_EVENTS':'1','GFS_KEPLER_FIXED_POLYGON':'1',
          'GFS_GPU_RAW_FACES':'1','GFS_GPU_FINAL_VOLUME':'1','GFS_GPU_TESS_DUMP':'0',
          'GFS_FAIL_NONPOSITIVE_U':'1','SEDOV_PHASE1':'1','EUNHA_IC':'kepler',
          'OMP_NUM_THREADS':'8','GFS_EXACT_STOP':'1','GFS_RK4_GAP_PROFILE':'1',
          'GFS_NATIVE_CONST_CONFIG':'0','HYDRO_TSTOP':'10.1','EUNHA_DUMP_DT':'1000000000'}
    with patch.dict(os.environ,good,clear=True):validate_environment(expected)
    for key,value in [('GFS_RK4_AUDIT_STEP_MIN','0'),('GFS_RK4_AUDIT_IDS','7'),
                      ('GFS_GPU_FACE_DUMP_FIRST','0'),('GFS_THERMAL64','1'),('GFS_DYN_W2','1')]:
        with patch.dict(os.environ,{**good,key:value},clear=True):
            rejected(lambda:validate_environment(expected),'Unbounded diagnostics')
    with patch.dict(os.environ,{**good,'GFS_DIAG_DT_MAX':'0.1'},clear=True):
        rejected(lambda:validate_environment(expected),'Unrequested timestep')
    with patch.dict(os.environ,{**good,'HYDRO_TSTOP':'0.1'},clear=True):
        rejected(lambda:validate_environment(expected),'Unexpected output target')
    with patch.dict(os.environ,{**good,'GFS_NATIVE_CONST_CONFIG':'1'},clear=True):
        rejected(lambda:validate_environment(expected),'Unexpected effective environment')
    # Real green t3 data must never substitute for an accepted t10 parent.
    actual=Path('/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/chain_h200_v2')
    fixture=Path(tempfile.mkdtemp(prefix='kepler-dt-parent-guard-',dir='/tmp'))
    (fixture/'chain_state.json').write_text(json.dumps({'status':'running'}))
    rejected(lambda:parent_checkpoint({'parent_root':str(fixture)}),'not complete')
    (fixture/'chain_state.json').write_text(json.dumps({'status':'finite_chain_complete','job':'439333'}))
    (fixture/'to10p0').symlink_to(actual/'to3p0')
    binary=actual/'to3p0/eunha2.exe'
    cfg={'parent_root':str(fixture),'parent_job':'439333','parent_params_sha256':hashlib.sha256((actual/'to3p0/params.dat').read_bytes()).hexdigest(),
         'baseline':str(binary),'binary_hashes':{str(binary):hashlib.sha256(binary.read_bytes()).hexdigest()}}
    rejected(lambda:parent_checkpoint(cfg),'Exact t10')
    print('PASS valid bounded environment; reject full audit, face dumps, incompatible physics, wrong cap/target/cache; reject incomplete parent and actual green t3 as t10')


if __name__=='__main__':main()
