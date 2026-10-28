"""Check the production factor API against an existing accepted full update."""
import argparse,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');fast=base/'Codex_full_disk_fast_modules_v1_20261019'
sys.path.insert(0,str(fast))
from full_disk_stage_memory import Backend
from full_disk_factor_memory import update
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
args=p.parse_args();args.output.mkdir(exist_ok=False)
read=lambda path:dict(np.load(base/path))
a=read('Codex_full_disk_fine_first_accepted_20261019/accepted.npz')
b=read('Codex_full_disk_fine_second_midpoint_predictor_20261019/candidate/state.npz')
saved=read('Codex_full_disk_fine_second_midpoint_stage_20261019/stage0.npz')
expected=read('Codex_full_disk_fine_final_accepted_20261019/accepted.npz')
actual,factor=update(a,b,saved,.0005,.001,Backend(base/'Codex_compiled_midpoint_build_20261018/liblagmoment.so','cuda',1))
checks={key:bool(np.array_equal(actual[key],value)) for key,value in expected.items()}
report=dict(status='passed_factor_memory_against_accepted' if all(checks.values()) else 'failed_factor_memory_against_accepted',
    exact_array_values=checks,factor=factor,provider_calls=0,evolution_steps=0,
    limits='Same saved h=.0005 full update. No new RHS or accepted evolution; binary80 padding bytes are not compared.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2));assert all(checks.values())
