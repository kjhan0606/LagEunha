"""One-factor-at-a-time audit of the failed N3 wave path."""
import json,os,sys
from pathlib import Path
from instantaneous_weight_path import run
CASES={'baseline':dict(tol=1e-10,max_step=1/16,root_xtol=1e-9),
       'root_tighter':dict(tol=1e-10,max_step=1/16,root_xtol=1e-11),
       'step_smaller':dict(tol=1e-10,max_step=1/64,root_xtol=1e-9),
       'time_tighter':dict(tol=1e-12,max_step=1/16,root_xtol=1e-9)}
if __name__=='__main__':
    root=Path.cwd()
    for name in sys.argv[1:]:
        folder=root/name;folder.mkdir();os.chdir(folder)
        try:r=run('wave',.05,n=3,analytic=True,**CASES[name])
        except Exception as e:r=dict(status='failed',error=repr(e),**CASES[name])
        r['case']=name
        with (root/'results.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
        print(json.dumps(r),flush=True);os.chdir(root)
