"""Isolate max-step effect for the previously failed shear case."""
import json
from instantaneous_weight_path import run
if __name__=='__main__':
    for step in (1/32,1/64):
        try:r=run('shear',.05,1e-8,n=2,analytic=False,max_step=step)
        except Exception as e:r=dict(kind='shear',duration=.05,rtol=1e-8,max_step=step,status='failed',error=repr(e))
        print(json.dumps(r),flush=True)
