"""Validate parallel/vector path samples against a preserved accepted step."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic
from predictive_evolution_batch import initialize_path_worker,assess


def run():
    saved=json.loads(Path('reference_step.json').read_text());source=json.loads(Path('inputs/n4_r100/checkpoint.json').read_text())
    record=saved['record'];mass=np.asarray(source['mass']);cfg=source['config'];samples={};started=time.monotonic()
    with ProcessPoolExecutor(max_workers=4,mp_context=mp.get_context('spawn'),initializer=initialize_path_worker,
        initargs=(cfg['n'],cfg['amplitude'],True)) as executor:
        result=assess(record,mass,1e-5,None,65,samples,Path('parallel_nodes.json'),'validation',executor)
    flux=max(float(np.max(abs(np.asarray(samples[k]['mass_rate'])-v['mass_rate']))) for k,v in saved['samples'].items())
    difference=max(float(np.max(abs(np.asarray(result[k])-np.asarray(saved['audit'][k])))) for k in ('rate','generator','pressure','factor_R','signed','absolute','quad_delta'))
    report=dict(passed=flux<1e-11 and difference<1e-11,max_flux_difference=flux,max_audit_difference=difference,seconds=time.monotonic()-started,workers=4,audit=result)
    atomic('parallel_report.json',report);print(report,flush=True)
    raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
